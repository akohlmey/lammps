// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_rigid_small_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain.h"
#include "domain_kokkos.h"
#include "error.h"
#include "kokkos_few.h"
#include "math_extra_kokkos.h"
#include "memory_kokkos.h"
#include "rigid_const.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace RigidConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidSmallKokkos<DeviceType>::FixRigidSmallKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigidSmall(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  commKK = (CommKokkos *) comm;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | F_MASK | TAG_MASK | TYPE_MASK | MASK_MASK |
                  IMAGE_MASK | RMASS_MASK;
  datamask_modify = X_MASK | V_MASK;

  // size of the per-atom exchange buffer (used once device exchange is enabled)

  maxexchange = 12 + bodysize;

  // The FixRigidSmall constructor already allocated the per-atom arrays with
  // the plain Memory class and populated bodytag (in create_bodies) and bodyown
  // for the owned atoms. Convert the arrays to DualViews (so the base-class host
  // pointers alias the host side of each view) while preserving bodytag/bodyown;
  // atom2body/xcmimage/displace are (re)computed later in setup_bodies_static().

  tagint *bodytag_tmp = bodytag;
  int *bodyown_tmp = bodyown;
  bodytag = nullptr;
  bodyown = nullptr;
  memory->destroy(atom2body);  atom2body = nullptr;
  memory->destroy(xcmimage);   xcmimage = nullptr;
  memory->destroy(displace);   displace = nullptr;

  grow_arrays(atom->nmax);

  for (int i = 0; i < atom->nlocal; i++) {
    bodytag[i] = bodytag_tmp[i];
    bodyown[i] = bodyown_tmp[i];
  }
  k_bodytag.modify_host();
  k_bodyown.modify_host();
  memory->destroy(bodytag_tmp);
  memory->destroy(bodyown_tmp);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidSmallKokkos<DeviceType>::~FixRigidSmallKokkos()
{
  if (copymode) return;

  // null the base-class pointers so ~FixRigidSmall does not double-free the
  // Kokkos-owned allocations

  memoryKK->destroy_kokkos(k_bodyown, bodyown);
  memoryKK->destroy_kokkos(k_bodytag, bodytag);
  memoryKK->destroy_kokkos(k_atom2body, atom2body);
  memoryKK->destroy_kokkos(k_xcmimage, xcmimage);
  memoryKK->destroy_kokkos(k_displace, displace);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::setmask()
{
  // need a pre_exchange hook to bring the device-resident body data back to
  // the host before atoms (and their bodies) migrate during exchange

  return FixRigidSmall::setmask() | PRE_EXCHANGE;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::init()
{
  FixRigidSmall::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use run_style respa with fix {}", style);
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays as DualViews
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::grow_arrays(int nmax)
{
  // The FixRigidSmall (host) base class writes these per-atom arrays directly
  // through the aliased host pointers (e.g. in create_bodies) without marking
  // the DualView host side modified. grow_kokkos() resizes the device view and
  // rebuilds the host mirror from it, so unless the host data is first pushed
  // to the device it would be lost. Mark host modified and sync to device so
  // the resize preserves the current (host) data.

  k_bodyown.modify_host();   k_bodyown.sync_device();
  k_bodytag.modify_host();   k_bodytag.sync_device();
  k_atom2body.modify_host(); k_atom2body.sync_device();
  k_xcmimage.modify_host();  k_xcmimage.sync_device();
  k_displace.modify_host();  k_displace.sync_device();

  memoryKK->grow_kokkos(k_bodyown,bodyown,nmax,"rigid/small:bodyown");
  memoryKK->grow_kokkos(k_bodytag,bodytag,nmax,"rigid/small:bodytag");
  memoryKK->grow_kokkos(k_atom2body,atom2body,nmax,"rigid/small:atom2body");
  memoryKK->grow_kokkos(k_xcmimage,xcmimage,nmax,"rigid/small:xcmimage");
  memoryKK->grow_kokkos(k_displace,displace,nmax,3,"rigid/small:displace");

  d_bodyown = k_bodyown.view<DeviceType>();
  d_bodytag = k_bodytag.view<DeviceType>();
  d_atom2body = k_atom2body.view<DeviceType>();
  d_xcmimage = k_xcmimage.view<DeviceType>();
  d_displace = k_displace.view<DeviceType>();

  // extended-particle arrays remain host-only for now

  if (extended) {
    memory->grow(eflags,nmax,"rigid/small:eflags");
    if (orientflag) memory->grow(orient,nmax,orientflag,"rigid/small:orient");
    if (dorientflag) memory->grow(dorient,nmax,3,"rigid/small:dorient");
  }

  // per-atom virial, grown the same way as in the base class

  if (nmax > maxvatom) {
    maxvatom = atom->nmax;
    memory->grow(vatom,maxvatom,6,"fix:vatom");
  }
}

/* ----------------------------------------------------------------------
   lifecycle hooks
   For this build stage all per-body integration still runs in the
   FixRigidSmall (host) base class. We only ensure the atom data the base
   methods read/write is synced to/from the host around each call so the
   fix is correct under a Kokkos run and coexists with other Kokkos styles.
   Subsequent build stages replace these host calls with device kernels.
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::setup(int vflag)
{
  atomKK->sync(Host, datamask_read);
  FixRigidSmall::setup(vflag);
  atomKK->modified(Host, datamask_modify);

  // the host work above modified per-atom data on the host, but ModifyKokkos
  // brackets this callback with modified(Device, datamask_modify) (the fix
  // execution space is the device). push the host changes to the device now so
  // the host side is not left "ahead", which would trip the DualView
  // concurrent host/device modification check on a GPU build.

  atomKK->sync(execution_space, datamask_modify);

  // body data and the per-atom body arrays become resident on the device for
  // the run loop (extended particles stay on the host)

  if (!extended) {
    copy_body_to_device();
    sync_peratom_to_device();
    body_resident_device = true;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::setup_pre_neighbor()
{
  atomKK->sync(Host, datamask_read);
  FixRigidSmall::setup_pre_neighbor();
  atomKK->modified(Host, datamask_modify);
  atomKK->sync(execution_space, datamask_modify);
}

/* ----------------------------------------------------------------------
   bring the device-resident body data back to the host before atoms (and
   their bodies) migrate during the (host) exchange
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::pre_exchange()
{
  if (body_resident_device) {
    copy_body_to_host();
    body_resident_device = false;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::pre_neighbor()
{
  atomKK->sync(Host, datamask_read);
  FixRigidSmall::pre_neighbor();
  atomKK->modified(Host, datamask_modify);
  atomKK->sync(execution_space, datamask_modify);

  // the host rebuilt the (ghost) body list and the per-atom body arrays;
  // make them resident on the device again

  if (!extended) {
    copy_body_to_device();
    sync_peratom_to_device();
    body_resident_device = true;
  }
}

/* ----------------------------------------------------------------------
   push the per-atom body arrays (written on the host in setup/pre_neighbor)
   to the device. They are only read by the device kernels, so this is done
   once per reneighbor rather than every step.
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::sync_peratom_to_device()
{
  k_atom2body.modify_host(); k_atom2body.sync_device();
  k_displace.modify_host();  k_displace.sync_device();
  k_xcmimage.modify_host();  k_xcmimage.sync_device();
  k_bodyown.modify_host();   k_bodyown.sync_device();
  d_atom2body = k_atom2body.view<DeviceType>();
  d_displace = k_displace.view<DeviceType>();
  d_xcmimage = k_xcmimage.view<DeviceType>();
  d_bodyown = k_bodyown.view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::initial_integrate(int vflag)
{
  // extended particles are not yet handled on the device: run the full base
  // (host) integration in that case

  if (extended) {
    atomKK->sync(Host, datamask_read);
    FixRigidSmall::initial_integrate(vflag);
    atomKK->modified(Host, datamask_modify);
    atomKK->sync(execution_space, datamask_modify);
    return;
  }

  // update body COM position/orientation on the device (body data resident)

  auto l_body = d_body;
  const double l_dtf = dtf, l_dtv = dtv, l_dtq = dtq;

  Kokkos::parallel_for("rigid/small/kk initial_integrate",
    Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
    KOKKOS_LAMBDA(const int ibody) {
      Body &b = l_body(ibody);

      const double dtfm = l_dtf / b.mass;
      b.vcm[0] += dtfm * b.fcm[0];
      b.vcm[1] += dtfm * b.fcm[1];
      b.vcm[2] += dtfm * b.fcm[2];

      b.xcm[0] += l_dtv * b.vcm[0];
      b.xcm[1] += l_dtv * b.vcm[1];
      b.xcm[2] += l_dtv * b.vcm[2];

      b.angmom[0] += l_dtf * b.torque[0];
      b.angmom[1] += l_dtf * b.torque[1];
      b.angmom[2] += l_dtf * b.torque[2];

      MathExtraKokkos::angmom_to_omega(b.angmom,b.ex_space,b.ey_space,
                                       b.ez_space,b.inertia,b.omega);
      MathExtraKokkos::richardson(b.quat,b.angmom,b.omega,b.inertia,l_dtq);
      MathExtraKokkos::q_to_exyz(b.quat,b.ex_space,b.ey_space,b.ez_space);
    });

  // virial setup before call to set_xv

  v_init(vflag);

  // forward communicate updated body info to ghost bodies (device)

  commflag = INITIAL;
  commKK->template forward_comm_device<DeviceType>(this,29);

  // set coords/orient and velocity/rotation of atoms in rigid bodies (device)

  set_xv_kokkos(1);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::post_force(int /*vflag*/)
{
  // langextra is computed on the host (RNG) from the current body vcm/omega,
  // so the device-resident body data must be synced to the host first; it is
  // added to fcm/torque (back on the device) in compute_forces_and_torques_kokkos()

  if (langflag) {
    if (body_resident_device) copy_body_to_host();
    apply_langevin_thermostat();
  }
  if (earlyflag) compute_forces_and_torques_kokkos();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::final_integrate()
{
  if (extended) {
    atomKK->sync(Host, datamask_read);
    FixRigidSmall::final_integrate();
    atomKK->modified(Host, datamask_modify);
    atomKK->sync(execution_space, datamask_modify);
    return;
  }

  if (!earlyflag) compute_forces_and_torques_kokkos();

  // update body vcm/angmom/omega on the device (body data resident).
  // for 2d, enforce2d() zeroing is folded in before the update.

  auto l_body = d_body;
  const double l_dtf = dtf;
  const int dim2 = (domain->dimension == 2);

  Kokkos::parallel_for("rigid/small/kk final_integrate",
    Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
    KOKKOS_LAMBDA(const int ibody) {
      Body &b = l_body(ibody);

      if (dim2) {
        b.xcm[2] = 0.0; b.vcm[2] = 0.0; b.fcm[2] = 0.0; b.xgc[2] = 0.0;
        b.torque[0] = 0.0; b.torque[1] = 0.0;
        b.angmom[0] = 0.0; b.angmom[1] = 0.0;
        b.omega[0] = 0.0; b.omega[1] = 0.0;
      }

      const double dtfm = l_dtf / b.mass;
      b.vcm[0] += dtfm * b.fcm[0];
      b.vcm[1] += dtfm * b.fcm[1];
      b.vcm[2] += dtfm * b.fcm[2];

      b.angmom[0] += l_dtf * b.torque[0];
      b.angmom[1] += l_dtf * b.torque[1];
      b.angmom[2] += l_dtf * b.torque[2];

      MathExtraKokkos::angmom_to_omega(b.angmom,b.ex_space,b.ey_space,
                                       b.ez_space,b.inertia,b.omega);
    });

  // forward communicate updated body info to ghost bodies (device)

  commflag = FINAL;
  commKK->template forward_comm_device<DeviceType>(this,10);

  // set velocity/rotation of atoms in rigid bodies (device)
  // virial is already setup from initial_integrate

  set_xv_kokkos(0);
}

/* ----------------------------------------------------------------------
   sum the per-atom forces/torques into the rigid bodies on the device.
   the ghost->owner reduction (reverse_comm) and the Langevin/gravity body
   contributions remain on the host.
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::compute_forces_and_torques_kokkos()
{
  // body data and the per-atom body arrays are resident on the device

  atomKK->sync(execution_space, X_MASK | F_MASK);

  auto l_x = atomKK->k_x.view<DeviceType>();
  auto l_f = atomKK->k_f.view<DeviceType>();
  auto l_body = d_body;
  auto l_atom2body = d_atom2body;
  auto l_xcmimage = d_xcmimage;

  const int nbody = nlocal_body + nghost_body;
  const int nlocal = atom->nlocal;

  Few<double,3> prd(domain->prd);
  Few<double,6> h(domain->h);
  const int triclinic = domain->triclinic;

  // zero forces/torques on all bodies (owned + ghost)

  Kokkos::parallel_for("rigid/small/kk zero fcm/torque",
    Kokkos::RangePolicy<DeviceType>(0, nbody),
    KOKKOS_LAMBDA(const int ibody) {
      Body &b = l_body(ibody);
      b.fcm[0] = b.fcm[1] = b.fcm[2] = 0.0;
      b.torque[0] = b.torque[1] = b.torque[2] = 0.0;
    });

  // accumulate atom forces/torques into the (owned or ghost) body

  Kokkos::parallel_for("rigid/small/kk sum fcm/torque",
    Kokkos::RangePolicy<DeviceType>(0, nlocal),
    KOKKOS_LAMBDA(const int i) {
      if (l_atom2body(i) < 0) return;
      Body &b = l_body(l_atom2body(i));

      Kokkos::atomic_add(&b.fcm[0], l_f(i,0));
      Kokkos::atomic_add(&b.fcm[1], l_f(i,1));
      Kokkos::atomic_add(&b.fcm[2], l_f(i,2));

      Few<double,3> xi;
      xi[0] = l_x(i,0); xi[1] = l_x(i,1); xi[2] = l_x(i,2);
      Few<double,3> unwrap = DomainKokkos::unmap(prd,h,triclinic,xi,l_xcmimage(i));
      const double dx = unwrap[0] - b.xcm[0];
      const double dy = unwrap[1] - b.xcm[1];
      const double dz = unwrap[2] - b.xcm[2];

      Kokkos::atomic_add(&b.torque[0], dy*l_f(i,2) - dz*l_f(i,1));
      Kokkos::atomic_add(&b.torque[1], dz*l_f(i,0) - dx*l_f(i,2));
      Kokkos::atomic_add(&b.torque[2], dx*l_f(i,1) - dy*l_f(i,0));
    });

  // reverse communicate ghost body fcm/torque into the owners (device)

  commflag = FORCE_TORQUE;
  commKK->template reverse_comm_device<DeviceType>(this,6);

  // add gravity force to the COM of each owned body (device)

  if (id_gravity) {
    auto l_body2 = d_body;
    const double gx = gvec[0], gy = gvec[1], gz = gvec[2];
    Kokkos::parallel_for("rigid/small/kk gravity",
      Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
      KOKKOS_LAMBDA(const int ibody) {
        Body &b = l_body2(ibody);
        b.fcm[0] += gx*b.mass;
        b.fcm[1] += gy*b.mass;
        b.fcm[2] += gz*b.mass;
      });
  }

  // include Langevin thermostat forces/torques (host fallback for the RNG)

  if (langflag) {
    copy_body_to_host();
    for (int ibody = 0; ibody < nlocal_body; ibody++) {
      body[ibody].fcm[0] += langextra[ibody][0];
      body[ibody].fcm[1] += langextra[ibody][1];
      body[ibody].fcm[2] += langextra[ibody][2];
      body[ibody].torque[0] += langextra[ibody][3];
      body[ibody].torque[1] += langextra[ibody][4];
      body[ibody].torque[2] += langextra[ibody][5];
    }
    copy_body_to_device();
  }
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   on the device. setx = 1 also updates positions (initial_integrate),
   setx = 0 updates velocities only (final_integrate).
   (per-atom virial and extended particles are handled on the host)
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::set_xv_kokkos(int setx)
{
  // per-atom (and centroid) virial accumulation is not yet on the device;
  // fall back to the host base routine when it is requested

  if (evflag && (vflag_atom || cvflag_atom)) {
    atomKK->sync(Host, X_MASK | V_MASK | F_MASK | TYPE_MASK | RMASS_MASK);
    if (setx) FixRigidSmall::set_xv();
    else FixRigidSmall::set_v();
    atomKK->modified(Host, X_MASK | V_MASK);
    atomKK->sync(execution_space, X_MASK | V_MASK);
    return;
  }

  const double xprd = domain->xprd;
  const double yprd = domain->yprd;
  const double zprd = domain->zprd;
  const double xy = domain->xy;
  const double xz = domain->xz;
  const double yz = domain->yz;
  const int triclinic = domain->triclinic;
  const int dim = domain->dimension;
  const int l_evflag = evflag;
  const double l_dtf = dtf;
  const int rmass_flag = (atom->rmass != nullptr);

  // body data and per-atom body arrays are resident on the device

  atomKK->sync(execution_space, X_MASK | V_MASK | F_MASK | TYPE_MASK | RMASS_MASK);
  atomKK->k_mass.template sync<DeviceType>();
  auto l_x = atomKK->k_x.view<DeviceType>();
  auto l_v = atomKK->k_v.view<DeviceType>();
  auto l_f = atomKK->k_f.view<DeviceType>();
  auto l_mass = atomKK->k_mass.view<DeviceType>();
  auto l_rmass = atomKK->k_rmass.view<DeviceType>();
  auto l_type = atomKK->k_type.view<DeviceType>();

  auto l_body = d_body;
  auto l_atom2body = d_atom2body;
  auto l_displace = d_displace;
  auto l_xcmimage = d_xcmimage;
  const int nlocal = atom->nlocal;

  EV_FLOAT ev;
  Kokkos::parallel_reduce("rigid/small/kk set_xv",
    Kokkos::RangePolicy<DeviceType>(0, nlocal),
    KOKKOS_LAMBDA(const int i, EV_FLOAT &ev) {
      if (l_atom2body(i) < 0) return;
      const Body &b = l_body(l_atom2body(i));

      const int xbox = (l_xcmimage(i) & IMGMASK) - IMGMAX;
      const int ybox = (l_xcmimage(i) >> IMGBITS & IMGMASK) - IMGMAX;
      const int zbox = (l_xcmimage(i) >> IMG2BITS) - IMGMAX;

      // old velocity / position (for the constraint virial)

      const double v0 = l_v(i,0), v1 = l_v(i,1), v2 = l_v(i,2);
      const double ox0 = l_x(i,0), ox1 = l_x(i,1), ox2 = l_x(i,2);

      double delta[3];
      delta[0] = b.ex_space[0]*l_displace(i,0) + b.ey_space[0]*l_displace(i,1) + b.ez_space[0]*l_displace(i,2);
      delta[1] = b.ex_space[1]*l_displace(i,0) + b.ey_space[1]*l_displace(i,1) + b.ez_space[1]*l_displace(i,2);
      delta[2] = b.ex_space[2]*l_displace(i,0) + b.ey_space[2]*l_displace(i,1) + b.ez_space[2]*l_displace(i,2);

      double nv0 = b.omega[1]*delta[2] - b.omega[2]*delta[1] + b.vcm[0];
      double nv1 = b.omega[2]*delta[0] - b.omega[0]*delta[2] + b.vcm[1];
      double nv2 = b.omega[0]*delta[1] - b.omega[1]*delta[0] + b.vcm[2];

      if (dim == 2) { nv2 = 0.0; delta[2] = 0.0; }

      l_v(i,0) = nv0; l_v(i,1) = nv1; l_v(i,2) = nv2;

      if (setx) {
        if (triclinic == 0) {
          l_x(i,0) = delta[0] + b.xcm[0] - xbox*xprd;
          l_x(i,1) = delta[1] + b.xcm[1] - ybox*yprd;
          l_x(i,2) = delta[2] + b.xcm[2] - zbox*zprd;
        } else {
          l_x(i,0) = delta[0] + b.xcm[0] - xbox*xprd - ybox*xy - zbox*xz;
          l_x(i,1) = delta[1] + b.xcm[1] - ybox*yprd - zbox*yz;
          l_x(i,2) = delta[2] + b.xcm[2] - zbox*zprd;
        }
      }

      if (l_evflag) {
        const double massone = rmass_flag ? l_rmass(i) : l_mass(l_type(i));
        const double fc0 = massone*(nv0 - v0)/l_dtf - l_f(i,0);
        const double fc1 = massone*(nv1 - v1)/l_dtf - l_f(i,1);
        const double fc2 = massone*(nv2 - v2)/l_dtf - l_f(i,2);

        double X0, X1, X2;
        if (triclinic == 0) {
          X0 = ox0 + xbox*xprd; X1 = ox1 + ybox*yprd; X2 = ox2 + zbox*zprd;
        } else {
          X0 = ox0 + xbox*xprd + ybox*xy + zbox*xz;
          X1 = ox1 + ybox*yprd + zbox*yz;
          X2 = ox2 + zbox*zprd;
        }

        ev.v[0] += 0.5*X0*fc0;
        ev.v[1] += 0.5*X1*fc1;
        ev.v[2] += 0.5*X2*fc2;
        ev.v[3] += 0.5*X0*fc1;
        ev.v[4] += 0.5*X0*fc2;
        ev.v[5] += 0.5*X1*fc2;
      }
    }, ev);

  if (evflag && vflag_global)
    for (int k = 0; k < 6; k++) virial[k] += ev.v[k];

  // update the position of the geometric center of each body

  if (setx) {
    const int nbody = nlocal_body + nghost_body;
    auto l_body2 = d_body;
    Kokkos::parallel_for("rigid/small/kk set_xv xgc",
      Kokkos::RangePolicy<DeviceType>(0, nbody),
      KOKKOS_LAMBDA(const int ibody) {
        Body &b = l_body2(ibody);
        b.xgc[0] = b.ex_space[0]*b.xgc_body[0] + b.ey_space[0]*b.xgc_body[1] + b.ez_space[0]*b.xgc_body[2] + b.xcm[0];
        b.xgc[1] = b.ex_space[1]*b.xgc_body[0] + b.ey_space[1]*b.xgc_body[1] + b.ez_space[1]*b.xgc_body[2] + b.xcm[1];
        b.xgc[2] = b.ex_space[2]*b.xgc_body[0] + b.ey_space[2]*b.xgc_body[1] + b.ez_space[2]*b.xgc_body[2] + b.xcm[2];
      });
  }

  if (setx) atomKK->modified(execution_space, X_MASK | V_MASK);
  else atomKK->modified(execution_space, V_MASK);
}

/* ----------------------------------------------------------------------
   copy the host body[] list to the device and back
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::copy_body_to_device()
{
  const int n = nlocal_body + nghost_body;
  if ((int) d_body.extent(0) < n) Kokkos::resize(d_body, n);
  auto h_body = Kokkos::create_mirror_view(d_body);
  for (int i = 0; i < n; i++) h_body(i) = body[i];
  Kokkos::deep_copy(d_body, h_body);
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::copy_body_to_host()
{
  const int n = nlocal_body + nghost_body;
  auto h_body = Kokkos::create_mirror_view(d_body);
  Kokkos::deep_copy(h_body, d_body);
  for (int i = 0; i < n; i++) body[i] = h_body(i);
}

/* ----------------------------------------------------------------------
   device forward/reverse communication of body data.
   Packing uses a fixed stride per atom (the upstream device fix comm sends
   nsize*sendnum and receives nsize*recvnum), so non-owner atom slots are
   simply skipped on both ends (their buffer contents are never read).
------------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                          DAT::tdual_double_1d &k_buf, int /*pbc_flag*/, int * /*pbc*/)
{
  auto d_sendlist = k_sendlist.view<DeviceType>();
  auto d_buf = k_buf.view<DeviceType>();
  auto l_body = d_body;
  auto l_bodyown = d_bodyown;

  if (commflag == INITIAL) {
    Kokkos::parallel_for("rigid/small/kk pack_forward INITIAL",
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int i) {
        const int j = d_sendlist(i);
        if (l_bodyown(j) < 0) return;
        const Body &b = l_body(l_bodyown(j));
        int m = 29*i;
        d_buf(m++) = b.xcm[0]; d_buf(m++) = b.xcm[1]; d_buf(m++) = b.xcm[2];
        d_buf(m++) = b.xgc[0]; d_buf(m++) = b.xgc[1]; d_buf(m++) = b.xgc[2];
        d_buf(m++) = b.vcm[0]; d_buf(m++) = b.vcm[1]; d_buf(m++) = b.vcm[2];
        d_buf(m++) = b.quat[0]; d_buf(m++) = b.quat[1]; d_buf(m++) = b.quat[2]; d_buf(m++) = b.quat[3];
        d_buf(m++) = b.omega[0]; d_buf(m++) = b.omega[1]; d_buf(m++) = b.omega[2];
        d_buf(m++) = b.ex_space[0]; d_buf(m++) = b.ex_space[1]; d_buf(m++) = b.ex_space[2];
        d_buf(m++) = b.ey_space[0]; d_buf(m++) = b.ey_space[1]; d_buf(m++) = b.ey_space[2];
        d_buf(m++) = b.ez_space[0]; d_buf(m++) = b.ez_space[1]; d_buf(m++) = b.ez_space[2];
        d_buf(m++) = b.conjqm[0]; d_buf(m++) = b.conjqm[1]; d_buf(m++) = b.conjqm[2]; d_buf(m++) = b.conjqm[3];
      });
    return 29*n;

  } else if (commflag == FINAL) {
    Kokkos::parallel_for("rigid/small/kk pack_forward FINAL",
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int i) {
        const int j = d_sendlist(i);
        if (l_bodyown(j) < 0) return;
        const Body &b = l_body(l_bodyown(j));
        int m = 10*i;
        d_buf(m++) = b.vcm[0]; d_buf(m++) = b.vcm[1]; d_buf(m++) = b.vcm[2];
        d_buf(m++) = b.omega[0]; d_buf(m++) = b.omega[1]; d_buf(m++) = b.omega[2];
        d_buf(m++) = b.conjqm[0]; d_buf(m++) = b.conjqm[1]; d_buf(m++) = b.conjqm[2]; d_buf(m++) = b.conjqm[3];
      });
    return 10*n;
  }
  return 0;
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first,
                                                                DAT::tdual_double_1d &k_buf)
{
  auto d_buf = k_buf.view<DeviceType>();
  auto l_body = d_body;
  auto l_bodyown = d_bodyown;
  const int l_first = first;

  if (commflag == INITIAL) {
    Kokkos::parallel_for("rigid/small/kk unpack_forward INITIAL",
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int i) {
        const int ii = l_first + i;
        if (l_bodyown(ii) < 0) return;
        Body &b = l_body(l_bodyown(ii));
        int m = 29*i;
        b.xcm[0] = d_buf(m++); b.xcm[1] = d_buf(m++); b.xcm[2] = d_buf(m++);
        b.xgc[0] = d_buf(m++); b.xgc[1] = d_buf(m++); b.xgc[2] = d_buf(m++);
        b.vcm[0] = d_buf(m++); b.vcm[1] = d_buf(m++); b.vcm[2] = d_buf(m++);
        b.quat[0] = d_buf(m++); b.quat[1] = d_buf(m++); b.quat[2] = d_buf(m++); b.quat[3] = d_buf(m++);
        b.omega[0] = d_buf(m++); b.omega[1] = d_buf(m++); b.omega[2] = d_buf(m++);
        b.ex_space[0] = d_buf(m++); b.ex_space[1] = d_buf(m++); b.ex_space[2] = d_buf(m++);
        b.ey_space[0] = d_buf(m++); b.ey_space[1] = d_buf(m++); b.ey_space[2] = d_buf(m++);
        b.ez_space[0] = d_buf(m++); b.ez_space[1] = d_buf(m++); b.ez_space[2] = d_buf(m++);
        b.conjqm[0] = d_buf(m++); b.conjqm[1] = d_buf(m++); b.conjqm[2] = d_buf(m++); b.conjqm[3] = d_buf(m++);
      });

  } else if (commflag == FINAL) {
    Kokkos::parallel_for("rigid/small/kk unpack_forward FINAL",
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int i) {
        const int ii = l_first + i;
        if (l_bodyown(ii) < 0) return;
        Body &b = l_body(l_bodyown(ii));
        int m = 10*i;
        b.vcm[0] = d_buf(m++); b.vcm[1] = d_buf(m++); b.vcm[2] = d_buf(m++);
        b.omega[0] = d_buf(m++); b.omega[1] = d_buf(m++); b.omega[2] = d_buf(m++);
        b.conjqm[0] = d_buf(m++); b.conjqm[1] = d_buf(m++); b.conjqm[2] = d_buf(m++); b.conjqm[3] = d_buf(m++);
      });
  }
}

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_reverse_comm_kokkos(int n, int first,
                                                             DAT::tdual_double_1d &k_buf)
{
  auto d_buf = k_buf.view<DeviceType>();
  auto l_body = d_body;
  auto l_bodyown = d_bodyown;
  const int l_first = first;

  // FORCE_TORQUE only
  Kokkos::parallel_for("rigid/small/kk pack_reverse",
    Kokkos::RangePolicy<DeviceType>(0, n),
    KOKKOS_LAMBDA(const int i) {
      const int ii = l_first + i;
      if (l_bodyown(ii) < 0) return;
      const Body &b = l_body(l_bodyown(ii));
      int m = 6*i;
      d_buf(m++) = b.fcm[0]; d_buf(m++) = b.fcm[1]; d_buf(m++) = b.fcm[2];
      d_buf(m++) = b.torque[0]; d_buf(m++) = b.torque[1]; d_buf(m++) = b.torque[2];
    });
  return 6*n;
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::unpack_reverse_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                                                DAT::tdual_double_1d &k_buf)
{
  auto d_sendlist = k_sendlist.view<DeviceType>();
  auto d_buf = k_buf.view<DeviceType>();
  auto l_body = d_body;
  auto l_bodyown = d_bodyown;

  // FORCE_TORQUE only: add ghost contributions into the owning body
  Kokkos::parallel_for("rigid/small/kk unpack_reverse",
    Kokkos::RangePolicy<DeviceType>(0, n),
    KOKKOS_LAMBDA(const int i) {
      const int j = d_sendlist(i);
      if (l_bodyown(j) < 0) return;
      Body &b = l_body(l_bodyown(j));
      int m = 6*i;
      b.fcm[0] += d_buf(m++); b.fcm[1] += d_buf(m++); b.fcm[2] += d_buf(m++);
      b.torque[0] += d_buf(m++); b.torque[1] += d_buf(m++); b.torque[2] += d_buf(m++);
    });
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   host-fallback diagnostics/operations: sync the device-resident body data
   to the host (and back, for the ones that modify it) around the base call
------------------------------------------------------------------------- */

template<class DeviceType>
double FixRigidSmallKokkos<DeviceType>::compute_scalar()
{
  if (body_resident_device) copy_body_to_host();
  return FixRigidSmall::compute_scalar();
}

template<class DeviceType>
void *FixRigidSmallKokkos<DeviceType>::extract(const char *str, int &dim)
{
  if (body_resident_device) copy_body_to_host();
  return FixRigidSmall::extract(str, dim);
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::deform(int flag)
{
  if (body_resident_device) copy_body_to_host();
  FixRigidSmall::deform(flag);
  if (body_resident_device) copy_body_to_device();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::zero_momentum()
{
  if (body_resident_device) copy_body_to_host();
  atomKK->sync(Host, datamask_read);
  FixRigidSmall::zero_momentum();
  atomKK->modified(Host, datamask_modify);
  atomKK->sync(execution_space, datamask_modify);
  if (body_resident_device) copy_body_to_device();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::zero_rotation()
{
  if (body_resident_device) copy_body_to_host();
  atomKK->sync(Host, datamask_read);
  FixRigidSmall::zero_rotation();
  atomKK->modified(Host, datamask_modify);
  atomKK->sync(execution_space, datamask_modify);
  if (body_resident_device) copy_body_to_device();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::write_restart_file(const char *file)
{
  if (body_resident_device) copy_body_to_host();
  atomKK->sync(Host, datamask_read);
  FixRigidSmall::write_restart_file(file);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixRigidSmallKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixRigidSmallKokkos<LMPHostType>;
#endif
}
