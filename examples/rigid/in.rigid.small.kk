# Rigid-body melt demonstrating fix rigid/small and its KOKKOS variant.
#
# Each molecule (a small rigid cluster built on the lattice) is integrated as
# an independent rigid body.  Run on the CPU with
#
#   lmp -in in.rigid.small.kk
#
# or on the device (GPU/OpenMP) with the KOKKOS accelerated rigid/small/kk:
#
#   lmp -k on g 1 -sf kk -in in.rigid.small.kk
#
# fix rigid/small integrates the rigid-body equations of motion (NVE for the
# bodies). The thermodynamic output is identical (to round-off) between the
# CPU and the KOKKOS rigid/small/kk runs and across MPI ranks.

units           lj
atom_style      molecular

lattice         sc 0.8
region          box block 0 12 0 12 0 12
create_box      1 box
create_atoms    1 box
mass            1 1.0

# group atoms into rigid clusters of 4 (two lattice cells along x) by molecule ID

variable        molid atom floor((id-1)/4)+1
set             atom * mol v_molid

velocity        all create 1.0 4928459 loop geom

pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0

neighbor        0.3 bin
neigh_modify    every 1 delay 0 check yes

timestep        0.001
fix             1 all rigid/small molecule

thermo          50
thermo_style    custom step temp pe etotal press
run             500
