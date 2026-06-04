/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// fix gcmc 'fast yes' evaluates each molecule insertion/deletion trial's
// energy change locally instead of recomputing the whole-system energy.
// Because it produces the same per-trial energy (to roundoff), the accept /
// reject decisions -- and therefore the entire trajectory -- must be identical
// to the default 'fast off' (full energy) path.  This test runs the same short
// GCMC simulation both ways and asserts the final loading and potential energy
// are bit-for-bit identical.

#include "../testing/core.h"

#include "atom.h"
#include "gtest/gtest.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "variable.h"

#include <cmath>
#include <cstdlib>
#include <mpi.h>
#include <string>

#ifndef TEST_INPUT_FOLDER
#define TEST_INPUT_FOLDER .
#endif
#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

bool verbose = false;

namespace LAMMPS_NS {

class GCMCFastTest : public LAMMPSTest {
protected:
    // run the short GCMC simulation with the given 'fast' setting and return
    // the final number of atoms and the final potential energy
    void run_gcmc(const std::string &fast, bigint &natoms, double &pe)
    {
        HIDE_OUTPUT([&] {
            command("clear");
            command("units real");
            command("atom_style full");
            command("boundary p p p");
            command("pair_style lj/cut/coul/long 8.0");
            command("pair_modify mix arithmetic tail yes");
            command("kspace_style pppm 1.0e-4");
            command("bond_style harmonic");
            command("special_bonds lj/coul 0.0 0.0 0.0");

            command("region box block 0 24 0 24 0 24 units box");
            command("create_box 3 box bond/types 1 extra/bond/per/atom 2 "
                    "extra/special/per/atom 4");

            command("lattice sc 6.0");
            command("create_atoms 1 region box");
            command("mass 1 28.0");
            command("mass 2 14.0");
            command("mass 3 14.0");
            command("set atom 1*32 charge 0.2");
            command("set atom 33*64 charge -0.2");

            command("pair_coeff 1 1 0.10 3.40");
            command("pair_coeff 2 2 0.07 3.30");
            command("pair_coeff 3 3 0.07 3.30");
            command("bond_coeff 1 100.0 1.1");

            command("molecule dimer " STRINGIFY(TEST_INPUT_FOLDER) "/mc_fast_dimer.mol");
            command("group gas type 2 3");

            command("fix mygcmc gas gcmc 25 25 0 0 12345 300.0 -0.5 0.5 mol dimer fast " +
                    fast);
            command("thermo 25");
            command("run 200 post no");
            command("variable mype equal pe");
        });

        natoms = lmp->atom->natoms;
        pe = lmp->input->variable->compute_equal(lmp->input->variable->find("mype"));
    }
};

TEST_F(GCMCFastTest, FastMatchesFullLoading)
{
    if (!info->has_style("kspace", "pppm")) GTEST_SKIP();
    if (!info->has_style("pair", "lj/cut/coul/long")) GTEST_SKIP();
    if (!info->has_package("MC")) GTEST_SKIP();

    bigint natoms_off, natoms_yes;
    double pe_off, pe_yes;
    run_gcmc("off", natoms_off, pe_off);
    run_gcmc("yes", natoms_yes, pe_yes);

    EXPECT_EQ(natoms_off, natoms_yes);
    EXPECT_DOUBLE_EQ(pe_off, pe_yes);
}

}    // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (const char *var = getenv("TEST_ARGS")) {
        if (std::string(var).find("-v") != std::string::npos) verbose = true;
    }
    if ((argc > 1) && (std::string(argv[1]) == "-v")) verbose = true;

    const int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
