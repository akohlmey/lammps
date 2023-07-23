Static LAMMPS universal binaries for MacOS X (arm64/x86_64)
===========================================================

This package provides static binaries of LAMMPS that should run on most
current MacOS X systems.  Note the binaries are serial and only enable a subset
of the available packages.

After copying the LAMMPS folder into your Applications folder, please follow
these steps:

1. Open the Terminal app

2. Type the following command and press ENTER:

   open ~/.zprofile

   This will open a text editor for modifying the .zprofile file in your home
   directory. 

3. Add the following lines to the end of the file, save it, and close the editor

   export LAMMPS_INSTALL_DIR=/Applications/LAMMPS
   export LAMMPS_POTENTIALS=/Applications/LAMMPS/share/lammps/potentials
   export LAMMPS_BENCH_DIR=$LAMMPS_INSTALL_DIR/bench
   export PATH=${LAMMPS_INSTALL_DIR}/bin:$PATH

4. In your existing terminal, type the following command make the settings active

   source ~/.zprofile

   Note, you don't have to type this in new terminals, since they will apply
   the changes from .zprofile automatically.

5. Try running LAMMPS

   lmp -i $LAMMPS_BENCH_DIR/in.lj
