# LAMMPS source -*- Makefile -*-

SHELL = /bin/bash
PYTHON = $(word 3,$(shell type python3))

help:
	@echo ' '
	@echo 'make check-whitespace    check source files for tabs and trailing whitespace'
	@echo 'make check-permissions   check source files for correct permissions'
	@echo 'make check-homepage      check that source files correctly refer to the LAMMPS homepage'
	@echo 'make check-errordocs     check whether source files contain legacy error docs'
	@echo 'make check-fmtlib        check source files for fmt::format features incompatible with std::format'
	@echo 'make check-docs          displace doc files status and check for common issues'
	@echo ' '
	@echo 'make check               run all check-<something> targets from above'
	@echo ' '
	@echo 'make fix-whitespace      remove tabs and trailing whitespace from source files'
	@echo 'make fix-permissions     adjust permissions of source files'
	@echo 'make fix-homepage        adjust LAMMPS homepage URL in source files'
	@echo 'make fix-errordocs       remove legacy error docs from source files'
	@echo 'make fix-fmtlib          replace fmt::format features with known alternatives for std::format compatibility'
	@echo ' '
	@echo 'make check-version       list source and doc files with pending version tags

# deprecation warning

warning clean-all clean-serial clean-mpi mpi-stubs install-python tar package package-status ps package-installed pi yes-all no-all yes-basic no-basic yes-most no-most package-update pu package-overwrite package-diff pd purge serial mpi:
	@echo ' '
	@echo ' Support for the legacy build system has been removed from LAMMPS'
	@echo ' Please use CMake instead. See https://docs.lammps.org/Build_cmake.html'
	@echo ''
	@echo ' Contact us at developers@lammps.org if you have concerns or questions'
	@echo ' '
	@exit 1

check: check-whitespace check-permissions check-homepage check-errordocs check-fmtlib check-docs

check-whitespace:
	@$(PYTHON) ../tools/coding_standard/whitespace.py ..

fix-whitespace:
	@$(PYTHON) ../tools/coding_standard/whitespace.py .. -f

check-permissions:
	@$(PYTHON) ../tools/coding_standard/permissions.py ..

fix-permissions:
	@$(PYTHON) ../tools/coding_standard/permissions.py .. -f

check-homepage:
	@$(PYTHON) ../tools/coding_standard/homepage.py ..

fix-homepage:
	@$(PYTHON) ../tools/coding_standard/homepage.py .. -f

check-errordocs:
	@$(PYTHON) ../tools/coding_standard/errordocs.py ..

fix-errordocs:
	@$(PYTHON) ../tools/coding_standard/errordocs.py .. -f

check-fmtlib:
	@$(PYTHON) ../tools/coding_standard/fmtlib.py ..

fix-fmtlib:
	@$(PYTHON) ../tools/coding_standard/fmtlib.py .. -f

check-docs:
	@$(MAKE) $(MFLAGS) -C ../doc anchor_check style_check package_check role_check

check-version:
	@$(PYTHON) ../tools/coding_standard/versiontags.py .. || echo

clean:
	@rm -rf Obj_* lmp_* liblammps*
	@rm -fv style_*.h style_*.cpp packages_*.h package_registry.cpp lmpgitversion.h lmpinstalledpkgs.h
	@rm -fv *.*~ *.*.rej *.*.orig */*.*~ */*.*.rej */*.*.orig

.PHONY: warning help check check-whitespace check-permissions check-homepage check-errordocs check-fmtlib check-docs fix-whitespace fix-permissions fix-homepage fix-errordocs fix-fmtlib serial mpi clean


