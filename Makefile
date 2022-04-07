PY = python3
MPIF = mpif90
#MPIF = /usr/lib64/openmpi/bin/mpif90 # Default binary path on Fedora
#MPIF = mpiifort # uncomment if using Intel complier
ifeq ($(strip $(MPIF)), mpiifort)
	MPIF_FLAG = -r8
else
	MPIF_FLAG = -fdefault-real-8
endif
#I_FLAG = -I${DIR_TO_HDF5_HEADER} # specify directory to HDF5 library if not in default PATH.
#I_FLAG = -I/usr/include/hdf5/openmpi # Example path for Deb-based Linux
#I_FLAG = -I/usr/lib64/gfortran/modules/openmpi # Example path for RHEL-based Linux
#I_FLAG = -I/usr/local/Cellar/hdf5-mpi/1.12.1/include # Example path for MacOS 
L_FLAG = -lhdf5_fortran -lhdf5
#L_FLAG = -L${DIR_TO_HDF5_LIB} -lhdf5_fortran -lhdf5 # specify directory to HDF5 library if not in default PATH.
#L_FLAG = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5_fortran -lhdf5 # Example path for Deb-based Linux
#L_FLAG = -L/usr/lib64/openmpi/lib -lhdf5_fortran -lhdf5 # Example path for RHEL-based Linux
#L_FLAG = -L/usr/local/Cellar/hdf5-mpi/1.12.1/lib -lhdf5_fortran -lhdf5 # Example path for MacOS
MPIRUN = mpirun
#MPIRUN = /usr/lib64/openmpi/bin/mpirun # Default binary path on Fedora
MPIRUN_FLAG = -np
NP = 4
INP = ./inputs/InputFile.py
OUT = ""
FSKIP = 100

export MPIF MPIF_FLAG I_FLAG L_FLAG MPIRUN MPIRUN_FLAG NP

.PHONY: all rk nb run clean reset

all: rk nb

rk:
	$(MAKE) --directory=./$@

nb:
	$(MAKE) --directory=./$@

run:
	@case $(INP) in \
		*/*) echo " >> Generating Fortran inputs from $(INP)"; $(PY) $(INP) ;;\
		*) echo " >> Generating Fortran inputs from ./inputs/$(INP)"; $(PY) ./inputs/$(INP);;\
	esac; 
	@echo -e " >> Fortran inputs are generated successfully.\n"
	@$(MAKE) skipIG

skipIG:
	@if [ `sed -n '2p' numMethod.inp` -lt 50 ]; \
		then $(MAKE) run --directory=./rk;\
	else $(MAKE) run --directory=./nb;\
	fi;

check:
	@$(PY) -c "from modules.postprocessing import NumStabilityChk; NumStabilityChk($(NP))";
#	$(PY) ./modules/NumStabilityCheck.py $(sed -n '2p' design.py)

movie:
	@case $(OUT) in \
		"" ) echo " >> ERROR: Specify the output file!!!" ;\
			echo " >> For example, make movie OUT=./outputs/P4_d04_F0.1_f8.0_RK4.h5" ;;\
		*/*) echo " >> Generating a video from $(OUT)";\
			$(PY) -c 'from modules.postprocessing import MovieGen; MovieGen("$(OUT)",$(FSKIP))' ;;\
		*) echo " >> Generating a video from ./outputs/$(OUT)";\
			$(PY) -c 'from modules.postprocessing import MovieGen; MovieGen("./outputs/$(OUT)",$(FSKIP))' ;;\
	esac;

clean:
	$(MAKE) clean --directory=./rk
	$(MAKE) clean --directory=./nb
	rm -f *.inp

reset:
	$(MAKE) reset --directory=./rk
	$(MAKE) reset --directory=./nb
	rm -f *.inp

