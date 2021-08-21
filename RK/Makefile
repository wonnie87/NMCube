EXE = main_RK.exe
#OBJ = phdf5_helpers.o sub_RK.o main_RK.o
OBJ = hdf5_helpers.o sub_RK.o main_RK.o
MPIF = mpif90
MPIF_FLAG = -fdefault-real-8
#MPIF_FLAG = -r8 # uncomment if using Intel complier
#I_FLAG = -I${PATH_TO_HDF5_HEADER} # specify path to HDF5 library if not in the default path.
L_FLAG = -lhdf5_fortran -lhdf5 -lz
#L_FLAG = -L${PATH_TO_HDF5_LIB} -lhdf5_fortran -lhdf5 -lz # specify path to HDF5 library if not in the default path.
MPIRUN = mpirun
MPIRUN_FLAG = -np
NP = 4
py = python3
INP = InputFileGeneration.py

$(EXE): $(OBJ)
	$(MPIF) -o $(EXE) $(OBJ) $(I_FLAG) $(L_FLAG)

#phdf5_helpers.mod: phdf5_helpers.o
#
#phdf5_helpers.o: hdf5_helpers.mod phdf5_helpers.f90 
#	$(MPIF) $(MPIF_FLAG) -c phdf5_helpers.f90 $(I_FLAG)

hdf5_helpers.mod: hdf5_helpers.o

hdf5_helpers.o: hdf5_helpers.f90
	$(MPIF) $(MPIF_FLAG) -c hdf5_helpers.f90 $(I_FLAG)

sub_RK.mod: sub_RK.o

sub_RK.o: sub_RK.f90
	$(MPIF) $(MPIF_FLAG) -c sub_RK.f90 $(I_FLAG)

#main_RK.o: phdf5_helpers.mod sub_RK.mod main_RK.f90
main_RK.o: hdf5_helpers.mod sub_RK.mod main_RK.f90
	$(MPIF) $(MPIF_FLAG) -c main_RK.f90 $(I_FLAG)

.PHONY: run clean reset

run:
	$(py) $(INP)
	$(MPIRUN) $(MPIRUN_FLAG) $(NP) ./$(EXE)

clean:
	rm -f *.o

reset:
	rm -f *.o *.mod $(EXE)
