EXE = main_NB.exe
OBJ = phdf5_helpers.o sub_NB.o main_NB.o
MPIF = mpif90
MPIF_FLAG = -fdefault-real-8
#MPIF_FLAG = -r8 ## uncomment if using Intel complier
#I_FLAG = -I${HDF5_DIR}/include/
L_FLAG = -lhdf5_fortran -lhdf5 -lz
#L_FLAG = -L${HDF5_DIR}/lib -lhdf5_fortran -lhdf5 -lz # specify path to HDF5 library if not in the default path.
MPIRUN = mpirun
MPIRUN_FLAG = -np
NP = 4
py = python3
INP = InputFileGeneration_NB.py

$(EXE): $(OBJ)
	$(MPIF) -o $(EXE) $(OBJ) $(I_FLAG) $(L_FLAG)

phdf5_helpers.mod: phdf5_helpers.o

phdf5_helpers.o: hdf5_helpers.mod phdf5_helpers.f90 
	$(MPIF) $(MPIF_FLAG) -c phdf5_helpers.f90 $(I_FLAG)

hdf5_helpers.mod: hdf5_helpers.o

hdf5_helpers.o: hdf5_helpers.f90
	$(MPIF) $(MPIF_FLAG) -c hdf5_helpers.f90 $(I_FLAG)

sub_NB.mod: sub_NB.o

sub_NB.o: sub_NB.f90
	$(MPIF) $(MPIF_FLAG) -c sub_NB.f90 $(I_FLAG)

main_NB.o: phdf5_helpers.mod sub_NB.mod main_NB.f90
	$(MPIF) $(MPIF_FLAG) -c main_NB.f90 $(I_FLAG)

.PHONY: run clean reset

run:
	$(py) $(INP)
	$(MPIRUN) $(MPIRUN_FLAG) $(NP) ./$(EXE)

clean:
	rm -f *.o

reset:
	rm -f *.o *.mod $(EXE)
