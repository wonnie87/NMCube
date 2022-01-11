EXE = main_rk.exe
#OBJ = phdf5_helpers.o sub_rk.o main_rk.o
OBJ = hdf5_helpers.o sub_rk.o main_rk.o

$(EXE): $(OBJ)
	$(MPIF) -o $(EXE) $(OBJ) $(I_FLAG) $(L_FLAG)

hdf5_helpers.mod: hdf5_helpers.o

hdf5_helpers.o: hdf5_helpers.f90
	$(MPIF) $(MPIF_FLAG) -c hdf5_helpers.f90 $(I_FLAG)

sub_rk.mod: sub_rk.o

sub_rk.o: sub_rk.f90
	$(MPIF) $(MPIF_FLAG) -c sub_rk.f90 $(I_FLAG)

#main_rk.o: phdf5_helpers.mod sub_rk.mod main_rk.f90
main_rk.o: hdf5_helpers.mod sub_rk.mod main_rk.f90
	$(MPIF) $(MPIF_FLAG) -c main_rk.f90 $(I_FLAG)

.PHONY: run clean reset

run:
#	$(py) $(INP)
	$(MPIRUN) $(MPIRUN_FLAG) $(NP) ./$(EXE)

clean:
	rm -f *.o

reset:
	rm -f *.o *.mod $(EXE)
