EXE = main_nb.exe
#OBJ = phdf5_helpers.o sub_nb.o main_nb.o
OBJ = hdf5_helpers.o sub_nb.o main_nb.o

$(EXE): $(OBJ)
	$(MPIF) -o $(EXE) $(OBJ) $(I_FLAG) $(L_FLAG)

#phdf5_helpers.mod: phdf5_helpers.o
#
#phdf5_helpers.o: hdf5_helpers.mod phdf5_helpers.f90 
#	$(MPIF) $(MPIF_FLAG) -c phdf5_helpers.f90 $(I_FLAG)

hdf5_helpers.mod: hdf5_helpers.o

hdf5_helpers.o: hdf5_helpers.f90
	$(MPIF) $(MPIF_FLAG) -c hdf5_helpers.f90 $(I_FLAG)

sub_nb.mod: sub_nb.o

sub_nb.o: sub_nb.f90
	$(MPIF) $(MPIF_FLAG) -c sub_nb.f90 $(I_FLAG)

#main_nb.o: phdf5_helpers.mod sub_nb.mod main_nb.f90
main_nb.o: hdf5_helpers.mod sub_nb.mod main_nb.f90
	$(MPIF) $(MPIF_FLAG) -c main_nb.f90 $(I_FLAG)

.PHONY: run clean reset

run:
#	$(py) $(INP)
	$(MPIRUN) $(MPIRUN_FLAG) $(NP) ./$(EXE)

clean:
	rm -f *.o *.stat

reset:
	rm -f *.o *.mod *.stat $(EXE)
