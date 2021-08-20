# **NM^3**: Nonlinear MetaMaterial MPI Solver
This is the repository for the source codes and example files used in paper:

**High-Performance Large-Scale Simulation of Multi-stable Metastructures**

M. Hwang, C. Scalo, and A. F. Arrieta, *To be published*.  

## Program overview
**NM^3** is massively parallel code for solving general multi-stable mechanical metamaterials and metastructures.
The code is written in **Fortran** and utilizes **MPI** (Message-passing interface) system calls to execute parallel computation.
The code supports both explicit and implicit numerical methods: up to the 4<sup>th</sup>-order Runge-Kutta (RK) methods are implemented for the explicit solver; and the Newmark-&beta; (NB) method with constant average acceleration is implemented for the implicit solver.
Currently, three types of multi-stable systems can be solved with **NM^3**: (1) 1D multi-stable lattice with coupled pendula (discrete sine-Gordon model), (2) 1D lattice with quartic on-site potentials (discerete &phi;-4 model), and (3) metabeam with a bi-stable microsctructure. 

## Before running the program
* Install **MPI**, **Fortran** complier, **GNU make**, **Python**, **git**, **hdf5**.
* Clone **NM^3** repository.
  > $ git clone https://github.com/wonnie87/NMCube  
* Compile and link the programs
  > $ cd NMCube/NB  
  > $ make  
  > $ cd ../RK  
  > $ make

## Running the program
For both NB and RK solvers, run the following command in their respective folders.
> $ make run NP=(number of processes) INP=(input file name)  

For example, to run the input file **InpGen_prob4.py** with RK solver:
> $ cd NMCube/RK
> $ make run NP=8 INP=InpGen_prob4.py
 
If parameters **NP** and **INP** are not specified, **make run** uses **NP**=4 and **INP**=InputFileGeneration.py by default. 

## Data storage
The simulated outputs are packaged in a **HDF5** file format (.h5).
The time displacements, velocities, and accelerations are stored in the datasets **/x**, **/xdot**, and **/xddot**, respectively.
* (Note:) For NB solver, the outputs at each time step are stored in separate datasets. The data structure will be modified to follow that of RK solver in the next revision.

Use **h5dump** or **h5ls** (default command line programs supplied with **hdf5** installation) for quick data retrieval.

For example, to view the entire time displacements of a file **d04_F0.1_8.0Hz.h5**:
> $ h5dump -d "/x" d04_F0.1_8.0Hz.h5

For more involved postprocessing, use more advanced tools, such as **Python**'s **h5py** library.
Example **Python** scripts (**postprocess.py**) can be found in each folder.
