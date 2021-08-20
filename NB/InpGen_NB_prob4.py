import numpy as np

##### Input parameters for numerical method #####
gamma_NB = 0.5
beta_NB = 0.25
tol_NR = 1.0e-16 # convergent criteria for NR iteration
tol_inv = 1.0e-16 # convergent criteria for Conjugate Gradient method
N_NR = 200
N_inv = 200
dt = 1.0e-3 # time step
tStart = 0.0 # start of simulation time
tEnd = 6.0 # end of simulation time
dtWrite = 1.0e-3 # file wirte frequency

##### Input parameters for system design #####
des = 'd04' # design name (will be the prefix of the output file)
N_global = 30 # Number of unit cells
prob_flag = 4 # 1: 1D pendula chain, 2: 1D phi4 chain, 3: N/A, 4: bistable metabeam

## BC - Boundary condition 
    # [Inner unit cells, Left end, Right end]
    # First digit (unit cell type): 1-normal 
    # Second digit (unit cell): 0-internal units; 1-leftmost unit; 2-rightmost unit
    # Third digit (condition): 0-free, 1-fixed
BC = np.array([100, 111, 120])

#G = 1. # [for prob 1] gravitational constant
#m = 1. # [for prob 1] pendulum mass

## k (stiffness) array
    # prob 1: [k_th (torsional spring), G, m] - [ToDo] Separate G and m
    # prob 2: [k, C1, C2, C3, C4] - refer to the manuscript
    # prob 4: [k1, k2, k3, k4, k5, k6, k7, k8] - refer to the manuscript
#k = np.array([1., G, m]) # for prob 1
#k = np.array([1.0, 0.0, -0.030, 0.0, 0.015]) # for prob 2
k = np.array([1.241, 0.6*1.793, 0.6, 100., 100., 100., 100., 100.]) # for prob 4

## L (length) array
    # prob 1: [pendulum lengh]
    # prob 2: [] - empty
    # prob 4: [L1, L2, L3, R] - Refer to the manuscript
#L = np.array([1.]) # for prob 1
L = np.array([20., 40., 20., 8.]) # for prob 4

## LC (Load condition) array: [[UC, local DoF, Loadcase], ... ]
    # Loadcases: 1-sine force, 2-modulated sine force,
    #   11-sine disp input, 12-modulated sine disp input, 
    #   13-long-stroke harmonic disp input
    # Accompanied by the corresponding LC_val array
    #   e.g., LC_val = [tForceStart, tForceEnd, amp, freqIn, phi0, c1, c2] for Loadcase 1
LC = np.array([[1,1,1]]) # Load cases: [[UC, local DoF, Loadcase], ... ]
amp = 0.1
freqIn = 8.0
phi0 = 0.0
t1 = 0.0
t2 = 0.0
LC_val = np.array([[tStart, tEnd, amp, freqIn, phi0, t1, t2]])

## DOF - number of DoFs for each unit cell
    # 1 for prob1; 1 for prob2; 1 for prob3; 6 for prob4;
if prob_flag == 1:
    DoF = 1
elif prob_flag == 2:
    DoF = 1
elif prob_flag == 4:
    DoF = 6

## mass array
    # prob 1: [I, m]
    # prob 2: [m]
    # prob 4: [m1, m1, m2, m2, m3, m3]
#m = 1. * np.ones((N_global,DoF)) # for prob 1
#m = 2. * np.ones((N_global,DoF)) # for prob 2
m = np.array([2.e-6, 2.e-6, 1.e-6, 1.e-6, 1.e-6, 1.e-6]) * np.ones((N_global,DoF)) # for prob 4

## on-site damping array
#b = 0.0 * np.ones((N_global,DoF))
zeta = 0.03
freq0 = 26.281
b = 2*zeta*2*np.pi*freq0*m

## Initial conditions
u = np.zeros((N_global,DoF))
udot = np.zeros((N_global,DoF))
#udot[0] = 2.0 # imposed initial velocity at the first unit

##### END of user inputs #####


########################################################
outF1 = open("numMethod_NB.inp", "w")
outF2 = open("design.inp", "w")

outF1.write("**gamma_NB\n")
outF1.write("{}\n".format(gamma_NB))
outF1.write("**beta_NB\n")
outF1.write("{}\n".format(beta_NB))
outF1.write("**tol_NR\n")
outF1.write("{}\n".format(tol_NR))
outF1.write("**tol_inv\n")
outF1.write("{}\n".format(tol_inv))
outF1.write("**N_NR\n")
outF1.write("{}\n".format(N_NR))
outF1.write("**N_inv\n")
outF1.write("{}\n".format(N_inv))

outF1.write("**dt\n")
outF1.write("{}\n".format(dt))
outF1.write("**tStart\n")
outF1.write("{}\n".format(tStart))
outF1.write("**tEnd\n")
outF1.write("{}\n".format(tEnd))
outF1.write("**dtWrite\n")
outF1.write("{}\n".format(dtWrite))

outF2.write("**filename\n")
if LC[0,2] == 1:
    outF2.write("{0:s}_F{1:}_{2:}Hz.h5\n".format(des,amp,freqIn))
elif LC[0,2] == 2:
    outF2.write("{0:s}_P{1:}_{2:}s.h5\n".format(des,amp,t2-t1))
elif LC[0,2] == 11 or LC[0,2] == 13:
    outF2.write("{0:s}_A{1:}_{2:}Hz.h5\n".format(des,amp,freqIn))

outF2.write("**N_global\n")
outF2.write("{}\n".format(N_global))
outF2.write("**prob_flag\n")
outF2.write("{}\n".format(prob_flag))
outF2.write("**BC(0), BC(1), BC(2)\n")
for n in range(len(BC)):
    if n == len(BC)-1:
        outF2.write("{}\n".format(BC[n]))
    else:
        outF2.write("{}, ".format(BC[n]))

if prob_flag == 1:
    outF2.write("**k(1), G, m\n")
elif prob_flag == 2:
    outF2.write("**k(1), k(2), k(3), k(4), k(5)\n")
elif prob_flag == 4:
    outF2.write("**k(1), k(2), k(3), k(4), k(5), k(6), k(7), k(8)\n")

for n in range(len(k)):
    if n == len(k)-1:
        outF2.write("{}\n".format(k[n]))
    else:
        outF2.write("{}, ".format(k[n]))

if prob_flag != 2:
    if prob_flag == 1:
        outF2.write("**L(1)\n")
    elif prob_flag ==4:
        outF2.write("**L(1), L(2), L(3), L(4)\n")

    for n in range(len(L)):
        if n == len(L)-1:
            outF2.write("{}\n".format(L[n]))
        else:
            outF2.write("{}, ".format(L[n]))

LC_dim2, LC_dim1 = LC.shape
outF2.write("**LC_dim2\n")
outF2.write("{}\n".format(LC_dim2))
outF2.write("**Load cases\n")
for i in range(LC_dim2):
    for j in range(LC_dim1):
        if j == LC_dim1-1:
            outF2.write("{}\n".format(LC[i,j]))
        else:
            outF2.write("{}, ".format(LC[i,j]))

outF2.write("**Corresponding LC parameters\n")
for i in range(LC_dim2):
    for j in range(LC_val.shape[1]):
        if j == LC_val.shape[1]-1:
            outF2.write("{}\n".format(LC_val[i,j]))
        else:
            outF2.write("{}, ".format(LC_val[i,j]))

outF2.write("**DoF\n")
outF2.write("{}\n".format(DoF))
if prob_flag == 1:
    outF2.write("**m(1), m(2)\n")
elif prob_flag == 2:
    outF2.write("**m(1)\n")
elif prob_flag == 4:
    outF2.write("**m(1), m(2), m(3), m(4), m(5), m(6)\n")

for i in range(N_global):
    for j in range(m.shape[1]):
        if j == m.shape[1]-1:
            outF2.write("{}\n".format(m[i,j]))
        else:
            outF2.write("{}, ".format(m[i,j]))

if prob_flag == 1:
    outF2.write("**b(1)\n")
elif prob_flag == 2:
    outF2.write("**b(1)\n")
elif prob_flag == 4:
    outF2.write("**b(1), b(2), b(3), b(4), b(5), b(6)\n")

for i in range(N_global):
    for j in range(b.shape[1]):
        if j == b.shape[1]-1:
            outF2.write("{}\n".format(b[i,j]))
        else:
            outF2.write("{}, ".format(b[i,j]))

outF2.write("**u0\n")
for i in range(N_global):
    for j in range(u.shape[1]):
        if j == u.shape[1]-1:
            outF2.write("{}\n".format(u[i,j]))
        else:
            outF2.write("{}, ".format(u[i,j]))

outF2.write("**udot0\n")
for i in range(N_global):
    for j in range(udot.shape[1]):
        if j == udot.shape[1]-1:
            outF2.write("{}\n".format(udot[i,j]))
        else:
            outF2.write("{}, ".format(udot[i,j]))

outF1.close()
outF2.close()
