import numpy as np
import os, sys
sys.path.insert(1, os.path.abspath('./'))
from modules import InpGen as ig

#################################################
##### Input parameters for numerical method #####
#################################################

#### Required inputs ####
dt = 1.0e-7 # time step
tStart = 0.0 # start of simulation time
tEnd = 0.01 # end of simulation time
dtWrite = 1.0e-7 # file wirte frequency
sol_flag = 'RK4' # choose among RK1, RK2, RK3, RK4, and NB

#### Optional inputs ####
#sol_params = {
#        'b1': 1.0
#        } # RK1
#sol_params = {
#        'c2': 1.0,
#        'a21': 1.0,
#        'b1': 1./2, 'b2': 1./2
#        } # RK2
#sol_params = {
#        'c2': 0.5, 'c3': 1.0,
#        'a21': 1.0, 'a31': -1.0, 'a32': 2.0,
#        'b1': 1./6, 'b2': 2./3, 'b3': 1./6
#        } # RK3
#sol_params = { 
#        'c2': 0.5, 'c3': 0.5, 'c4': 1.0,
#        'a21': 0.5, 'a31': 0.0, 'a32': 0.5, 'a41': 0.0, 'a42': 0.0, 'a43': 1.0,
#        'b1': 1./6, 'b2': 1./3, 'b3': 1./3, 'b4': 1./6
#        } # RK4
#sol_params = {
#        'gamma_NB': 0.5, 'beta_NB': 0.25
#        'tol_NR': 1.0e-16, 'tol_inv': 1.0e-16
#        'N_NR': 200, 'N_inv': 200
#        } # NB


##############################################
##### Input parameters for system design #####
##############################################

des = 'FSI3' # design name (will be the prefix of the output file)
N_global = 200 # Number of unit cells
prob = ig.ProbType(4) # 1: 1D pendula chain, 2: 1D phi4 chain, 3: N/A, 4: bistable metabeam

### BC - Boundary condition ###
    ## [Inner unit cells, Left end, Right end]
    ## First digit (unit cell type): 1-normal; 2-bottom-fixed;
    ## Second digit (unit cell): 0-internal units; 1-leftmost unit; 2-rightmost unit
    ## Third digit (condition): 0-free; 1-fixed; 
        #00(=0)-free-free; 01-free-fixed; 02-free-HFixed; 03-free-VFixed;
        #10-fixed-free; 11(=1)-fixed-fixed; 12-fixed-HFixed; 13-fixed-VFixed
BC = np.array([100, 2112, 2212])

#G = 1. # [for prob 1 only] gravitational constant

### k (stiffness) array ###
    ## prob 1: [k_th (stiffness of torsional spring)] 
    ## prob 2: [k, C1, C2, C3, C4] - refer to the manuscript
    ## prob 4: [k1, k2, k3, k4, k5, k6, k7, k8] - refer to the manuscript
#k = np.array([1.]) # for prob 1
#k = np.array([1.0, 0.0, -0.030, 0.0, 0.015]) # for prob 2
k = np.array([20_000, 100_000.*1.451, 100_000., 50_000., 50_000., 50_000., 50_000., 50_000.]) # for prob 4

### L (length) array ###
    # prob 1: [pendulum lengh]
    # prob 2: [] - empty
    # prob 4: [L1, L2, L3, R] - Refer to the manuscript
#L = np.array([1.]) # for prob 1
L = np.array([5., 12., 8., 2.]) # for prob 4

### LC (Load condition) array: [[UC, local DoF, Loadcase], ... ] ###
    ## Loadcases: 1-sine force, 2-modulated sine force,
    ##   11-sine disp input, 12-modulated sine disp input, 
    ##   13-long-stroke harmonic disp input
    ## Accompanied by the corresponding LC_val array
    ##   e.g., LC_val = [tForceStart, tForceEnd, amp, freqIn, phi0, c1, c2] for Loadcase 1
#LC = np.array([[1,4,1]]) # Load cases: [[UC, local DoF, Loadcase], ... ]
LC = np.empty((N_global,3),dtype=np.int8)
LC_val = np.empty((N_global,7))
amp = 1_000.
freqIn = 45_000.0
vFlow = 10_000
for n in range(0,N_global):
    LC[n,] = [1,4,1]
    phi0 = -2*np.pi*freqIn*L[0]*n/vFlow
    t1 = L[0]*n/vFlow
    t2 = tEnd
    LC_val[n,] = [t1, t2, amp, freqIn, phi0, t1, t2]

### mass array ###
    ## prob 1: [I, m]
    ## prob 2: [m]
    ## prob 4: [m1, m1, m2, m2, m3, m3]
#m = np.array([1., 1.]) # for prob 1
#m = np.array([2.]) # for prob 2
m = np.array([1.e-7, 1.e-7, 1.e-7, 1.e-7, 1.e-7, 1.e-7]) # for prob 4
    ## If non-periodic, properties can be assigned for each element. For example:
        #m = m * np.ones((N_global,prob.DoF))
        #m[0,:] = np.array([4.e-6, 4.e-6, 2.e-6, 2.e-6, 2.e-6, 2.e-6])

### on-site damping array ###
#b = 0.0 
zeta = 0.003
freq0 = 1000
b = 2*zeta*2*np.pi*freq0*m # for prob 4

### Initial conditions ###
u0 = np.zeros((N_global, prob.DoF))
udot0 = np.zeros((N_global, prob.DoF))
#udot0[0,0] = 1.0

##### END of user inputs #####

########################################################
########################################################

try: 
    ig.MethodGen(dt, tStart, tEnd, dtWrite, sol_flag, sol_params)
except NameError:
    ig.MethodGen(dt, tStart, tEnd, dtWrite, sol_flag)
    print(f" >> Parameters for {sol_flag.upper()} method are not provided. Default parameters will be used.")

if prob.flag == 1:
    ig.DesignGen(sol_flag,des,N_global,prob,BC,k,m,b,LC,LC_val,u0,udot0,L,G)
elif prob.flag == 2:
    ig.DesignGen(sol_flag,des,N_global,prob,BC,k,m,b,LC,LC_val,u0,udot0)
else:
    ig.DesignGen(sol_flag,des,N_global,prob,BC,k,m,b,LC,LC_val,u0,udot0,L)
