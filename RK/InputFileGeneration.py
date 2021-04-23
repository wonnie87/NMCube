import numpy as np

### Input parameters for numerical method ###
sol_flag = 40
if sol_flag == 10:
    b1=1.0
elif sol_flag == 20:
    c2=1.0
    a21=1.0
    b1=1./2; b2=1./2
elif sol_flag == 30:
    c2=0.5; c3=1.0
    a21=0.5; a31=-1.0; a32=2.0
    b1=1./6; b2=2./3; b3=1./6
elif sol_flag == 40:
    c2=0.5; c3=0.5; c4=1.0
    a21=0.5; a31=0.0; a32=0.5; a41=0.0; a42=0.0; a43=1.0
    b1=1./6; b2=1./3; b3=1./3; b4=1./6

dt = 1.0e-4
tStart = 0.0
tEnd = 10.0
dtWrite = 0.1

### Input parameters for system design ###
N_global = 4 # Number of unit cells
prob_flag = 1 # 1: 1D pendula chain, 2: 1D phi4 chain, 3: N/A, 4: bistable metabeam
BC = np.array([100, 110, 120])
G = 1. # for prob 1
L = np.array([1.]) # for prob 1
#L = np.array([20., 40., 20., 8.]) # for prob 4
k = np.array([1.]) # for prob 1
#k = np.array([1.0, 0.0, -0.030, 0.0, 0.015]) # for prob 2
#k = np.array([1.241, 1.0758, 0.6, 1500., 1500., 1500., 1500., 1500.]) # for prob 4
LC = np.array([[1,1,1]]) # [[UC, direction, Loadcase], ... ]
amp = 0.0
freqIn = 0.0
phi0 = 0.0
t1 = 0.0
t2 = 0.0
## LC_val = [tForceStart, tForceEnd, amp, freqIn, phi0, c1, c2] for LC_case 1
LC_val = np.array([[t1, t2, amp, freqIn, phi0, t1, t2]])
LC_dim2, LC_dim1 = LC.shape
DoF = 1
noState = 2*DoF
m = np.array([1., 1.]) * np.ones((N_global,DoF)) # for prob 1
#m = 1. * np.ones((N_global,DoF)) # for prob 2
#m = np.array([2.e-6, 2.e-6, 1.e-6, 1.e-6, 1.e-6, 1.e-6]) * np.ones((N_global,DoF)) # for prob 4
b = 0.0 * np.ones((N_global,DoF))
x = np.zeros((N_global,noState))
x[0,1] = 2.0

outF1 = open("numMethod.inp", "w")
outF2 = open("design.inp", "w")

outF1.write("**sol_flag\n")
outF1.write("{}\n".format(sol_flag))
if sol_flag == 10:
    outF1.write("**b1\n")
    outF1.write("{}\n".format(b1))
elif sol_flag == 20:
    outF1.write("**c2\n")
    outF1.write("{}\n".format(c2))
    outF1.write("**a21\n")
    outF1.write("{} \n".format(a21))
    outF1.write("**b1, b2\n")
    outF1.write("{}, {}\n".format(b1, b2))
elif sol_flag == 30:
    outF1.write("**c2, c3\n")
    outF1.write("{}, {}\n".format(c2, c3))
    outF1.write("**a21, a31, a32\n")
    outF1.write("{}, {}, {}\n".format(a21, a31, a32))
    outF1.write("**b1, b2, b3\n")
    outF1.write("{}, {}, {}\n".format(b1, b2, b3))
elif sol_flag == 40:
    outF1.write("**c2, c3, c4\n")
    outF1.write("{}, {}, {}\n".format(c2, c3, c4))
    outF1.write("**a21, a31, a32, a41, a42, a43\n")
    outF1.write("{}, {}, {}, {}, {}, {}\n".format(a21, a31, a32, a41, a42, a43))
    outF1.write("**b1, b2, b3, b4\n")
    outF1.write("{}, {}, {}, {}\n".format(b1, b2, b3, b4))

outF1.write("**dt\n")
outF1.write("{}\n".format(dt))
outF1.write("**tStart\n")
outF1.write("{}\n".format(tStart))
outF1.write("**tEnd\n")
outF1.write("{}\n".format(tEnd))
outF1.write("**dtWrite\n")
outF1.write("{}\n".format(dtWrite))

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
    outF2.write("**G\n")
    outF2.write("{}\n".format(G))

if prob_flag == 1:
    outF2.write("**k(1)\n")
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

outF2.write("**x0\n")
for i in range(N_global):
    for j in range(x.shape[1]):
        if j == x.shape[1]-1:
            outF2.write("{}\n".format(x[i,j]))
        else:
            outF2.write("{}, ".format(x[i,j]))

outF1.close()
outF2.close()
