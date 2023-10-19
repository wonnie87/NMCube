import numpy as np

#####################################################################
class ProbType:
    """..."""

    def __init__(self, prob_flag):
        """..."""
        self.flag = prob_flag

        if self.flag == 1:
            self.note = '1D pendula chain'
            self.DoF = 1
        elif self.flag == 2:
            self.note = '1D phi4 chain'
            self.DoF = 1
        elif self.flag == 4:
            self.note = 'Bistable metabeam'
            self.DoF = 6
        elif self.flag == 45:
            self.note = 'Cont6 approximation of bistable metabeam'
            self.DoF = 6

        self.noState = 2*self.DoF


#####################################################################
def MethodGen(dt, tStart, tEnd, dtWrite, sol_flag, params=''):
    """..."""
    file = open("numMethod.inp", 'w')
    file.write("**sol_flag\n") 

    if sol_flag.upper() == 'RK1':
        file.write("10\n")
    elif sol_flag.upper() == 'RK2':
        file.write("20\n")
    elif sol_flag.upper() == 'RK3':
        file.write("30\n")
    elif sol_flag.upper() == 'RK4':
        file.write("40\n")
    elif sol_flag.upper() == 'NB':
        file.write("50\n")

    if params:
        if sol_flag.upper() == 'RK1':
            file.write("**b1\n")
            file.write(f"{params['b1']}")
            file.write("**tol_inv\n")
            file.write(f"{params['tol_inv']}\n")
            file.write("**N_inv\n")
            file.write(f"{params['N_inv']}\n")
        elif sol_flag.upper() == 'RK2':
            file.write("**c2\n")
            file.write(f"{params['c2']}\n")
            file.write("**a21\n")
            file.write(f"{params['a21']}\n")
            file.write("**b1, b2\n")
            file.write(f"{params['b1']}, {params['b2']}\n")
            file.write("**tol_inv\n")
            file.write(f"{params['tol_inv']}\n")
            file.write("**N_inv\n")
            file.write(f"{params['N_inv']}\n")
        elif sol_flag.upper() == 'RK3':
            file.write("**c2, c3\n")
            file.write(f"{params['c2']}, {params['c3']}\n")
            file.write("**a21, a31, a32\n")
            file.write(f"{params['a21']}, {params['a31']}, {params['a32']}\n")
            file.write("**b1, b2, b3\n")
            file.write(f"{params['b1']}, {params['b2']}, {params['b3']}\n")
            file.write("**tol_inv\n")
            file.write(f"{params['tol_inv']}\n")
            file.write("**N_inv\n")
            file.write(f"{params['N_inv']}\n")
        elif sol_flag.upper() == 'RK4':
            file.write("**c2, c3, c4\n")
            file.write(f"{params['c2']}, {params['c3']}, {params['c4']}\n")
            file.write("**a21, a31, a32, a41, a42, a43\n")
            file.write(f"{params['a21']}, {params['a31']}, {params['a32']}, \
{params['a41']}, {params['a42']}, {params['a43']}\n")
            file.write("**b1, b2, b3, b4\n")
            file.write(f"{params['b1']}, {params['b2']}, {params['b3']}, {params['b4']}\n")
            file.write("**tol_inv\n")
            file.write(f"{params['tol_inv']}\n")
            file.write("**N_inv\n")
            file.write(f"{params['N_inv']}\n")
        elif sol_flag.upper() == 'NB':
            file.write("**gamma_NB\n")
            file.write(f"{params['gamma_NB']}\n")
            file.write("**beta_NB\n")
            file.write(f"{params['beta_NB']}\n")
            file.write("**tol_NR\n")
            file.write(f"{params['tol_NR']}\n")
            file.write("**tol_inv\n")
            file.write(f"{params['tol_inv']}\n")
            file.write("**N_NR\n")
            file.write(f"{params['N_NR']}\n")
            file.write("**N_inv\n")
            file.write(f"{params['N_inv']}\n")

    else:
        if sol_flag.upper() == 'RK1':
            file.write("**b1\n")
            file.write(f"{1.0}\n")
            file.write("**tol_inv\n")
            file.write(f"{1.0e-16}\n")
            file.write("**N_inv\n")
            file.write(f"{200}\n")
        elif sol_flag.upper() == 'RK2':
            file.write("**c2\n")
            file.write(f"{1.0}\n")
            file.write("**a21\n")
            file.write(f"{1.0}\n")
            file.write("**b1, b2\n")
            file.write(f"{1./2}, {1./2}\n")
            file.write("**tol_inv\n")
            file.write(f"{1.0e-16}\n")
            file.write("**N_inv\n")
            file.write(f"{200}\n")
        elif sol_flag.upper() == 'RK3':
            file.write("**c2, c3\n")
            file.write(f"{0.5}, {1.0}\n")
            file.write("**a21, a31, a32\n")
            file.write(f"{1.0}, {-1.0}, {2.0}\n")
            file.write("**b1, b2, b3\n")
            file.write(f"{1./6}, {2./3}, {1./6}\n")
            file.write("**tol_inv\n")
            file.write(f"{1.0e-16}\n")
            file.write("**N_inv\n")
            file.write(f"{200}\n")
        elif sol_flag.upper() == 'RK4':
            file.write("**c2, c3, c4\n")
            file.write(f"{0.5}, {0.5}, {1.0}\n")
            file.write("**a21, a31, a32, a41, a42, a43\n")
            file.write(f"{0.5}, {0.0}, {0.5}, {0.0}, {0.0}, {1.0}\n")
            file.write("**b1, b2, b3, b4\n")
            file.write(f"{1./6}, {1./3}, {1./3}, {1./6}\n")
            file.write("**tol_inv\n")
            file.write(f"{1.0e-16}\n")
            file.write("**N_inv\n")
            file.write(f"{200}\n")
        elif sol_flag.upper() == 'NB':
            file.write("**gamma_NB\n")
            file.write(f"{0.5}\n")
            file.write("**beta_NB\n")
            file.write(f"{0.25}\n")
            file.write("**tol_NR\n")
            file.write(f"{1.0e-16}\n")
            file.write("**tol_inv\n")
            file.write(f"{1.0e-16}\n")
            file.write("**N_NR\n")
            file.write(f"{200}\n")
            file.write("**N_inv\n")
            file.write(f"{200}\n")

    file.write("**dt\n")
    file.write("{}\n".format(dt))
    file.write("**tStart\n")
    file.write("{}\n".format(tStart))
    file.write("**tEnd\n")
    file.write("{}\n".format(tEnd))
    file.write("**dtWrite\n")
    file.write("{}\n".format(dtWrite))

    file.close()


#####################################################################
def DesignGen(sol_flag,des,N_global,prob,BC,s,m,b,LC,LC_val,u0,udot0,NFE_perElem=2):
#def DesignGen(sol_flag,des,N_global,prob,BC,k,m,b,LC,LC_val,u0,udot0,L=[],G=[]):
    """..."""
    import sys

    file = open("design.inp", "w")
    file.write("**filename\n")
    if LC[0,2] == 1:
        file.write(f"P{prob.flag}_{des}_F{LC_val[0,2]}_f{LC_val[0,3]}\n")
    elif LC[0,2] == 2:
        file.write(f"P{prob.flag}_{des}_F{LC_val[0,2]}_f{LC_val[0,3]}_\
t{LC_val[0,6]-LC_val[0,5]}\n")
    elif LC[0,2] == 11 or LC[0,2] == 13:
        file.write(f"P{prob.flag}_{des}_A{LC_val[0,2]}_f{LC_val[0,3]}\n")

    file.write("**N_global\n")
    if prob.flag == 45:
        N_eq = (N_global-1)*NFE_perElem+1
    else:
        N_eq = N_global

    file.write("{}\n".format(N_eq))

    file.write("**prob_flag\n")
    file.write("{}\n".format(prob.flag))
    file.write("**BC(0), BC(1), BC(2)\n")
    for n in range(len(BC)):
        if n == len(BC)-1:
            file.write("{}\n".format(BC[n]))
        else:
            file.write("{}, ".format(BC[n]))

    if prob.flag == 1:
        file.write("**k_th, m, g, l\n")
    elif prob.flag == 2:
        file.write("**k(1), k(2), k(3), k(4), k(5)\n")
    elif prob.flag == 4:
        file.write("**k1, k2, k3, k4, k5, k6, k7, k8, L1, L2, L3, R\n")
    elif prob.flag == 45:
        file.write("**Vbar1st2, Vbar1st4, Vbar1st6, Vbar2nd11, Vbar2nd13, Vbar2nd15, \
Vbar2nd22, Vbar2nd24, Vbar2nd26, Vbar2nd33, Vbar2nd35, Vbar2nd310, \
Vbar2nd312, Vbar2nd44, Vbar2nd46, Vbar2nd49, Vbar2nd411, Vbar2nd55, \
Vbar2nd510, Vbar2nd512, Vbar2nd66, Vbar2nd69, Vbar2nd611, Vbar2nd77, \
Vbar2nd99, Vbar2nd1010, Vbar2nd1111, Vbar2nd1212, Vbar4th1111, \
Vbar4th1113, Vbar4th1115, Vbar4th1133, Vbar4th1155, Vbar4th1333, \
Vbar4th1555, Vbar4th3333, Vbar4th3335, Vbar4th33310, Vbar4th33312, \
Vbar4th3355, Vbar4th33510, Vbar4th33512, Vbar4th331010, \
Vbar4th331212, Vbar4th3555, Vbar4th35510, Vbar4th35512, \
Vbar4th351010, Vbar4th351212, Vbar4th3101010, Vbar4th3121212, \
Vbar4th5555, Vbar4th55510, Vbar4th55512, Vbar4th551010, \
Vbar4th551212, Vbar4th5101010, Vbar4th5121212, Vbar4th10101010, \
Vbar4th12121212, L1, h\n")

    if prob.flag == 45:
        file.write("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n".format(\
(s[1]*(-s[9] + np.sqrt(s[9]**2 + s[11]**2)) + s[2]*(s[10] - np.sqrt(s[10]**2 + s[11]**2)))/s[8],\
(s[1]*(s[9] - np.sqrt(s[9]**2 + s[11]**2)))/s[8],\
(s[2]*(-s[10] + np.sqrt(s[10]**2 + s[11]**2)))/s[8],\
(s[1]*s[10]*(s[9] - np.sqrt(s[9]**2 + s[11]**2)) + s[2]*s[9]*(s[10] - np.sqrt(s[10]**2 + s[11]**2)))/(s[8]*s[9]*s[10]),\
(s[1]*(-s[9] + np.sqrt(s[9]**2 + s[11]**2)))/(s[8]*s[9]),\
(s[2]*(-s[10] + np.sqrt(s[10]**2 + s[11]**2)))/(s[8]*s[10]),\
(s[1] + s[2])/s[8],\
-(s[1]/s[8]),\
-(s[2]/s[8]),\
(s[1] + (2*s[6]*s[8]**2)/(s[8]**2 + (s[9] + s[10])**2) - (s[1]*np.sqrt(s[9]**2 + s[11]**2))/s[9])/s[8],\
(-2*s[6]*s[8])/(s[8]**2 + (s[9] + s[10])**2),\
(s[6]*s[8]*(s[9] + s[10]))/(s[8]**2 + (s[9] + s[10])**2),\
(s[6]*s[8]*(s[9] + s[10]))/(s[8]**2 + (s[9] + s[10])**2),\
(s[1] + s[5] + (2*s[6]*(s[9] + s[10])**2)/(s[8]**2 + (s[9] + s[10])**2))/s[8],\
-((s[5] + (2*s[6]*(s[9] + s[10])**2)/(s[8]**2 + (s[9] + s[10])**2))/s[8]),\
(s[6]*s[8]*(s[9] + s[10]))/(s[8]**2 + (s[9] + s[10])**2),\
(s[6]*s[8]*(s[9] + s[10]))/(s[8]**2 + (s[9] + s[10])**2),\
(s[2] + (2*s[6]*s[8]**2)/(s[8]**2 + (s[9] + s[10])**2) - (s[2]*np.sqrt(s[10]**2 + s[11]**2))/s[10])/s[8],\
-((s[6]*s[8]*(s[9] + s[10]))/(s[8]**2 + (s[9] + s[10])**2)),\
-((s[6]*s[8]*(s[9] + s[10]))/(s[8]**2 + (s[9] + s[10])**2)),\
(s[2] + s[5] + (2*s[6]*(s[9] + s[10])**2)/(s[8]**2 + (s[9] + s[10])**2))/s[8],\
-((s[6]*s[8]*(s[9] + s[10]))/(s[8]**2 + (s[9] + s[10])**2)),\
-((s[6]*s[8]*(s[9] + s[10]))/(s[8]**2 + (s[9] + s[10])**2)),\
s[0]*s[8],\
s[3]*s[8] + (s[6]*s[8]**3)/(s[8]**2 + (s[9] + s[10])**2),\
(s[6]*s[8]*(s[9] + s[10])**2)/(s[8]**2 + (s[9] + s[10])**2),\
s[4]*s[8] + (s[6]*s[8]**3)/(s[8]**2 + (s[9] + s[10])**2),\
(s[6]*s[8]*(s[9] + s[10])**2)/(s[8]**2 + (s[9] + s[10])**2),\
(3*((s[1]*np.sqrt(s[9]**2 + s[11]**2))/s[9]**3 + (s[2]*np.sqrt(s[10]**2 + s[11]**2))/s[10]**3))/s[8],\
(-3*s[1]*np.sqrt(s[9]**2 + s[11]**2))/(s[8]*s[9]**3),\
(-3*s[2]*np.sqrt(s[10]**2 + s[11]**2))/(s[8]*s[10]**3),\
(3*s[1]*np.sqrt(s[9]**2 + s[11]**2))/(s[8]*s[9]**3),\
(3*s[2]*np.sqrt(s[10]**2 + s[11]**2))/(s[8]*s[10]**3),\
(-3*s[1]*np.sqrt(s[9]**2 + s[11]**2))/(s[8]*s[9]**3),\
(-3*s[2]*np.sqrt(s[10]**2 + s[11]**2))/(s[8]*s[10]**3),\
(3*(-8*s[6]*s[8]**2*s[9]**5*(s[9] + s[10])**2 + 2*s[6]*s[9]**7*(s[9] + s[10])**2 - 16*s[6]*s[8]**2*s[9]**4*s[10]*(s[9] + s[10])**2 + 8*s[6]*s[9]**6*s[10]*(s[9] + s[10])**2 - 8*s[6]*s[8]**2*s[9]**3*s[10]**2*(s[9] + s[10])**2 + 12*s[6]*s[9]**5*s[10]**2*(s[9] + s[10])**2 + 8*s[6]*s[9]**4*s[10]**3*(s[9] + s[10])**2 + 2*s[6]*s[9]**3*s[10]**4*(s[9] + s[10])**2 + s[5]*s[9]**3*(s[8]**2 + (s[9] + s[10])**2)**3 + s[1]*(s[9] + s[10])**2*(s[8]**2 + (s[9] + s[10])**2)**3*np.sqrt(s[9]**2 + s[11]**2)))/(s[8]*s[9]**3*(s[9] + s[10])**2*(s[8]**2 + (s[9] + s[10])**2)**3),\
(3*(-s[5] - (2*s[6]*(-4*s[8]**2*(s[9] + s[10])**4 + (s[9] + s[10])**6))/(s[8]**2 + (s[9] + s[10])**2)**3))/(s[8]*(s[9] + s[10])**2),\
(3*s[6]*s[8]*(s[9] + s[10])*(2*s[8]**2 - 3*(s[9] + s[10])**2))/(s[8]**2 + (s[9] + s[10])**2)**3,\
(3*s[6]*s[8]*(s[9] + s[10])*(2*s[8]**2 - 3*(s[9] + s[10])**2))/(s[8]**2 + (s[9] + s[10])**2)**3,\
(3*(s[5] + (2*s[6]*(-4*s[8]**2*(s[9] + s[10])**4 + (s[9] + s[10])**6))/(s[8]**2 + (s[9] + s[10])**2)**3))/(s[8]*(s[9] + s[10])**2),\
(3*s[6]*s[8]*(s[9] + s[10])*(-2*s[8]**2 + 3*(s[9] + s[10])**2))/(s[8]**2 + (s[9] + s[10])**2)**3,\
(3*s[6]*s[8]*(s[9] + s[10])*(-2*s[8]**2 + 3*(s[9] + s[10])**2))/(s[8]**2 + (s[9] + s[10])**2)**3,\
-((s[6]*s[8]*(2*s[8]**4 - 11*s[8]**2*(s[9] + s[10])**2 + 2*(s[9] + s[10])**4))/(s[8]**2 + (s[9] + s[10])**2)**3),\
-((s[6]*s[8]*(2*s[8]**4 - 11*s[8]**2*(s[9] + s[10])**2 + 2*(s[9] + s[10])**4))/(s[8]**2 + (s[9] + s[10])**2)**3),\
(3*(-s[5] - (2*s[6]*(-4*s[8]**2*(s[9] + s[10])**4 + (s[9] + s[10])**6))/(s[8]**2 + (s[9] + s[10])**2)**3))/(s[8]*(s[9] + s[10])**2),\
(3*s[6]*s[8]*(s[9] + s[10])*(2*s[8]**2 - 3*(s[9] + s[10])**2))/(s[8]**2 + (s[9] + s[10])**2)**3,\
(3*s[6]*s[8]*(s[9] + s[10])*(2*s[8]**2 - 3*(s[9] + s[10])**2))/(s[8]**2 + (s[9] + s[10])**2)**3,\
(s[6]*s[8]*(2*s[8]**4 - 11*s[8]**2*(s[9] + s[10])**2 + 2*(s[9] + s[10])**4))/(s[8]**2 + (s[9] + s[10])**2)**3,\
(s[6]*s[8]*(2*s[8]**4 - 11*s[8]**2*(s[9] + s[10])**2 + 2*(s[9] + s[10])**4))/(s[8]**2 + (s[9] + s[10])**2)**3,\
(3*s[6]*s[8]**3*(s[9] + s[10])*(-3*s[8]**2 + 2*(s[9] + s[10])**2))/(s[8]**2 + (s[9] + s[10])**2)**3,\
(3*s[6]*s[8]**3*(s[9] + s[10])*(-3*s[8]**2 + 2*(s[9] + s[10])**2))/(s[8]**2 + (s[9] + s[10])**2)**3,\
(3*(-8*s[6]*s[8]**2*s[9]**2*s[10]**3*(s[9] + s[10])**2 + 2*s[6]*s[9]**4*s[10]**3*(s[9] + s[10])**2 - 16*s[6]*s[8]**2*s[9]*s[10]**4*(s[9] + s[10])**2 + 8*s[6]*s[9]**3*s[10]**4*(s[9] + s[10])**2 - 8*s[6]*s[8]**2*s[10]**5*(s[9] + s[10])**2 + 12*s[6]*s[9]**2*s[10]**5*(s[9] + s[10])**2 + 8*s[6]*s[9]*s[10]**6*(s[9] + s[10])**2 + 2*s[6]*s[10]**7*(s[9] + s[10])**2 + s[5]*s[10]**3*(s[8]**2 + (s[9] + s[10])**2)**3 + s[2]*(s[9] + s[10])**2*(s[8]**2 + (s[9] + s[10])**2)**3*np.sqrt(s[10]**2 + s[11]**2)))/(s[8]*s[10]**3*(s[9] + s[10])**2*(s[8]**2 + (s[9] + s[10])**2)**3),\
(3*s[6]*s[8]*(s[9] + s[10])*(-2*s[8]**2 + 3*(s[9] + s[10])**2))/(s[8]**2 + (s[9] + s[10])**2)**3,\
(3*s[6]*s[8]*(s[9] + s[10])*(-2*s[8]**2 + 3*(s[9] + s[10])**2))/(s[8]**2 + (s[9] + s[10])**2)**3,\
-((s[6]*s[8]*(2*s[8]**4 - 11*s[8]**2*(s[9] + s[10])**2 + 2*(s[9] + s[10])**4))/(s[8]**2 + (s[9] + s[10])**2)**3),\
-((s[6]*s[8]*(2*s[8]**4 - 11*s[8]**2*(s[9] + s[10])**2 + 2*(s[9] + s[10])**4))/(s[8]**2 + (s[9] + s[10])**2)**3),\
(3*s[6]*s[8]**3*(s[9] + s[10])*(3*s[8]**2 - 2*(s[9] + s[10])**2))/(s[8]**2 + (s[9] + s[10])**2)**3,\
(3*s[6]*s[8]**3*(s[9] + s[10])*(3*s[8]**2 - 2*(s[9] + s[10])**2))/(s[8]**2 + (s[9] + s[10])**2)**3,\
(3*s[8]*(s[6]*s[8]**6 - 4*s[6]*s[8]**4*s[9]**2 - 8*s[6]*s[8]**4*s[9]*s[10] - 4*s[6]*s[8]**4*s[10]**2 + s[3]*(s[8]**2 + (s[9] + s[10])**2)**3))/(s[8]**2 + (s[9] + s[10])**2)**3,\
(3*s[8]*(s[6]*s[8]**6 - 4*s[6]*s[8]**4*s[9]**2 - 8*s[6]*s[8]**4*s[9]*s[10] - 4*s[6]*s[8]**4*s[10]**2 + s[4]*(s[8]**2 + (s[9] + s[10])**2)**3))/(s[8]**2 + (s[9] + s[10])**2)**3,\
s[8],\
s[8]/NFE_perElem \
))

    else:
        for n in range(len(s)):
            if n == len(s)-1:
                file.write("{}\n".format(s[n]))
            else:
                file.write("{}, ".format(s[n]))

# The following force arrays need to be updated for FE problem
    LC_dim2, LC_dim1 = LC.shape
    file.write("**LC_dim2\n")
    file.write("{}\n".format(LC_dim2))
    file.write("**Load cases\n")
    for i in range(LC_dim2):
        for j in range(LC_dim1):
            if j == LC_dim1-1:
                file.write("{}\n".format(LC[i,j]))
            else:
                file.write("{}, ".format(LC[i,j]))

    file.write("**Corresponding LC parameters\n")
    for i in range(LC_dim2):
        for j in range(LC_val.shape[1]):
            if j == LC_val.shape[1]-1:
                file.write("{}\n".format(LC_val[i,j]))
            else:
                file.write("{}, ".format(LC_val[i,j]))

    file.write("**DoF\n")
    file.write("{}\n".format(prob.DoF))

    if prob.flag == 1:
        file.write("**I\n")
    elif prob.flag == 2:
        file.write("**m(1)\n")
    elif prob.flag == 4:
        file.write("**m1, m1, m2, m2, m3, m3\n")
    elif prob.flag == 45:
        file.write("**Tbar2nd11, Tbar2nd22, Tbar2nd33, Tbar2nd44, Tbar2nd55, Tbar2nd66, \
Tbar2nd77, Tbar2nd88, Tbar2nd99, Tbar2nd1010, Tbar2nd1111, Tbar2nd1212\n")
 
#    if prob.flag == 1:
#        if m.shape == (N_global, 2*prob.DoF):
#            m_global = m.copy()
#        elif m.shape == (2*prob.DoF,):
#            m_global = m * np.ones((N_global,2*prob.DoF))
#        else:
#            sys.exit(f" >> Shape of mass array is incompatible. It needs to be either ({N_global},2x{prob.DoF}) or (2x{prob.DoF},).")
#
#    else:
#        if m.shape == (N_global, prob.DoF):
#            m_global = m.copy()
#        elif m.shape == (prob.DoF,):
#            m_global = m * np.ones((N_global,prob.DoF))
#        else:
#            sys.exit(f" >> Shape of mass array is incompatible. It needs to be either ({N_global},{prob.DoF}) or ({prob.DoF},).")

    if m.shape == (N_global, prob.DoF):
        if prob.flag == 45:
            # Do something (e.g., linear interpolation)
            pass
        else:
            m_global = m.copy()
    elif m.shape == (prob.DoF,):
        if prob.flag == 45:
            m_global = np.array([m[0]/s[8], m[0]/s[8], m[2]/s[8], m[2]/s[8], m[4]/s[8], m[4]/s[8],\
(m[0]*s[8])/12., (m[0]*s[8])/12., (m[2]*s[8])/12., (m[2]*s[8])/12., (m[4]*s[8])/12., (m[4]*s[8])/12.])\
* np.ones((N_eq, 12))
        else:
            m_global = m * np.ones((N_eq,prob.DoF))
    else:
        sys.exit(f" >> Shape of mass array is incompatible. It needs to be either ({N_global},{prob.DoF}) or ({prob.DoF},).")

    for i in range(N_eq):
        for j in range(m_global.shape[1]):
            if j == m_global.shape[1]-1:
                file.write("{}\n".format(m_global[i,j]))
            else:
                file.write("{}, ".format(m_global[i,j]))

    if prob.flag == 1:
        file.write("**b(1)\n")
    elif prob.flag == 2:
        file.write("**b(1)\n")
    elif prob.flag == 4:
        file.write("**b(1), b(2), b(3), b(4), b(5), b(6)\n")
    elif prob.flag == 45:
        file.write("**b(1)\n")

    if b.shape == (N_global, prob.DoF):
        if prob.flag == 45:
            # Do something
            pass
        else:
            b_global = b.copy()
    elif b.shape == (prob.DoF,):
        if prob.flag == 45:
            b_global = b[1]/m[1] * np.ones((N_eq, 1))
        else:
            b_global = b * np.ones((N_eq,prob.DoF))
    else:
        sys.exit(f" >> Shape of on-site damping array is incompatible. It needs to be either ({N_global},{prob.DoF}) or ({prob.DoF},).")

    for i in range(N_eq):
        for j in range(b_global.shape[1]):
            if j == b_global.shape[1]-1:
                file.write("{}\n".format(b_global[i,j]))
            else:
                file.write("{}, ".format(b_global[i,j]))

#    if sol_flag.upper() in ['RK1', 'RK2', 'RK3', 'RK4']:
#        file.write("**x0\n")
#        for i in range(N_global):
#            for j in range(prob.DoF):
#                if j == prob.DoF-1:
#                    file.write("{}, ".format(u0[i,j]))
#                    file.write("{}\n".format(udot0[i,j]))
#                else:
#                    file.write("{}, ".format(u0[i,j]))
#                    file.write("{}, ".format(udot0[i,j]))
#
#    elif sol_flag.upper() in ['NB']:
#        file.write("**u0\n")
#        for i in range(N_global):
#            for j in range(prob.DoF):
#                if j == prob.DoF-1:
#                    file.write(f"{u0[i,j]}\n")
#                else:
#                    file.write(f"{u0[i,j]}, ")
#
#        file.write("**udot0\n")
#        for i in range(N_global):
#            for j in range(prob.DoF):
#                if j == prob.DoF-1:
#                    file.write(f"{udot0[i,j]}\n")
#                else:
#                    file.write(f"{udot0[i,j]}, ")

    file.write("**u0\n")
    for i in range(N_eq):
        for j in range(prob.DoF):
            #### Modify s.t. values at mid points linearly interpolate.
            if i == N_eq-1:
                u0_val = u0[i//NFE_perElem,j]
            else:
                u0_val = u0[i//NFE_perElem,j] + (i%NFE_perElem)/NFE_perElem * (u0[i//NFE_perElem+1,j]-u0[i//NFE_perElem,j])

            if j == prob.DoF-1:
                file.write(f"{u0_val}\n")
            else:
                file.write(f"{u0_val}, ")

    file.write("**udot0\n")
    for i in range(N_eq):
        for j in range(prob.DoF):
            #### Modify s.t. values at mid points linearly interpolate.
            if i == N_eq-1:
                udot0_val = udot0[i//NFE_perElem,j]
            else:
                udot0_val = udot0[i//NFE_perElem,j] + (i%NFE_perElem)/NFE_perElem * (udot0[i//NFE_perElem+1,j]-udot0[i//NFE_perElem,j])

            if j == prob.DoF-1:
                file.write(f"{udot0_val}\n")
            else:
                file.write(f"{udot0_val}, ")

    file.close()
