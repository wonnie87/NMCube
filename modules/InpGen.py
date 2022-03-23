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
        elif sol_flag.upper() == 'RK2':
            file.write("**c2\n")
            file.write(f"{params['c2']}\n")
            file.write("**a21\n")
            file.write(f"{params['a21']}\n")
            file.write("**b1, b2\n")
            file.write(f"{params['b1']}, {params['b2']}\n")
        elif sol_flag.upper() == 'RK3':
            file.write("**c2, c3\n")
            file.write(f"{params['c2']}, {params['c3']}\n")
            file.write("**a21, a31, a32\n")
            file.write(f"{params['a21']}, {params['a31']}, {params['a32']}\n")
            file.write("**b1, b2, b3\n")
            file.write(f"{params['b1']}, {params['b2']}, {params['b3']}\n")
        elif sol_flag.upper() == 'RK4':
            file.write("**c2, c3, c4\n")
            file.write(f"{params['c2']}, {params['c3']}, {params['c4']}\n")
            file.write("**a21, a31, a32, a41, a42, a43\n")
            file.write(f"{params['a21']}, {params['a31']}, {params['a32']}, \
{params['a41']}, {params['a42']}, {params['a43']}\n")
            file.write("**b1, b2, b3, b4\n")
            file.write(f"{params['b1']}, {params['b2']}, {params['b3']}, {params['b4']}\n")
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
        elif sol_flag.upper() == 'RK2':
            file.write("**c2\n")
            file.write(f"{1.0}\n")
            file.write("**a21\n")
            file.write(f"{1.0}\n")
            file.write("**b1, b2\n")
            file.write(f"{1./2}, {1./2}\n")
        elif sol_flag.upper() == 'RK3':
            file.write("**c2, c3\n")
            file.write(f"{0.5}, {1.0}\n")
            file.write("**a21, a31, a32\n")
            file.write(f"{1.0}, {-1.0}, {2.0}\n")
            file.write("**b1, b2, b3\n")
            file.write(f"{1./6}, {2./3}, {1./6}\n")
        elif sol_flag.upper() == 'RK4':
            file.write("**c2, c3, c4\n")
            file.write(f"{0.5}, {0.5}, {1.0}\n")
            file.write("**a21, a31, a32, a41, a42, a43\n")
            file.write(f"{0.5}, {0.0}, {0.5}, {0.0}, {0.0}, {1.0}\n")
            file.write("**b1, b2, b3, b4\n")
            file.write(f"{1./6}, {1./3}, {1./3}, {1./6}\n")
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
def DesignGen(sol_flag,des,N_global,prob,BC,s,m,b,LC,LC_val,u0,udot0):
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
    file.write("{}\n".format(N_global))
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

    for n in range(len(s)):
        if n == len(s)-1:
            file.write("{}\n".format(s[n]))
        else:
            file.write("{}, ".format(s[n]))

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
        m_global = m.copy()
    elif m.shape == (prob.DoF,):
        m_global = m * np.ones((N_global,prob.DoF))
    else:
        sys.exit(f" >> Shape of mass array is incompatible. It needs to be either ({N_global},{prob.DoF}) or ({prob.DoF},).")

    for i in range(N_global):
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

    if b.shape == (N_global, prob.DoF):
        b_global = b.copy()
    elif b.shape == (prob.DoF,):
        b_global = b * np.ones((N_global,prob.DoF))
    else:
        sys.exit(f" >> Shape of on-site damping array is incompatible. It needs to be either ({N_global},{prob.DoF}) or ({prob.DoF},).")

    for i in range(N_global):
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
    for i in range(N_global):
        for j in range(prob.DoF):
            if j == prob.DoF-1:
                file.write(f"{u0[i,j]}\n")
            else:
                file.write(f"{u0[i,j]}, ")

    file.write("**udot0\n")
    for i in range(N_global):
        for j in range(prob.DoF):
            if j == prob.DoF-1:
                file.write(f"{udot0[i,j]}\n")
            else:
                file.write(f"{udot0[i,j]}, ")

    file.close()
