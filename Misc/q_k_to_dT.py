import numpy as np
from scipy import optimize

def q_minimise_K(q, LMTD_out = False):
    K = 90
    A = 1.4

    C_w = 0.6944 * 4200

    C_a = 0.197 * 1005

    T_a_in = 70
    T_w_in = 24.5
    T_a_out = T_a_in - q/C_a
    T_w_out = T_w_in + q/C_w
    dT1 = T_a_in - T_w_out
    dT2 = T_a_out - T_w_in
    if dT2 < 0 or dT1 < 0:
        return(1e10)
    
    LMTD = (dT2 - dT1) / np.log(dT2/dT1)
    if LMTD_out:
        print(LMTD)

    q_calc = A * K * LMTD
    return((q_calc - q)**2)

#Q = optimize.minimize_scalar(q_minimise_K).x
#print(Q)
#q_minimise_K(Q, LMTD_out=True)

def LMTD_q(q):
    C_w = 0.264 * 4200

    C_a = 0.15092175* 1005

    T_a_in = 20
    T_w_in = 40
    T_a_out = T_a_in - q/C_a
    T_w_out = T_w_in + q/C_w
    dT1 = T_a_in - T_w_out
    dT2 = T_a_out - T_w_in

    
    LMTD = (dT2 - dT1) / np.log(dT2/dT1)
    print(LMTD)

LMTD_q(2210)