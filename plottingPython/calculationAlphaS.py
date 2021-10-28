#calculation of alpha s
import numpy as np

qref = 91.1876

def alphas(Q, nf, order):
    beta_0 = 11 - 2*nf/3
    beta_1 = 102 - 38*nf/3
    lambda_arr = {0:{3:0.232, 4:0.2, 5:0.153, 6:0.082},
                  1:{3:0.248, 4:0.2, 5:0.131, 6:0.053}}

    lnQ2L2 = np.log(Q**2/(lambda_arr[order][nf]**2))
    
    if (order == 0):
        return 4*np.pi / (beta_0*lnQ2L2)
    if (order == 1):
        return 4*np.pi * (1/(beta_0*lnQ2L2) - beta_1/(beta_0**3) * np.log(lnQ2L2)/(lnQ2L2**2))


    
print(alphas(qref, 5, 1))