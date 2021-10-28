###########################
# imports
###########################
import numpy as np



###########################
# definitions of exp data  
###########################

# ALEPH (Nisius, Appendix D, Table D.1, p. 177)
ALEPH_DATA = {"F2Gamma"     : [ np.array([0.30, 0.40, 0.41, 0.27]),             # Q2 = 9.9
                                np.array([0.36, 0.34, 0.56, 0.45]),             # Q2 = 20.7
                                np.array([0.65, 0.70, 1.28])],                  # Q2 = 284

              "y_error"     : [ np.array([0.03, 0.07, 0.10, 0.16]),             # Q2 = 9.9
                                np.array([0.05, 0.12, 0.11, 0.12]),             # Q2 = 20.7
                                np.array([0.14, 0.25, 0.37])],                  # Q2 = 284

              "x_data"      : [ np.array([]),                                   # Q2 = 9.9
                                np.array([]),                                   # Q2 = 20.7
                                np.array([])],                                  # Q2 = 284

              "x_error"     : [ np.array([]),                                   # Q2 = 9.9
                                np.array([]),                                   # Q2 = 20.7
                                np.array([])],                                  # Q2 = 284

              "intervals"   : [ np.array([0.005, 0.080, 0.200, 0.400, 0.800]),  # Q2 = 9.9
                                np.array([0.009, 0.120, 0.270, 0.500, 0.890]),  # Q2 = 20.7
                                np.array([0.030, 0.350, 0.650, 0.970])]}        # Q2 = 284

# AMY (Nisius, Appendix D, Table D.2, p. 178)
AMY_DATA =   {"F2Gamma"     : [ np.array([0.337, 0.302, 0.322]),            # Q2 = 6.8
                                np.array([0.650, 0.600, 0.650]),            # Q2 = 73
                                np.array([0.94, 0.82])],                    # Q2 = 390

              "y_error"     : [ np.array([0.053, 0.049, 0.097]),            # Q2 = 6.8
                                np.array([0.100, 0.160, 0.140]),            # Q2 = 73
                                np.array([0.250, 0.190])],                  # Q2 = 390

              "x_data"      : [ np.array([]),                               # Q2 = 6.8
                                np.array([]),                               # Q2 = 73
                                np.array([])],                              # Q2 = 390

              "x_error"     : [ np.array([]),                               # Q2 = 6.8
                                np.array([]),                               # Q2 = 73
                                np.array([])],                              # Q2 = 390

              "intervals"   : [ np.array([0.015, 0.125, 0.375, 0.625]),     # Q2 = 6.8
                                np.array([0.125, 0.375, 0.625, 0.875]),     # Q2 = 73
                                np.array([0.120, 0.500, 0.800])]}           # Q2 = 390
