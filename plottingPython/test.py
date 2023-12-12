import numpy as np
import matplotlib.pyplot as plt
import os
from experimentalData import *





DataSetsLess  = [ALEPH1_DATA, ALEPH2_DATA, AMY_DATA, DELPHI_DATA, JADE_DATA, L3_DATA, OPAL1_DATA, OPAL2_less_DATA, PLUTO_DATA, TASSO_DATA, TOPAZ_DATA]
NamesSetsLess = ["ALEPH1", "ALEPH2", "AMY", "DELPHI", "JADE", "L3", "OPAL1", "OPAL2_less", "PLUTO", "TASSO", "TOPAZ"]

#DataSetsOther = [ALEPH1_DATA, ALEPH2_DATA, AMY_DATA, DELPHI_DATA, JADE_DATA, L3_DATA, OPAL1_DATA, OPAL2_DATA, PLUTO_DATA, TASSO_DATA, TOPAZ_DATA]


n = 0
for i in range(len(DataSetsLess)):
    DataSet = DataSetsLess[i]
    nTemp = len(DataSet["Q2Data"])

    n = n + nTemp
    print(str(nTemp) + "  - " + NamesSetsLess[i])

print("\nTotal: " + str(n))