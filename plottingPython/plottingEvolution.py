import numpy as np
import matplotlib.pyplot as plt

# import the data
x, apfel, lhapdf = np.loadtxt(open("data_Evolution.txt"), delimiter="	", unpack=True)

# delete the zeros out of the 
apfel2 = np.array([])
x2 = np.array([])
for i in range(len(apfel)):
    if apfel[i]!=0:
        apfel2 = np.append(apfel2, apfel[i])
        x2 = np.append(x2, x[i])
    
plt.plot(x2, apfel2, label="apfel values")
plt.plot(x, lhapdf, label="lhapdf values")

plt.xlabel("x")
plt.legend(loc="upper right")
plt.ylim(top=1, bottom=0)
plt.show()