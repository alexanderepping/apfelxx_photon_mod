import numpy as np
import matplotlib.pyplot as plt

# import the data
particle = np.genfromtxt('data_Evolution.txt',dtype='str', max_rows=1)
ylabel = str(particle) + " PDF"
x, apfel, lhapdf = np.loadtxt(open("data_Evolution.txt"), delimiter="	", unpack=True, skiprows=1)

# delete the zeros out of the 
apfel2 = np.array([])
x2 = np.array([])
for i in range(len(apfel)):
    if apfel[i]!=0:
        apfel2 = np.append(apfel2, apfel[i])
        x2 = np.append(x2, x[i])

plt.plot(x2, apfel2, label="apfel values")
plt.plot(x, lhapdf, label="lhapdf values")

plt.xscale('log')

plt.xlabel("x")
plt.ylabel(ylabel)
plt.legend(loc="upper right")
#plt.ylim(top=1, bottom=0)
#plt.xlim(right=1, left=10**(-5))
plt.show()
