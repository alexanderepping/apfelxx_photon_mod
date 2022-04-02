import numpy as np
import matplotlib.pyplot as plt

mode = "data"

# import the data
particle = np.genfromtxt('data_Evolution.txt',dtype='str', max_rows=1)
ylabel = "x * " + str(particle) + " / alpha  PDF"
x, apfel, lhapdf, ratio = np.loadtxt(open("data_Evolution.txt"), delimiter="	", unpack=True, skiprows=1)





if (mode == "data"):
    # delete the zeros out of the apfel values
    # apfel2 = np.array([])
    # lhapdf2 = np.array([])
    # x2 = np.array([])
    # for i in range(len(apfel)):
    #     if apfel[i]!=0:
    #         apfel2 = np.append(apfel2, apfel[i])
    #         lhapdf2 = np.append(lhapdf2, lhapdf[i])
    #         x2 = np.append(x2, x[i])
    
    plt.plot(x, apfel, label="apfel values")
    plt.plot(x, lhapdf, label="lhapdf values")
    #plt.plot(x2, apfel2, label="apfel values")
    #plt.plot(x2, lhapdf2, label="lhapdf values")
    
    plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel("x")
    plt.ylabel(ylabel)
    plt.legend(loc="upper right")
    plt.xlim(right=1, left=10**(-4))
    plt.show()

if (mode == "ratio"):
    # delete the zeros out of the ratio values
    # ratio2 = np.array([])
    # x2 = np.array([])
    # for i in range(len(ratio)):
    #     if ratio[i]!=0:
    #         ratio2 = np.append(ratio2, ratio[i])
    #         x2 = np.append(x2, x[i])
    # plt.plot(x2,ratio2, label="ratio of apfel/lhapdf")
    # plt.plot(x2,x2/x2, label="ratio = 1")
    plt.plot(x,ratio, label="ratio of apfel/lhapdf")
    plt.plot(x,x/x, label="ratio = 1")

    #plt.xscale('log')
    plt.xlabel("x")
    plt.ylabel(ylabel)
    plt.legend(loc="upper left")
    plt.show()

