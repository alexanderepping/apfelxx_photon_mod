import numpy as np
import matplotlib.pyplot as plt

# import the data
x_vadim4, lhapdf_vadim4, = np.loadtxt(open("/home/alexander/Documents/apfelxx_photon_mod/plottingPython/GRV_param_Q2_4_vadim.dat"), delimiter=";", unpack=True)
x_vadim100, lhapdf_vadim100, = np.loadtxt(open("/home/alexander/Documents/apfelxx_photon_mod/plottingPython/GRV_param_Q2_100_vadim.dat"), delimiter=";", unpack=True)
x_alex4, apfel_alex4, lhapdf_alex4, ratio_alex4 = np.loadtxt(open("/home/alexander/Documents/apfelxx_photon_mod/plottingPython/GRV_param_Q2_4_alex.txt"), delimiter="	", unpack=True)
x_alex100, apfel_alex100, lhapdf_alex100, ratio_alex100 = np.loadtxt(open("/home/alexander/Documents/apfelxx_photon_mod/plottingPython/GRV_param_Q2_100_alex.txt"), delimiter="	", unpack=True)

# delete the zeros out of the apfel values
apfel2_alex4 = np.array([])
x2_alex4 = np.array([])
for i in range(len(apfel_alex4)):
    if apfel_alex4[i]!=0:
        apfel2_alex4 = np.append(apfel2_alex4, apfel_alex4[i])
        x2_alex4 = np.append(x2_alex4, x_alex4[i])


# delete the zeros out of the apfel values
apfel2_alex100 = np.array([])
x2_alex100 = np.array([])
for i in range(len(apfel_alex100)):
    if apfel_alex100[i]!=0:
        apfel2_alex100 = np.append(apfel2_alex100, apfel_alex100[i])
        x2_alex100 = np.append(x2_alex100, x_alex100[i])

    
plt.plot(x_vadim4, lhapdf_vadim4, label="vadim4")
plt.plot(x_alex4, lhapdf_alex4, label="lhapdf_alex4")
plt.plot(x2_alex4, apfel2_alex4, label="apfel_alex4")
plt.xscale('log')
plt.xlabel("x")
plt.legend(loc="upper right")
plt.ylim(top=1.5*np.amax(lhapdf_vadim4), bottom=0)
plt.xlim(right=1, left=x_vadim4[0])
plt.show()

    
plt.plot(x_vadim100, lhapdf_vadim100, label="vadim100")
plt.plot(x_alex100, lhapdf_alex100, label="lhapdf_alex100")
plt.plot(x2_alex100, apfel2_alex100, label="apfel_alex100")
plt.xscale('log')
plt.xlabel("x")
plt.legend(loc="upper right")
plt.ylim(top=1.5*np.amax(lhapdf_vadim100), bottom=0)
plt.xlim(right=1, left=x_vadim100[0])
plt.show()



