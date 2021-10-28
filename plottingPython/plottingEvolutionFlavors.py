###################
# imports
###################
import numpy as np
import matplotlib.pyplot as plt
import math
import sys



###################
# definitions
###################
# variables changeable by user
mode = "ratio"
input_file = "data_EvolutionFlavors.txt"
range_particles = 6
c = 1.75


# order and name of particles
array_particles = ["gluon", "d", "u", "s", "c", "b", "t"]


# setting up the layout of the plot
fig, axs = plt.subplots(2, 3, sharex=True, figsize=(12*c,7*c))
subplt = [axs[0,2], axs[0,0], axs[1,0], axs[0,1], axs[1,1]]


# loading basic information on the following data
mu_vals = np.loadtxt(open(input_file), delimiter=",", unpack=False, skiprows=1, max_rows=1)
num_x_vals = int(np.loadtxt(open(input_file), delimiter=",", unpack=True, skiprows=3, max_rows=1))


# setting up arrays to save the data, val[particle, index of mu] gives the array of values
x = np.zeros((range_particles, len(mu_vals), num_x_vals))
apfel = np.zeros((range_particles, len(mu_vals), num_x_vals))
lhapdf = np.zeros((range_particles, len(mu_vals), num_x_vals))
ratio = np.zeros((range_particles, len(mu_vals), num_x_vals))


# setting up arrays used to calculate the singlet pdfs
a_singlet_lhapdf = np.zeros((len(mu_vals), num_x_vals))
a_singlet_apfel = np.zeros((len(mu_vals), num_x_vals))



###################
# main program
###################
# import the data
for i_mu in range(len(mu_vals)):

    for particle in range(range_particles):

        x[particle, i_mu], apfel[particle, i_mu], lhapdf[particle, i_mu], ratio[particle, i_mu] = np.loadtxt(open(input_file), delimiter=",", unpack=True, skiprows=((i_mu * len(array_particles) + particle) * num_x_vals + 5), max_rows=num_x_vals)


# calculating the values of the singlet pdfs
for particle in range(range_particles):
    # gluon is not included
    if (particle != 0):
        a_singlet_lhapdf += 2 * lhapdf[particle]
        a_singlet_apfel += 2 * apfel[particle]


# plot the values for the imported pdfs/ratios
for particle in range(range_particles):
    # bottom quark is not plotted
    if (particle != 5):
        for i_mu in range(len(mu_vals)):

            if (mode == "data"):
                subplt[particle].plot(x[particle, i_mu], lhapdf[particle, i_mu], label="lhapdf values at Q="+str(mu_vals[i_mu])+"GeV")
                subplt[particle].plot(x[particle, i_mu], apfel[particle, i_mu], label="apfel values at Q="+str(mu_vals[i_mu])+"GeV")
                subplt[particle].set_title(array_particles[particle]+" PDF")
                subplt[particle].tick_params('x', labelbottom=True)

            elif (mode == "ratio"):
                subplt[particle].plot(x[particle, i_mu], apfel[particle, i_mu]/lhapdf[particle, i_mu], label="apfel/lhapdf at Q="+str(mu_vals[i_mu])+"GeV")
                subplt[particle].set_title(array_particles[particle]+" PDF")
                subplt[particle].tick_params('x', labelbottom=True)


# plot singlet
for i_mu in range(len(mu_vals)):
    if (mode == "data"):
        axs[1,2].plot(x[0, i_mu], a_singlet_apfel[i_mu], label="apfel values at Q="+str(mu_vals[i_mu])+"GeV")
        axs[1,2].plot(x[0, i_mu], a_singlet_lhapdf[i_mu], label="lhapdf values at Q="+str(mu_vals[i_mu])+"GeV") 
        axs[1,2].set_title("singlet PDF")
        axs[1,2].tick_params('x', labelbottom=True)

    elif (mode == "ratio"):
        axs[1,2].plot(x[0, i_mu], a_singlet_apfel[i_mu]/a_singlet_lhapdf[i_mu], label="apfel/lhapdf at Q="+str(mu_vals[i_mu])+"GeV")
        axs[1,2].set_title("singlet PDF")
        axs[1,2].tick_params('x', labelbottom=True)
    


plt.xscale('log')
plt.xlim(left=10**(-4), right=1)
plt.legend(loc="upper right")
#plt.savefig('plotAllFlavors_LO.png', bbox_inches='tight')
plt.show()