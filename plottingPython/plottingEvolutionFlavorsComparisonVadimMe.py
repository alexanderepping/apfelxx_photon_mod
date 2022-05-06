###################
# imports
###################
import numpy as np
import matplotlib.pyplot as plt
import os



###################
# definitions
###################
# variables changeable by user
mode = "ratio"
ylim_bool = True
ylim = 0.02
input_file_Vadim = "/home/alexander/Uni/apfelxx_photon_mod/results/FinalEvolution01_2022_04_14/Vadim_data_LO_Q4_2022_04_13.dat"
input_file = "data_EvolutionFlavors.txt"
range_particles = 6
range_particles_vadim = 8
c = 1.75

# make path to the correct file
if os.path.dirname(__file__) != "":
    input_file = os.path.dirname(__file__) + "/" + input_file

# order and name of particles
array_particles = ["gluon", "down quark", "up quark", "strange quark", "charm quark", "b", "t"]
array_particles_short = ["g", "d", "u", "s", "c", "b", "t"]


# setting up the layout of the plot
fig, axs = plt.subplots(2, 3, sharex=True, figsize=(12*c,7*c))
subplt = [axs[0,2], axs[0,0], axs[1,0], axs[0,1], axs[1,1]]


# loading basic information on the following data
mu_vals = np.array([4.0])
num_x_vals = int(np.loadtxt(open(input_file), delimiter=",", unpack=True, skiprows=3, max_rows=1))


# setting up arrays to save the data, val[particle, index of mu] gives the array of values
x               = np.zeros((range_particles, len(mu_vals), num_x_vals))
x_Vadim         = np.zeros((range_particles, len(mu_vals), num_x_vals))
apfel           = np.zeros((range_particles, len(mu_vals), num_x_vals))
qcdnum_vadim    = np.zeros((range_particles_vadim, len(mu_vals), num_x_vals))
lhapdf          = np.zeros((range_particles, len(mu_vals), num_x_vals))
ratio_alex      = np.zeros((range_particles, len(mu_vals), num_x_vals))
ratio_VadimAlex = np.zeros((range_particles, len(mu_vals), num_x_vals))

# setting up arrays used to calculate the singlet pdfs
a_singlet_lhapdf= np.zeros((len(mu_vals), num_x_vals))
a_singlet_apfel = np.zeros((len(mu_vals), num_x_vals))



###################
# main program
###################
# import the data Alex
for i_mu in range(len(mu_vals)):

    for particle in range(range_particles):

        x[particle, i_mu], apfel[particle, i_mu], lhapdf[particle, i_mu], ratio_alex[particle, i_mu] = np.loadtxt(open(input_file), delimiter=",", unpack=True, skiprows=((i_mu * len(array_particles) + particle) * num_x_vals + 5), max_rows=num_x_vals)

# import the data Vadim
x_Vadim[0,0], qcdnum_vadim[5, 0], qcdnum_vadim[6, 0], qcdnum_vadim[2, 0], qcdnum_vadim[1, 0], qcdnum_vadim[3, 0], qcdnum_vadim[4, 0], qcdnum_vadim[0, 0], qcdnum_vadim[7, 0] = np.loadtxt(open(input_file_Vadim), delimiter=",", unpack=True, skiprows=1, max_rows=num_x_vals)

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

            subplt[particle].plot(x[particle, i_mu], apfel[particle, i_mu]/qcdnum_vadim[particle, i_mu], label="Apfel++/QCDNUM at Q="+str(mu_vals[i_mu])+"GeV")
            subplt[particle].plot(x[particle, i_mu], apfel[particle, i_mu]/apfel[particle, i_mu])
            subplt[particle].set_title(array_particles[particle]+" PDF")
            subplt[particle].tick_params('x', labelbottom=True)
            subplt[particle].legend(loc="upper left")
            #if (particle != 0):
            #    subplt[particle].legend(loc="upper left")
            #else:
            #    subplt[particle].legend(loc="upper right")
            subplt[particle].set_xlabel("$x$")
            if ylim_bool:
                if (particle != 4):
                    #subplt[particle].set_ylim([1-ylim, 1+ylim])
                    subplt[particle].set_ylim([1-0.01, 1+ylim])
                else:
                    subplt[particle].set_ylim([1-0.005, 1+0.025])


# plot singlet
for i_mu in range(len(mu_vals)):
    axs[1,2].plot(x[0, i_mu], a_singlet_apfel[i_mu]/qcdnum_vadim[7,i_mu], label="Apfel++/QCDNUM at Q="+str(mu_vals[i_mu])+"GeV")
    axs[1,2].plot(x[0, i_mu], qcdnum_vadim[7,i_mu]/qcdnum_vadim[7,i_mu])
    axs[1,2].set_title("singlet PDF")
    axs[1,2].tick_params('x', labelbottom=True)
    axs[1,2].legend(loc="upper left")
    axs[1,2].set_xlabel("$x$")
    if ylim_bool:
        #axs[1,2].set_ylim([1-ylim, 1+ylim])
        axs[1,2].set_ylim([1-0.01, 1+ylim])
    


plt.xscale('log')
plt.xlim(left=10**(-4), right=1)
#plt.savefig('/home/alexander/Uni/apfelxx_photon_mod/results/FinalEvolution01_2022_04_14/plot_LO_ratioVadim_2022_04_14.pdf', bbox_inches='tight')
plt.show()