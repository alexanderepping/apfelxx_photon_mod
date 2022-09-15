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
OutputFileName = "allData_2022_06_23.txt"
range_particles = 6
c = 1.75

# make path to the correct file
filepath = ""
if os.path.dirname(__file__) != "":
    filepath = os.path.dirname(__file__) + "/"

# order and name of particles
array_particles = ["gluon", "down quark", "up quark", "strange quark", "charm quark", "b", "t"]
array_particles_short = ["g", "d", "u", "s", "c", "b", "t"]


# setting up the layout of the plot
fig, axs = plt.subplots(2, 3, sharex=True, figsize=(12*c,7*c))
subplt = [axs[0,2], axs[0,0], axs[1,0], axs[0,1], axs[1,1]]



###################
# loading the data
###################

# loading basic information on the following data
num_x_vals = int(np.loadtxt(open(filepath+"data_GRVHOSqrt2.txt"), delimiter=",", unpack=True, skiprows=3, max_rows=1))

x           = np.zeros((range_particles, num_x_vals))
GRVEvolved  = np.zeros((range_particles, num_x_vals))
GRVExact    = np.zeros((range_particles, num_x_vals))
SAL         = np.zeros((range_particles, num_x_vals))
SAL3Evolved = np.zeros((range_particles, num_x_vals))
SAL5Evolved = np.zeros((range_particles, num_x_vals))
ratio       = np.zeros((range_particles, num_x_vals))

for particle in range(range_particles):
    # GRV HO data at Q=sqrt(2)
    x[particle], GRVEvolved[particle],  GRVExact[particle], ratio[particle] = np.loadtxt(open(filepath+"data_GRVHOSqrt2.txt"),    delimiter=",", unpack=True, skiprows=(particle * num_x_vals + 5), max_rows=num_x_vals)
    # SAL3HOEvolved at Q=sqrt(2)
    x[particle], SAL3Evolved[particle], SAL[particle],      ratio[particle] = np.loadtxt(open(filepath+"data_SAL3HOEvolved.txt"), delimiter=",", unpack=True, skiprows=(particle * num_x_vals + 5), max_rows=num_x_vals)
    # SAL5HOEvolved at Q=sqrt(2)
    x[particle], SAL5Evolved[particle], SAL[particle],      ratio[particle] = np.loadtxt(open(filepath+"data_SAL5HOEvolved.txt"), delimiter=",", unpack=True, skiprows=(particle * num_x_vals + 5), max_rows=num_x_vals)



# setting up arrays used to calculate the singlet pdfs
Singlet_GRVEvolved  = np.zeros(num_x_vals)
Singlet_GRVExact    = np.zeros(num_x_vals)
Singlet_SAL         = np.zeros(num_x_vals)
Singlet_SAL3Evolved = np.zeros(num_x_vals)
Singlet_SAL5Evolved = np.zeros(num_x_vals)

# calculating the values of the singlet pdfs
for particle in range(range_particles):
    # gluon is not included
    if (particle != 0):
        Singlet_GRVEvolved  += 2 * GRVEvolved[particle]
        Singlet_GRVExact    += 2 * GRVExact[particle]   
        Singlet_SAL         += 2 * SAL[particle]        
        Singlet_SAL3Evolved += 2 * SAL3Evolved[particle]
        Singlet_SAL5Evolved += 2 * SAL5Evolved[particle]



###################
# write data to file
###################
OutputFile = open(filepath+OutputFileName, "w")

# mu values:
1.414214
# num_x_vals
36
# x, apfel++, lhapdf, ratio
OutputFile.write("# mu values:\n")
OutputFile.write("1.414214\n")
OutputFile.write("# num_x_vals\n")
OutputFile.write(str(num_x_vals)+"\n")
OutputFile.write("# x, GRVEvolved, GRVExact, SAL, SAL3Evolved, SAL5Evolved\n")

for particle in range(range_particles):
    for i in range(num_x_vals):
        OutputFile.write(str(x[particle][i])+", ")
        OutputFile.write(str(GRVEvolved[particle][i])+", ")
        OutputFile.write(str(GRVExact[particle][i])+", ")
        OutputFile.write(str(SAL[particle][i])+", ")
        OutputFile.write(str(SAL3Evolved[particle][i])+", ")
        OutputFile.write(str(SAL5Evolved[particle][i])+"\n")

for i in range(num_x_vals):
    OutputFile.write(str(x[particle][i])+", ")
    OutputFile.write(str(Singlet_GRVEvolved[i])+", ")
    OutputFile.write(str(Singlet_GRVExact[i])+", ")
    OutputFile.write(str(Singlet_SAL[i])+", ")
    OutputFile.write(str(Singlet_SAL3Evolved[i])+", ")
    OutputFile.write(str(Singlet_SAL5Evolved[i])+"\n")

OutputFile.close()



###################
# main program
###################


# plot the values for the imported pdfs/ratios
for particle in range(range_particles):
    # bottom quark is not plotted
    if (particle != 5):
        subplt[particle].plot(x[particle], GRVEvolved[particle],  label="GRV, taken at 1.3 GeV,\nevolved to sqrt(2) GeV")
        subplt[particle].plot(x[particle], GRVExact[particle],    label="GRV, taken at sqrt(2) GeV")
        subplt[particle].plot(x[particle], SAL[particle],         label="SAL at sqrt(2) GeV")
        subplt[particle].plot(x[particle], SAL3Evolved[particle], label="calculated InitialPDF SAL3,\nevolved from 1.3 GeV to sqrt(2) GeV")
        subplt[particle].plot(x[particle], SAL5Evolved[particle], label="calculated InitialPDF SAL5,\nevolved from 1.3 GeV to sqrt(2) GeV")

        subplt[particle].set_title(array_particles[particle]+" PDF")
        subplt[particle].tick_params('x', labelbottom=True)
        subplt[particle].legend(loc="upper right")
        subplt[particle].set_xlabel("$x$")
        subplt[particle].set_ylabel("$x"+array_particles_short[particle ]+"/\\alpha_{QED}$")



# plot singlet
axs[1,2].plot(x[particle], Singlet_GRVEvolved,            label="GRV, taken at 1.3 GeV,\n evolved to sqrt(2) GeV")
axs[1,2].plot(x[particle], Singlet_GRVExact,              label="GRV, taken at sqrt(2) GeV")
axs[1,2].plot(x[particle], Singlet_SAL,                   label="SAL at sqrt(2) GeV")
axs[1,2].plot(x[particle], Singlet_SAL3Evolved,           label="calculated InitialPDF SAL3,\nevolved from 1.3 GeV to sqrt(2) GeV")
axs[1,2].plot(x[particle], Singlet_SAL5Evolved,           label="calculated InitialPDF SAL5,\nevolved from 1.3 GeV to sqrt(2) GeV")
axs[1,2].set_title("singlet PDF")
axs[1,2].tick_params('x', labelbottom=True)
axs[1,2].legend(loc="upper right")
axs[1,2].set_xlabel("$x$")
axs[1,2].set_ylabel("$x\Sigma(x)/\\alpha_{QED}$")
    


plt.xscale('log')
plt.xlim(left=10**(-4), right=1)
plt.savefig(filepath+"plotAllDataHOatSqrt2.pdf", bbox_inches='tight')
plt.savefig(filepath+"plotAllDataHOatSqrt2.png", bbox_inches='tight')
plt.show()