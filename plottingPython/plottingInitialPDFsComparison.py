###################
# imports
###################
import numpy as np
import matplotlib.pyplot as plt



###################
# definitions
###################
c = 1.75
show_all_legends = False

order = "LO"

# save figure options
save_fig = False
pltdir = "/home/alexander/Documents/apfelxx_photon_mod/plots/mainPlots/" 
pltname = pltdir + "PlotInitialPdfs_LO_SAL3_Ks03_pto0_pl0.pdf"
dpi = 200 # default is 100

# x values to plot
a1 = np.array([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5])
a2 = np.array([1e-1, 1e-2, 1e-3, 1e-4])#, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9])
x = np.sort(np.outer(a2, a1).flatten())

# order and name of particles
array_particles = ["gluon", "d", "u", "s"]

# parameters if the SAL input PDFs (taken from photon_pdfs_v3, eq. (5.3)), at Q^2=4
SAL_Parameters       = [0.3, -0.57, 3, 0.065, -0.16, 1, 4.45, 1.9, 0.027]

# calculated parameters, using ALEPH1, ALEPH2, AMY, DELPHI, JADE, L3, OPAL1, PLUTO, TASSO, TOPAZ (134 points), at Q^2=4

if order == "LO":
    SAL6_Parameters      = [4.61032e-09, -0.183213, 3, 0.630064, 0.758714, 1, 0.382165, 0.482942, 0.61562]
    SAL4_Parameters      = [4.81049e-10, -0.314553, 3, 0.497431, 0.271455, 1, 0, 1, 0.31245] # = SAL6_had
    SAL5_KS03_Parameters = [0.3, -0.185909, 3, 0.565416, 0.734498, 1, 0.381047, 0.399418, 0.518411] # same as SAL6 but fixed KS, comparable to SAL
    SAL3_KS03_Parameters = [0.3, -0.298791, 3, 0.464334, 0.265692, 1, 0, 1, 0.214826] # = SAL5_had
    SAL3_KS05_Parameters = [0.5, -0.290716, 3, 0.444219, 0.261989, 1, 0, 1, 0.151961] # = SAL5_had
    array_chi2PerDP      = [0, 0.690308, 0.890824, 0.719782, 0.928777, 0.956312]

if order == "HO":
    SAL6_Parameters      = [2.66005e-09, -0.172879, 3, 0.640524, 0.614294, 1, 0.319267, 0.992109, 0.497544]
    SAL4_Parameters      = [0.3, -0.196667, 3, 0.584508, 0.606934, 1, 0.322353, 0.825024, 0.364456] # = SAL6_had
    SAL5_KS03_Parameters = [1.4222e-10, -0.112081, 3, 0.52892, 0.282323, 1, 0, 1, 0.407849] # same as SAL6 but fixed KS, comparable to SAL
    SAL3_KS03_Parameters = [0.3, -0.128374, 3, 0.493618, 0.277336, 1, 0, 1, 0.218647] # = SAL5_had
    SAL3_KS05_Parameters = [0.5, -0.136891, 3, 0.472198, 0.274209, 1, 0, 1, 0.113126] # = SAL5_had
    array_chi2PerDP      = [0, 0.768736, 0.797023, 0.860964, 0.898902, 0.924739]

# array of parameters
array_parameters = [SAL_Parameters, SAL6_Parameters, SAL4_Parameters, SAL5_KS03_Parameters, SAL3_KS03_Parameters, SAL3_KS05_Parameters]


# array of parameter labels
array_labels = ["original SAL input PDFs", "calculated SAL6 input PDFs", "calculated SAL4=SAL6_had",
                "calculated SAL5 input PDFs", "calculated SAL3=SAL5_had, $K_S=0.3$", "calculated SAL3=SAL5_had, $K_S=0.5$"]

# array of linestyles
array_linestyles = ["-", "-", "-.", "-", "-.", "-."]

# array of colors
array_colors = ["tab:blue", "tab:green", "tab:green", "tab:orange", "tab:orange", "tab:olive"]



###################
# loading GRV data
###################
if order == "LO":
    GRV_file = "/home/alexander/Documents/apfelxx_photon_mod/results/Bestandsaufnahme_2022_05_01/dataGRVInitialPDFsLO.txt" 
    GRV_label = 'GRVLO "Initial PDFs"'
if order == "HO":
    GRV_file = "/home/alexander/Documents/apfelxx_photon_mod/results/Bestandsaufnahme_2022_05_01/dataGRVInitialPDFsHO.txt" 
    GRV_label = 'GRVHO "Initial PDFs"'
range_particles = 6

# loading basic information on the following data
num_x_vals = int(np.loadtxt(open(GRV_file), delimiter=",", unpack=True, skiprows=3, max_rows=1))


# setting up arrays to save the data, val[particle, index of mu] gives the array of values
x_GRV = np.zeros((range_particles, num_x_vals))
apfel = np.zeros((range_particles, num_x_vals))
GRV_values = np.zeros((range_particles, num_x_vals))
ratio = np.zeros((range_particles, num_x_vals))

# import the data
for i_mu in range(1):

    for particle in range(range_particles):

        x_GRV[particle], apfel[particle], GRV_values[particle], ratio[particle] = np.loadtxt(open(GRV_file), delimiter=",", unpack=True, skiprows=((i_mu * len(array_particles) + particle) * num_x_vals + 5), max_rows=num_x_vals)


###################
# functions
###################
def InitialPDFgluonSAL(x, A, B, C):
    return A * x**B * (1-x)**C

def InitialPDFdownSAL(x, A_had, B_had, C_had, A_pl, B_pl):
    return A_had * x**B_had * (1-x)**C_had + 1/9 * A_pl * x * (x**2 + (1-x)**2)/(1-B_pl*np.log(1-x))

def InitialPDFupSAL(x, A_had, B_had, C_had, A_pl, B_pl):
    return A_had * x**B_had * (1-x)**C_had + 4/9 * A_pl * x * (x**2 + (1-x)**2)/(1-B_pl*np.log(1-x))

def InitialPDFstrangeSAL(x, K, A_had, B_had, C_had, A_pl, B_pl):
    return K * A_had * x**B_had * (1-x)**C_had + 1/9 * A_pl * x * (x**2 + (1-x)**2)/(1-B_pl*np.log(1-x))



###################
# plotting
###################
# setting up the layout of the plot
fig, axs = plt.subplots(2, 2, sharex=True, figsize=(12*c,7*c))
subplt = [axs[1,1], axs[1,0], axs[0,1], axs[0,0]]


# plot GRV values
if True:
    # plot gluon
    subplt[0].plot(x_GRV[0], GRV_values[0], label=GRV_label, linestyle="--", color="tab:blue")
    subplt[0].set_title("gluon InitialPDF")

    subplt[0].tick_params('x', labelbottom=True)
    subplt[0].legend(loc="upper right")

    # plot down quark
    subplt[2].plot(x_GRV[1], GRV_values[1], label=GRV_label, linestyle="--", color="tab:blue")
    subplt[2].set_title("down quark InitialPDF")

    subplt[2].tick_params('x', labelbottom=True)
    if show_all_legends:
        subplt[2].legend(loc="upper left")

    # plot up quark
    subplt[3].plot(x_GRV[2], GRV_values[2], label=GRV_label, linestyle="--", color="tab:blue")
    subplt[3].set_title("up quark InitialPDF")

    subplt[3].tick_params('x', labelbottom=True)
    if show_all_legends:
        subplt[3].legend(loc="upper left")

    # plot strange quark
    subplt[1].plot(x_GRV[3], GRV_values[3], label=GRV_label, linestyle="--", color="tab:blue")
    subplt[1].set_title("strange quark InitialPDF")

    subplt[1].tick_params('x', labelbottom=True)
    if show_all_legends:
        subplt[1].legend(loc="upper left")


# plotting the SAL values
for i in range(len(array_labels)):
    Parameters = array_parameters[i]

    # plot gluon
    subplt[0].plot(x, InitialPDFgluonSAL(x, Parameters[8], Parameters[1], Parameters[2]), label=array_labels[i], linestyle=array_linestyles[i], color=array_colors[i])
    subplt[0].set_title("gluon InitialPDF")

    subplt[0].tick_params('x', labelbottom=True)
    subplt[0].legend(loc="upper right")

    # plot down quark
    subplt[2].plot(x, InitialPDFdownSAL(x, Parameters[3], Parameters[4], Parameters[5], Parameters[6], Parameters[7]), label=array_labels[i], linestyle=array_linestyles[i], color=array_colors[i])
    subplt[2].set_title("down quark InitialPDF")

    subplt[2].tick_params('x', labelbottom=True)
    if show_all_legends:
        subplt[2].legend(loc="upper left")

    # plot up quark
    subplt[3].plot(x, InitialPDFupSAL(x, Parameters[3], Parameters[4], Parameters[5], Parameters[6], Parameters[7]), label=array_labels[i], linestyle=array_linestyles[i], color=array_colors[i])
    subplt[3].set_title("up quark InitialPDF")

    subplt[3].tick_params('x', labelbottom=True)
    if show_all_legends:
        subplt[3].legend(loc="upper left")

    # plot strange quark
    subplt[1].plot(x, InitialPDFstrangeSAL(x, Parameters[0], Parameters[3], Parameters[4], Parameters[5], Parameters[6], Parameters[7]), label=array_labels[i], linestyle=array_linestyles[i], color=array_colors[i])
    subplt[1].set_title("strange quark InitialPDF")

    subplt[1].tick_params('x', labelbottom=True)
    if show_all_legends:
        subplt[1].legend(loc="upper left")


#title = "InitialPDFs: " + str(usedInitialPDFs)
#title += ",\nchi2: " + str(chi2) + ", chi2/NumberOfDataPoints: " + str(chi2PerDP) + ",\nexperimentalData: "
#for text in usedExperimentalData:
#    title += text + "," 
#title = title[:len(title)-1] + ",\nParameters: "
#for i in range(len(Parameters)):
#    title += ParametersNames[i] + "=" + str(Parameters[i]) + ", " 
#title = title[:len(title)-2]

#plt.suptitle(title)
plt.xscale('log')
plt.xlim(left=10**(-4), right=1)
if save_fig:
    plt.savefig(pltname, bbox_inches='tight', dpi=dpi)
plt.show()