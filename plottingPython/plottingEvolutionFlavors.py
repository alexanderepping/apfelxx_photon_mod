import numpy as np
import matplotlib.pyplot as plt
import math

mode = "ratio"
arr_mu = [1.295, 10]


num_x_vals = 36



array_particles = ["gluon", "d", "u", "s", "c", "b", "t"]

fig, axs = plt.subplots(2, 3, sharex=True, figsize=(12,7))
subplt = [axs[0,2], axs[0,0], axs[1,0], axs[0,1], axs[1,1]]

a_singlet_lhapdf = np.array([])
a_singlet_apfel = np.array([])

if (mode == "data"):
    for i_mu in range(2):
    #for particle in range(5):

        for particle in range(6):
        #for i_mu in range(2):
            
            # import the data
            x, apfel, lhapdf, ratio = np.loadtxt(open("data_EvolutionFlavors.txt"), delimiter="	", unpack=True, skiprows=((i_mu*7 + particle)*num_x_vals), max_rows=36)
            lbl = "ratio of apfel/lhapdf for " + array_particles[particle] 
            ylabel = "x * " + array_particles[particle] + " / alpha  PDF"


            if (particle==0):
                a_singlet_lhapdf = x*0
                a_singlet_apfel = x*0
            else:
                a_singlet_lhapdf += 2*lhapdf
                a_singlet_apfel += 2*apfel
            
            if (particle!=5):
                subplt[particle].plot(x, apfel, label="apfel values at Q="+str(arr_mu[i_mu])+"GeV")
                subplt[particle].plot(x, lhapdf, label="lhapdf values at Q="+str(arr_mu[i_mu])+"GeV")
                subplt[particle].set_title(array_particles[particle]+" PDF")
                subplt[particle].tick_params('x', labelbottom=True)
                

        axs[1,2].plot(x, a_singlet_apfel, label="apfel values at Q="+str(arr_mu[i_mu])+"GeV")
        axs[1,2].plot(x, a_singlet_lhapdf, label="lhapdf values at Q="+str(arr_mu[i_mu])+"GeV") 
        axs[1,2].set_title("singlet PDF")
        axs[1,2].tick_params('x', labelbottom=True)

    plt.xscale('log')
    plt.xlim(left=10**(-4), right=1)
    plt.legend(loc="upper right")
    #plt.savefig('plotAllFlavors.png', bbox_inches='tight')
    plt.show()

if (mode=="ratio"):
    for i_mu in range(2):
    #for particle in range(5):

        for particle in range(6):
        #for i_mu in range(2):
            
            # import the data
            x, apfel, lhapdf, ratio = np.loadtxt(open("data_EvolutionFlavors.txt"), delimiter="	", unpack=True, skiprows=((i_mu*7 + particle)*num_x_vals), max_rows=36)
            lbl = "ratio of apfel/lhapdf for " + array_particles[particle] 
            ylabel = "x * " + array_particles[particle] + " / alpha  PDF"


            if (particle==0):
                a_singlet_lhapdf = x*0
                a_singlet_apfel = x*0
            else:
                a_singlet_lhapdf += 2*lhapdf
                a_singlet_apfel += 2*apfel
            
            if (particle!=5):
                subplt[particle].plot(x, apfel/lhapdf, label="apfel/lhapdf at Q="+str(arr_mu[i_mu])+"GeV")
                subplt[particle].set_title(array_particles[particle]+" PDF")
                subplt[particle].tick_params('x', labelbottom=True)
                

        axs[1,2].plot(x, a_singlet_apfel/a_singlet_lhapdf, label="apfel/lhapdf at Q="+str(arr_mu[i_mu])+"GeV")
        axs[1,2].set_title("singlet PDF")
        axs[1,2].tick_params('x', labelbottom=True)

    plt.xscale('log')
    plt.xlim(left=10**(-4), right=1)
    plt.legend(loc="upper right")
    #plt.savefig('plotAllFlavors.png', bbox_inches='tight')
    plt.show()


