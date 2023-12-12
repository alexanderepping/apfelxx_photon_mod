###################
# Imports
###################
import numpy as np


###################
# InitialPDFs
###################
def InitialPDFgluonSAL(x, params):
    A = params[8]
    B = params[1]
    C = params[2]
    return A * x**B * (1-x)**C

def InitialPDFdownSAL(x, params):
    A_had = params[3]
    B_had = params[4]
    C_had = params[5]
    A_pl  = params[6]
    B_pl  = params[7]
    return A_had * x**B_had * (1-x)**C_had + 1/9 * A_pl * x * (x**2 + (1-x)**2)/(1-B_pl*np.log(1-x))
    #return A_had * x**B_had * (1-x)**C_had

def InitialPDFupSAL(x, params):
    A_had = params[3]
    B_had = params[4]
    C_had = params[5]
    A_pl  = params[6]
    B_pl  = params[7]
    return A_had * x**B_had * (1-x)**C_had + 4/9 * A_pl * x * (x**2 + (1-x)**2)/(1-B_pl*np.log(1-x))
    #return A_had * x**B_had * (1-x)**C_had

def InitialPDFstrangeSAL(x, params):
    K     = params[0]
    A_had = params[3]
    B_had = params[4]
    C_had = params[5]
    A_pl  = params[6]
    B_pl  = params[7]
    return K * A_had * x**B_had * (1-x)**C_had + 1/9 * A_pl * x * (x**2 + (1-x)**2)/(1-B_pl*np.log(1-x))
    #return K * A_had * x**B_had * (1-x)**C_had


###################
# ErrorPDFs
###################
def DeltaErrorPDF(InitialPDF, x, errorParamsPlus, errorParamsMinus):
    result = 0.0
    for k in range(len(errorParamsPlus)):
        result += (InitialPDF(x, errorParamsPlus[k]) - InitialPDF(x, errorParamsMinus[k]))**2
    result = 0.5 * np.sqrt(result)
    return result

def LowerErrorPDF(InitialPDF, x, params, errorParamsPlus, errorParamsMinus):
    return InitialPDF(x, params) - DeltaErrorPDF(InitialPDF, x, errorParamsPlus, errorParamsMinus)

def UpperErrorPDF(InitialPDF, x, params, errorParamsPlus, errorParamsMinus):
    return InitialPDF(x, params) + DeltaErrorPDF(InitialPDF, x, errorParamsPlus, errorParamsMinus)

###################
# Other
###################
def makeTitle(usedInitialPDFs, chi2, chi2PerDP, usedExperimentalData, Parameters, ParametersNames, ErrorPDFs=False, deltaChi2=1):
    title = "InitialPDFs: " + str(usedInitialPDFs)
    title += ",\nchi2: " + str(chi2) + ", chi2/NumberOfDataPoints: " + str(chi2PerDP) 
    if ErrorPDFs:
        title += ", delta chi2: " + str(deltaChi2)
    title += ",\nexperimentalData: "
    for text in usedExperimentalData:
        title += text + "," 
    title = title[:len(title)-1] + ",\nParameters: "
    for i in range(len(Parameters)):
        title += ParametersNames[i] + "=" + str(Parameters[i]) + ", " 
    title = title[:len(title)-2]
    return title