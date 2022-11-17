import os



#######################
# paths to directories
#######################
dirThisFile = os.path.dirname(__file__) + "/"
dirApfel = "/home/alexander/Uni/apfelxx_photon_mod/"
dirReference = dirApfel+"results/ReferenceResults/"

dirBestandsaufnahme = dirApfel+"results/Bestandsaufnahme_2022_11_13/"
dirBestandsaufnahme2 = dirApfel+"results/Calculations_2022_10_13/"
dirResults = dirApfel+"results/"
dirCurrent = dirBestandsaufnahme



##############################
# single data sets
##############################
data_SAL3LO                   = {"FilePath":   dirCurrent+"dataEvolvedPDFsSAL3LO.txt",
                                 "DataColumn": 1}                   
data_SAL4VadimLO              = {"FilePath":   dirCurrent+"dataEvolvedPDFsSAL4VADIMLO.txt",
                                 "DataColumn": 1}                   
data_SAL5LO                   = {"FilePath":   dirCurrent+"dataEvolvedPDFsSAL5LO.txt",
                                 "DataColumn": 1}                   
data_SAL3HO                   = {"FilePath":   dirCurrent+"dataEvolvedPDFsSAL3HO.txt",
                                 "DataColumn": 1}                   
data_SAL4VadimHO              = {"FilePath":   dirCurrent+"dataEvolvedPDFsSAL4VADIMHO.txt",
                                 "DataColumn": 1}                   
data_SAL5HO                   = {"FilePath":   dirCurrent+"dataEvolvedPDFsSAL5HO.txt",
                                 "DataColumn": 1}                   

error_SAL3LO                   = {"FilePath":   dirCurrent+"dataEvolvedPDFsSAL3LO.txt",
                                 "DataColumn": 4}                   
error_SAL4VadimLO              = {"FilePath":   dirCurrent+"dataEvolvedPDFsSAL4VADIMLO.txt",
                                 "DataColumn": 4}                   
error_SAL5LO                   = {"FilePath":   dirCurrent+"dataEvolvedPDFsSAL5LO.txt",
                                 "DataColumn": 4}                   
error_SAL3HO                   = {"FilePath":   dirCurrent+"dataEvolvedPDFsSAL3HO.txt",
                                 "DataColumn": 4}                   
error_SAL4VadimHO              = {"FilePath":   dirCurrent+"dataEvolvedPDFsSAL4VADIMHO.txt",
                                 "DataColumn": 4}                   
error_SAL5HO                   = {"FilePath":   dirCurrent+"dataEvolvedPDFsSAL5HO.txt",
                                 "DataColumn": 4}                   

data_EvolvedPDFs1             = {"FilePath":   dirThisFile+"dataEvolvedPDFs.txt",
                                 "DataColumn": 1}                   
data_EvolvedPDFs2             = {"FilePath":   dirThisFile+"dataEvolvedPDFs.txt",
                                 "DataColumn": 2}                   
data_EvolvedPDFs3             = {"FilePath":   dirThisFile+"dataEvolvedPDFs.txt",
                                 "DataColumn": 3}                   
data_EvolvedPDFs4             = {"FilePath":   dirThisFile+"dataEvolvedPDFs.txt",
                                 "DataColumn": 4}                   

data_EvolvedPDFs5             = {"FilePath":   dirBestandsaufnahme2+"dataEvolvedPDFsSAL5HOErrors_withMSR.txt",
                                 "DataColumn": 1}                   
data_EvolvedPDFs5E            = {"FilePath":   dirBestandsaufnahme2+"dataEvolvedPDFsSAL5HOErrors_withMSR.txt",
                                 "DataColumn": 4}                   
data_EvolvedPDFs6             = {"FilePath":   dirBestandsaufnahme2+"dataEvolvedPDFsSAL5HOErrors_woutMSR.txt",
                                 "DataColumn": 1}                   
data_EvolvedPDFs6E            = {"FilePath":   dirBestandsaufnahme2+"dataEvolvedPDFsSAL5HOErrors_woutMSR.txt",
                                 "DataColumn": 4}                   


#######################
# data set directories
#######################




##########################################################################################
# REFERENCE DATA SETS
##########################################################################################

##############################
# single reference data sets
##############################
# GRVLO, taken at 1.3 GeV and evolved to sqrt2 GeV
data_GRVLOSqrt2_Evolved       = {"FilePath":   dirReference+"dataPDFsGRVLOSqrt2.txt",
                                 "DataColumn": 1}                      
# GRVLO, taken at sqrt2 GeV
data_GRVLOSqrt2_Exact         = {"FilePath":   dirReference+"dataPDFsGRVLOSqrt2.txt",
                                 "DataColumn": 2}                   
# GRVHO, taken at 1.3 GeV and evolved to sqrt2 GeV
data_GRVHOSqrt2_Evolved       = {"FilePath":   dirReference+"dataPDFsGRVHOSqrt2.txt",
                                 "DataColumn": 1}                      
# GRVHO, taken at sqrt2 GeV
data_GRVHOSqrt2_Exact         = {"FilePath":   dirReference+"dataPDFsGRVHOSqrt2.txt",
                                 "DataColumn": 2}                   
# SAL data, taken at sqrt2 GeV
data_SAL                      = {"FilePath":   dirReference+"dataPDFsSAL.txt",
                                 "DataColumn": 1}                   



###################################
# reference data set directories
###################################
DataComparisonEvolutionGRVLO = { #calculated using LHAPDF-Set GRVLO_GRVParameters and precompiler-setting asGRV
            "GRVLO13_Evolved":  {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVLO13_GRVParameters.txt",
                                 "Label":      "Apfel++ at 1.3 GeV",   
                                 "DataColumn": 1,                 
                                 "XColumn":    0},               
            "GRVLO13_Exact":    {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVLO13_GRVParameters.txt",
                                 "Label":      "GRV at 1.3 GeV",   
                                 "DataColumn": 2,                 
                                 "XColumn":    0},               
            "GRVLO40_Evolved":  {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVLO40_GRVParameters.txt",
                                 "Label":      "Apfel++ at 4.0 GeV",   
                                 "DataColumn": 1,                 
                                 "XColumn":    0},               
            "GRVLO40_Exact":    {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVLO40_GRVParameters.txt",
                                 "Label":      "GRV at 4.0 GeV",   
                                 "DataColumn": 2,                 
                                 "XColumn":    0},
            "GRVLO100_Evolved":  {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVLO100_GRVParameters.txt",
                                 "Label":      "Apfel++ at 10 GeV",   
                                 "DataColumn": 1,                 
                                 "XColumn":    0},               
            "GRVLO100_Exact":    {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVLO100_GRVParameters.txt",
                                 "Label":      "GRV at 10 GeV",   
                                 "DataColumn": 2,                 
                                 "XColumn":    0},               
            "GRVLO200_Evolved":  {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVLO200_GRVParameters.txt",
                                 "Label":      "Apfel++ at 20 GeV",   
                                 "DataColumn": 1,                 
                                 "XColumn":    0},               
            "GRVLO200_Exact":    {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVLO200_GRVParameters.txt",
                                 "Label":      "GRV at 20 GeV",   
                                 "DataColumn": 2,                 
                                 "XColumn":    0}               
            }

DataComparisonEvolutionGRVLORatio = { #calculated using LHAPDF-Set GRVLO_GRVParameters and precompiler-setting asGRV
            "GRVLO13":          {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVLO13_GRVParameters.txt",
                                 "Label":      "Apfel++/GRV at 1.3 GeV",   
                                 "DataColumn": 3,                 
                                 "XColumn":    0},               
            "GRVLO40":          {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVLO40_GRVParameters.txt",
                                 "Label":      "Apfel++/GRV at 4.0 GeV",   
                                 "DataColumn": 3,                 
                                 "XColumn":    0},
            "GRVLO100":         {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVLO100_GRVParameters.txt",
                                 "Label":      "Apfel++/GRV at 10 GeV",   
                                 "DataColumn": 3,                 
                                 "XColumn":    0},               
            "GRVLO200":         {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVLO200_GRVParameters.txt",
                                 "Label":      "Apfel++/GRV at 20 GeV",   
                                 "DataColumn": 3,                 
                                 "XColumn":    0},               
            }

DataComparisonEvolutionGRVHO = { #calculated using LHAPDF-Set GRVHO_GRVParameters and precompiler-setting asGRV
            "GRVHO13_Evolved":  {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVHO13_GRVParameters.txt",
                                 "Label":      "Apfel++ at 1.3 GeV",   
                                 "DataColumn": 1,                 
                                 "XColumn":    0},               
            "GRVHO13_Exact":    {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVHO13_GRVParameters.txt",
                                 "Label":      "GRV at 1.3 GeV",   
                                 "DataColumn": 2,                 
                                 "XColumn":    0},               
            "GRVHO40_Evolved":  {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVHO40_GRVParameters.txt",
                                 "Label":      "Apfel++ at 4.0 GeV",   
                                 "DataColumn": 1,                 
                                 "XColumn":    0},               
            "GRVHO40_Exact":    {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVHO40_GRVParameters.txt",
                                 "Label":      "GRV at 4.0 GeV",   
                                 "DataColumn": 2,                 
                                 "XColumn":    0},
            "GRVHO100_Evolved":  {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVHO100_GRVParameters.txt",
                                 "Label":      "Apfel++ at 10 GeV",   
                                 "DataColumn": 1,                 
                                 "XColumn":    0},               
            "GRVHO100_Exact":    {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVHO100_GRVParameters.txt",
                                 "Label":      "GRV at 10 GeV",   
                                 "DataColumn": 2,                 
                                 "XColumn":    0},               
            "GRVHO200_Evolved":  {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVHO200_GRVParameters.txt",
                                 "Label":      "Apfel++ at 20 GeV",   
                                 "DataColumn": 1,                 
                                 "XColumn":    0},               
            "GRVHO200_Exact":    {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVHO200_GRVParameters.txt",
                                 "Label":      "GRV at 20 GeV",   
                                 "DataColumn": 2,                 
                                 "XColumn":    0}               
            }

DataComparisonEvolutionGRVHORatio = { #calculated using LHAPDF-Set GRVHO_GRVParameters and precompiler-setting asGRV
            "GRVHO13":          {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVHO13_GRVParameters.txt",
                                 "Label":      "Apfel++/GRV at 1.3 GeV",   
                                 "DataColumn": 3,                 
                                 "XColumn":    0},               
            "GRVHO40":          {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVHO40_GRVParameters.txt",
                                 "Label":      "Apfel++/GRV at 4.0 GeV",   
                                 "DataColumn": 3,                 
                                 "XColumn":    0},
            "GRVHO100":         {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVHO100_GRVParameters.txt",
                                 "Label":      "Apfel++/GRV at 10 GeV",   
                                 "DataColumn": 3,                 
                                 "XColumn":    0},               
            "GRVHO200":         {"FilePath":   dirReference+"ComparisonEvolutionGRV/dataPDFsGRVHO200_GRVParameters.txt",
                                 "Label":      "Apfel++/GRV at 20 GeV",   
                                 "DataColumn": 3,                 
                                 "XColumn":    0},               
            }