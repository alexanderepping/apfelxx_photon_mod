import os



#######################
# paths to directories
#######################
dirThisFile = os.path.dirname(__file__) + "/"
dirApfel = "/home/alexander/Uni/apfelxx_photon_mod/"
dirBestandsaufnahme = dirApfel+"results/Bestandsaufnahme_2023_10_16/"


def DataComparisonEvolutionGRV(order):
    output = { # For Plot 1 & 2
            "GRV13_Evolved":      {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFsGRV"+order+"asGRV28_13.txt",
                                   "Label":      "APFEL++, (1.3)² GeV²",   
                                   "DataColumn": 1,                 
                                   "XColumn":    0},               
            "GRV13_Exact":        {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFsGRV"+order+"asGRV28_13.txt",
                                   "Label":      "GRV, (1.3)² GeV²",   
                                   "DataColumn": 2,                 
                                   "XColumn":    0},               
            "GRVsqrt10_Evolved":  {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFsGRV"+order+"asGRV28_sqrt10.txt",
                                   "Label":      "APFEL++, 10 GeV²",   
                                   "DataColumn": 1,                 
                                   "XColumn":    0},               
            "GRVsqrt10_Exact":    {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFsGRV"+order+"asGRV28_sqrt10.txt",
                                   "Label":      "GRV, 10 GeV²",   
                                   "DataColumn": 2,                 
                                   "XColumn":    0} }
    return output

def DataInitialPDFs(order, namePrint):
    output= { # For Plot 3 & 4
            "SAL3":      {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFsSAL3"+order+"_13.txt",
                          "Label":      namePrint[0]+order+", (1.3)² GeV²",   
                          "DataColumn": 1,                 
                          "XColumn":    0},               
            "SAL4":      {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFsSAL4"+order+"_13.txt",
                          "Label":      namePrint[1]+order+", (1.3)² GeV²",   
                          "DataColumn": 1,                 
                          "XColumn":    0},               
            "SAL5":      {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFsSAL5"+order+"_13.txt",
                          "Label":      namePrint[2]+order+", (1.3)² GeV²",   
                          "DataColumn": 1,                 
                          "XColumn":    0} }
    return output

def DataEvolvedPDFs(name, order, namePrint): 
    output = { # For Plot 5 - 10
            name+"_13":          {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFs"+name+order+"_13.txt",
                                  "Label":      namePrint+order+", (1.3)² GeV²",   
                                  "DataColumn": 1,                 
                                  "XColumn":    0},               
            name+"_sqrt10":      {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFs"+name+order+"_sqrt10.txt",
                                  "Label":      namePrint+order+", 10 GeV²",   
                                  "DataColumn": 1,                 
                                  "XColumn":    0} }
    return output

def DataEvolvedPDFsErrors(name, order, namePrint): 
    output = { # For Plot 5 - 10
            name+"_13":          {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFs"+name+order+"_13.txt",
                                    "Label":      namePrint+order+", (1.3)² GeV²",   
                                    "DataColumn": 4,                 
                                    "XColumn":    0},               
            name+"_sqrt10":      {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFs"+name+order+"_sqrt10.txt",
                                    "Label":      namePrint+order+", 10 GeV²",   
                                    "DataColumn": 4,                 
                                    "XColumn":    0} }
    return output

def DataComparisonPDFs(energySq, namePrint): 
    output = { # For Plot 15 & 16
            namePrint:  {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFsSAL4HO_sqrt"+energySq+".txt",
                         "Label":      namePrint+"HO, "+energySq+" GeV²",   
                         "DataColumn": 1,                 
                         "XColumn":    0},               
            "SAL":      {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFsSAL4HO_sqrt"+energySq+".txt",
                         "Label":      "SAL, "+energySq+" GeV²",   
                         "DataColumn": 2,                 
                         "XColumn":    0},               
            "GRVHO":    {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFsGRVHOasGRV28_sqrt"+energySq+".txt",
                         "Label":      "GRVHO, "+energySq+" GeV²",   
                         "DataColumn": 2,                 
                         "XColumn":    0} }               
    return output

def DataComparisonPDFsErrors(energySq, namePrint): 
    output = { # For Plot 15 & 16
            namePrint:  {"FilePath":   dirBestandsaufnahme+"dataEvolvedPDFsSAL4HO_sqrt"+energySq+".txt",
                         "Label":      namePrint+"HO, "+energySq+" GeV²",   
                         "DataColumn": 4,                 
                         "XColumn":    0} }               
    return output




DataSets01 = {"Data":    DataComparisonEvolutionGRV("LO"), 
              "Errors":  {},
              "pltName": "plotEvolvedPDFsGRVLOasGRV28"}
DataSets02 = {"Data":    DataComparisonEvolutionGRV("HO"), 
              "Errors":  {},
              "pltName": "plotEvolvedPDFsGRVHOasGRV28"}
DataSets03 = {"Data":    DataInitialPDFs("LO", ["Fit3", "Fit4", "Fit5"]), 
              "Errors":  {},
              "pltName": "plotInitialPDFsLO"}
DataSets04 = {"Data":    DataInitialPDFs("HO", ["Fit3", "Fit4", "Fit5"]), 
              "Errors":  {},
              "pltName": "plotInitialPDFsHO"}
DataSets05 = {"Data":    DataEvolvedPDFs("SAL3", "LO", "Fit3"), 
              "Errors":  DataEvolvedPDFsErrors("SAL3", "LO", "Fit3"),
              "pltName": "plotEvolvedPDFsSAL3LOonly"}
DataSets06 = {"Data":    DataEvolvedPDFs("SAL3", "HO", "Fit3"), 
              "Errors":  DataEvolvedPDFsErrors("SAL3", "HO", "Fit3"),
              "pltName": "plotEvolvedPDFsSAL3HOonly"}
DataSets07 = {"Data":    DataEvolvedPDFs("SAL4", "LO", "Fit4"), 
              "Errors":  DataEvolvedPDFsErrors("SAL4", "LO", "Fit4"),
              "pltName": "plotEvolvedPDFsSAL4LOonly"}
DataSets08 = {"Data":    DataEvolvedPDFs("SAL4", "HO", "Fit4"), 
              "Errors":  DataEvolvedPDFsErrors("SAL4", "HO", "Fit4"),
              "pltName": "plotEvolvedPDFsSAL4HOonly"}
DataSets09 = {"Data":    DataEvolvedPDFs("SAL5", "LO", "Fit5"), 
              "Errors":  DataEvolvedPDFsErrors("SAL5", "LO", "Fit5"),
              "pltName": "plotEvolvedPDFsSAL5LOonly"}
DataSets10 = {"Data":    DataEvolvedPDFs("SAL5", "HO", "Fit5"), 
              "Errors":  DataEvolvedPDFsErrors("SAL5", "HO", "Fit5"),
              "pltName": "plotEvolvedPDFsSAL5HOonly"}
DataSets15 = {"Data":    DataComparisonPDFs("2", "Fit4"), 
              "Errors":  DataComparisonPDFsErrors("2", "Fit4"),
              "pltName": "plotEvolvedPDFsSAL4HO_sqrt2"}
DataSets16 = {"Data":    DataComparisonPDFs("20", "Fit4"), 
              "Errors":  DataComparisonPDFsErrors("20", "Fit4"),
              "pltName": "plotEvolvedPDFsSAL4HO_sqrt20"}

DataSetsDict      = {"DataSets01": DataSets01,
                     "DataSets02": DataSets02,
                     "DataSets03": DataSets03,
                     "DataSets04": DataSets04,
                     "DataSets05": DataSets05,
                     "DataSets06": DataSets06,
                     "DataSets07": DataSets07,
                     "DataSets08": DataSets08,
                     "DataSets09": DataSets09,
                     "DataSets10": DataSets10,
                     "DataSets15": DataSets15,
                     "DataSets16": DataSets16}


