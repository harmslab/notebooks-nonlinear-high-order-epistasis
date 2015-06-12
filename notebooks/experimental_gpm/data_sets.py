# This is a major hack ... 
#
# ------------------------------------------------
# Define global variables
# ------------------------------------------------

PATH_TO_DATA = "../../datasets/"

# ------------------------------------------------
# Initial imports
# ------------------------------------------------

import numpy as np
import math


# ------------------------------------------------
# Methods
# ------------------------------------------------

def load_hall_data(filename):
    """ loads data from `filename` and returns genotypes-phenotypes """
    data = np.loadtxt('datasets/' + filename, dtype="S20", delimiter=' ', ndmin=2)
    sequences = data[:,0].astype(str)
    fitness = np.array(data[:,1].astype(str), dtype=float)
    try:
        std = np.array(data[:,2].astype(str), dtype=float)
        return sequences, fitness, std
    except:
        return sequences, fitness
    
def seq_dict(sequences, fitness):
    seq_dict = dict()
    for i in range(len(sequences)):
        if fitness[i] == 0  :
            seq_dict[sequences[i]] = 1e-3
        else:
            seq_dict[sequences[i]] = fitness[i]
    return seq_dict

def sd_dict(sequences, sd, n):
    mean_error = np.mean(sd)
    sd_dict = dict()
    for i in range(len(sequences)):
        if sd[i] == 0  :
            sd_dict[sequences[i]] = mean_error/np.sqrt(n)
        else:
            sd_dict[sequences[i]] = sd[i]/np.sqrt(n)
    return sd_dict, mean_error


sequences1, fitness1, std1 = load_hall_data("hall_haploid_growth.txt")
pheno1 = seq_dict(sequences1,fitness1)
std1_dict, mean_error1 = sd_dict(sequences1,std1,10)
title1 = "Haploid Growth"
sequences2, fitness2, std2 = load_hall_data("hall_diploid_growth.txt")
pheno2 = seq_dict(sequences2,fitness2)
std2_dict, mean_error2 = sd_dict(sequences2,std2,10)
title2 = "Diploid Growth"
sequences3, fitness3, std3 = load_hall_data("hall_sporulation_efficiency.txt")
pheno3 = seq_dict(sequences3,fitness3)
std3_dict, mean_error3 = sd_dict(sequences3,std3,3)
title3 = "Sporulation Efficiency"
sequences4, fitness4, std4 = load_hall_data("hall_dating_efficiency.txt")
pheno4 = seq_dict(sequences4,fitness4)
std4_dict, mean_error4 = sd_dict(sequences4,std4,3)
title4 = "Mating Efficiency"    
