import os
import glob
import pandas as pd
import numpy as np

def create_data():    
    os.system("cd ../tcsm/experiments/simulated-data && snakemake -R simulate_dmr_data --cores 1")
    data_list = glob.glob("../tcsm/experiments/simulated-data/processed/*mutation_count.tsv")
    d = {}
    for i in data_list:
        num_samples = i.split("_")[0]
        d[int(num_samples)] = i

    return d

def get_signatures_3_5():
    df = pd.read_csv("../tcsm/experiments/simulated-data/COSMIC/cosmic-signatures.tsv",sep="\t")
    all_signatures = df.values
    return all_signatures[[2,4],:]

def factor_somatic(data_location):
    # Clear out any current signatures
    os.system("rm -rf ../tcsm/SomaticSignatures/processed/*")

    # Create the new signatures
    os.system("cp {} ../tcsm/SomaticSignatures/processed/mutation-count-all_TCGA-BRCA_0.tsv".format(data_location))
    os.system("cd ../tcsm/SomaticSignatures && snakemake -R run_somatic_signatures_all --cores 1")

    # Read in the signatures
    file_name = glob.glob("../tcsm/SomaticSignatures/processed/*signatures*.tsv")[0]
    df = pd.read_csv(file_name, sep="\t")  

    # Returns a 4x96 numpy array 
    return df.values

def factor_tcsm(data_location):
    os.system("../tcsm/src/run_tcsm.R {} 4".format(data_location))

    return pd.read_csv("../tcsm/signatures.tsv",sep="\t")

def factor_tcsm_with_HR_covariate(data_location):    
    covariate_location = data_location.replace("mutation_count","features")
    os.system("../tcsm/src/run_tcsm.R {} 4 -c={}".format(data_location,covariate_location))

    return pd.read_csv("../tcsm/signatures.tsv",sep="\t")




        
                     
        
