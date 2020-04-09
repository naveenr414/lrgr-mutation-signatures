import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import os
import numpy as np
import matplotlib.pyplot as plt

def subsample(n_samples):
        files_to_subsample = "raw/250_44_features.tsv"
        output_locations = "raw/250_{}_44_features.tsv".format(n_samples)

        f_df = pd.read_csv(files_to_subsample, sep="\t", index_col=0)
        f_df.iloc[0:n_samples].to_csv(output_locations, sep="\t")

def create_data(n_samples):
	command = """python sim_dmr_data.py -tf COSMIC/cosmic-signatures.tsv -ff  raw/{total_n_samples}_{num_samples}_{iter_num}_features.tsv -l  raw/lambdas.csv -enm raw/BRCA-nmutations.tsv -at "1" "2" "3" "5" -n {total_n_samples} -od processed/{total_n_samples}_{num_samples}_{iter_num}_documents.tsv  -om processed/{total_n_samples}_{num_samples}_{iter_num}_model.npz  -omf processed/{total_n_samples}_{num_samples}_{iter_num}_mutation_count.tsv -os processed/{total_n_samples}_{num_samples}_{iter_num}_signatures.tsv -oe processed/{total_n_samples}_{num_samples}_{iter_num}_exposures.tsv -ni {iter_num}"""
	command = command.replace("{iter_num}","44")
	command = command.replace("{total_n_samples}","250")
	command = command.replace("{num_samples}",str(n_samples))	
	os.system(command)

def factor_SS(n_samples):
	data_location = "processed/250_{}_44_mutation_count.tsv".format(n_samples)
	signature_location = "processed/250_{}_44_SS_signatures.tsv".format(n_samples)
	exposure_location = "processed/250_{}_44_SS_exposures.tsv".format(n_samples)
	os.system("Rscript run_somatic_signatures.R {} {} {}".format(data_location,signature_location,exposure_location))

def best_similarity(n_samples):
	cosmic_location = "COSMIC/cosmic-signatures.tsv"
	our_location = "processed/250_{}_44_SS_signatures.tsv".format(n_samples)

	cosmic_data = pd.read_csv(cosmic_location,sep="\t").values[:,1:]
	our_data = pd.read_csv(our_location,sep="\t").values

	similarity_matrix = cosine_similarity(cosmic_data,our_data)

	return np.max(similarity_matrix[2]),np.max(similarity_matrix[4])

def full_pipeline(n_samples):
	subsample(n_samples)
	create_data(n_samples)
	factor_SS(n_samples)
	return best_similarity(n_samples)

def run_one_trial():
	return [full_pipeline(i) for i in [50,100,150,200,250]]


def run_multiple_trials(num_trial):
	all_trials = [[] for i in [50,100,150,200,250]]
	for i in range(num_trial):
		print("On trial {}".format(i))
		current_trial = run_one_trial()
		for j in range(len(current_trial)):
			all_trials[j].append(current_trial[j])
	return all_trials

data = run_multiple_trials(100)
signatures = [3,5]
for i in range(2):
	temp_data = [[j[i] for j in k] for k in data]
	plt.boxplot(temp_data)
	plt.title("Signature number {}".format(signatures[i]))
	plt.xticks(list(range(1,len(data)+1)),["50","100","150","200","250"])
	plt.ylabel("Cosine similarity")
	plt.savefig("plots/signature_{}.pdf".format(signatures[i]))
	plt.clf()
