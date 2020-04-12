import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import os
import numpy as np
import matplotlib.pyplot as plt
from nmf_methods import *

def subsample(n_samples):
        files_to_subsample = "data/raw/250_44_features.tsv"
        output_locations = "data/raw/250_{}_44_features.tsv".format(n_samples)

        f_df = pd.read_csv(files_to_subsample, sep="\t", index_col=0)
        f_df.iloc[0:n_samples].to_csv(output_locations, sep="\t")

def create_data(n_samples):
	command = """python sim_dmr_data.py -tf data/COSMIC/cosmic-signatures.tsv -ff  data/raw/{total_n_samples}_{num_samples}_{iter_num}_features.tsv -l  data/raw/lambdas.csv -enm data/raw/BRCA-nmutations.tsv -at "1" "2" "3" "5" -n {total_n_samples} -od data/processed/{total_n_samples}_{num_samples}_{iter_num}_documents.tsv  -om data/processed/{total_n_samples}_{num_samples}_{iter_num}_model.npz  -omf data/processed/{total_n_samples}_{num_samples}_{iter_num}_mutation_count.tsv -os data/processed/{total_n_samples}_{num_samples}_{iter_num}_signatures.tsv -oe data/processed/{total_n_samples}_{num_samples}_{iter_num}_exposures.tsv -ni {iter_num}"""
	command = command.replace("{iter_num}","44")
	command = command.replace("{total_n_samples}","250")
	command = command.replace("{num_samples}",str(n_samples))	
	os.system(command)

def best_similarity(n_samples,nmf_function):
	cosmic_location = "data/COSMIC/cosmic-signatures.tsv"
	cosmic_data = pd.read_csv(cosmic_location,sep="\t").values[:,1:]
	our_data = nmf_function(n_samples)

	similarity_matrix = cosine_similarity(cosmic_data,our_data)

	return np.max(similarity_matrix[2]),np.max(similarity_matrix[4])

def full_pipeline(n_samples,nmf_function):
	subsample(n_samples)
	create_data(n_samples)
	return best_similarity(n_samples,nmf_function)

def run_one_trial(nmf_function):
	return [full_pipeline(i,nmf_function) for i in [50,100,150,200,250]]


def run_multiple_trials(num_trial,nmf_function):
	all_trials = [[] for i in [50,100,150,200,250]]
	for i in range(num_trial):
		print("On trial {}".format(i))
		current_trial = run_one_trial(nmf_function)
		for j in range(len(current_trial)):
			all_trials[j].append(current_trial[j])
	return all_trials

def plot_results(trials,nmf_function,nmf_name):
	data = run_multiple_trials(trials,nmf_function)
	signatures = [3,5]
	for i in range(2):
		temp_data = [[j[i] for j in k] for k in data]
		plt.boxplot(temp_data)
		plt.title("Signature number {}".format(signatures[i]))
		plt.xticks(list(range(1,len(data)+1)),["50","100","150","200","250"])
		plt.ylabel("Cosine similarity")
		plt.savefig("plots/signature_{}_{}.pdf".format(signatures[i],nmf_name))
		plt.clf()

def plot_different_methods(trials,function_list,name_list):
	data_list = [run_multiple_trials(trials,i) for i in function_list]
	signatures = [3,5]
	for i in range(2):
		for k in range(len(function_list)):
			temp_data = [[j[i] for j in r] for r in data_list[k]]
			temp_data = [np.mean(np.array(j)) for j in temp_data]
			plt.plot([50,100,150,200,250],temp_data,label=name_list[k])
		plt.ylabel("Cosine Similarity")
		plt.legend()
		plt.savefig("plots/compare_signature_{}.pdf".format(signatures[i]))
		plt.clf()

def curry_netnmf(alpha):
	def f(x):
		return netnmf_method(x,alpha=alpha)
	return f

nmf_functions = [netnmf_method,factor_SS,sklearn_method]
nmf_names = ["MU_netnmf","SS","Sklearn"]

results = []
alpha_values = [1,1.5,2,2.5,3,3.5,4]
for value in alpha_values:
	results.append((value,run_multiple_trials(5,curry_netnmf(value))))

for i in results:
	print(i)
