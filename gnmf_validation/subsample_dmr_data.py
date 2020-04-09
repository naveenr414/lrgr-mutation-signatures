import pandas as pd
import sys


def subsample(n_samples):
	files_to_subsample = "raw/250_44_features.tsv"
	output_locations = "raw/250_{}_44_features.tsv".format(n_samples)

	f_df = pd.read_csv(files_to_subsample, sep="\t", index_col=0)
	f_df.iloc[0:n_samples].to_csv(output_locations, sep="\t")

