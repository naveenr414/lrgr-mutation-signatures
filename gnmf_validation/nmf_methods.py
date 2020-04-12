import pandas as pd
import os
import numpy as np
import tensorflow as tf
from sklearn.decomposition import NMF
import numpy as np
import othercode

def factor_SS(n_samples):
	data_location = "data/processed/250_{}_44_mutation_count.tsv".format(n_samples)
	signature_location = "data/processed/250_{}_44_SS_signatures.tsv".format(n_samples)
	exposure_location = "data/processed/250_{}_44_SS_exposures.tsv".format(n_samples)
	os.system("Rscript othercode/run_somatic_signatures.R {} {} {}".format(data_location,signature_location,exposure_location))

	df = pd.read_csv(signature_location,sep="\t").values
	return df

def kl_divergence(p, q): 
    return tf.reduce_sum((p+.00001) * tf.log((p+.0001)/(q+.0001)))

def gradient_descent(n_samples,lr=.0001,trials=8000):
	tf.reset_default_graph()
	data_location = "data/processed/250_{}_44_mutation_count.tsv".format(n_samples)
	data = pd.read_csv(data_location,sep="\t").values[:,1:]

	V_tf = tf.constant(data,dtype=tf.float32)
	W = np.random.random((n_samples,4))*0
	H = np.random.random((4,96))*10
	W_tf = tf.Variable(W,dtype=tf.float32)
	H_tf = tf.Variable(H,dtype=tf.float32)

	cost = kl_divergence(tf.matmul(W_tf,H_tf),V_tf)
	train_step = tf.train.AdamOptimizer(lr).minimize(cost)	

	init = tf.global_variables_initializer()

	clipW = W_tf.assign(tf.maximum(tf.zeros_like(W_tf),W_tf))
	clipH = H_tf.assign(tf.maximum(tf.zeros_like(H_tf),H_tf))
	clip = tf.group(clipW,clipH)

	with tf.Session() as sess:
		sess.run(init)
		for i in range(trials):
			sess.run(train_step)
			sess.run(clip)
		return sess.run(H_tf)

def sklearn_method(n_samples):
	data_location = "data/processed/250_{}_44_mutation_count.tsv".format(n_samples)
	data = pd.read_csv(data_location,sep="\t").values[:,1:]
	model = NMF(n_components=4,init='random',solver='mu',beta_loss='kullback-leibler')
	model.fit_transform(data)
	return model.components_

def get_graph(n_samples):
	data_location = "data/raw/250_{}_44_features.tsv".format(n_samples)
	data = pd.read_csv(data_location,sep="\t").values[:,1]
	return data.reshape(len(data),1).dot(data.T.reshape(1,len(data)))+.0001
                                                                   
def netnmf_method(n_samples,alpha=0):
	dimensions = 4
	max_iters = 100000
	n_jobs = 1

	operator = othercode.netNMFMU(n_components=4,init='random',solver='mu',beta_loss='KL',l1_ratio=0.0, alpha=alpha,tol=0.0001, max_iter=20000)
	data_location = "data/processed/250_{}_44_mutation_count.tsv".format(n_samples)
	data = pd.read_csv(data_location,sep="\t").values[:,1:]
	operator.X = data

	operator.N = get_graph(data.shape[0])
	print(np.sum(operator.N))
	operator.fit_transform(alpha,operator.N,operator.X)
	return operator.components_
