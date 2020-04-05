from get_signatures import *
from scipy.spatial.distance import cosine 
import matplotlib.pyplot as plt

def plot_signatures(signature_number):
    all_locations = create_data()
    methods = {'SomaticSignatures':factor_somatic,'TCSM':factor_tcsm,'TCSM+Covariate':factor_tcsm_with_HR_covariate}
    true_signatures = get_signatures_3_5()

    for i in methods:
        f = methods[i]
        x_loc = []
        y_loc = []

        for loc in all_locations:
            factored_signatures = f(loc)
            num_samples = int(loc.split("/")[-1].split("_")[0])

            pair_similarities = []

            for i in range(len(factored_signatures)):
                for j in range(len(factored_signatures)):
                    if i!=j:
                        signature_3_similarity = 1-cosine(factored_signatures[i],true_signatures[0])
                        signature_5_similarity = 1-cosine(factored_signatures[j],true_signatures[1])

                        pair_similarities[(i,j)] = signature_3_similarity+signature_5_similarity
                        
            i,j = max(pair_similarities)
            x_loc.append(num_samples)
            if(signature_number == 3):
                y_loc.append(1-cosine(factored_signatures[i],true_signatures[0]))
            else:
                y_loc.append(1-cosine(factored_signatures[j],true_signatures[1]))

        plt.plot(x,y,label=i)

    plt.legend()
    plt.savefig("signature_{}_similarity.jpg".format(signature_number))

plot_signatures(3)

