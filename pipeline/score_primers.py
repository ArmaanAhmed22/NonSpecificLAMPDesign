import pandas as pd
import json
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import yaml

def score_error(x,weight,y_int):
	return 1 / (weight * x + 1/y_int)

def score_good(x_normalized,y_int):
	return -(1 - y_int) * x_normalized + 1

def normalize(x, minimum, maximum):
	return (x - minimum) / (maximum - minimum)

def score_with_ideal_min_max(minimum,maximum,ideal,elem,weight=1,y_int=1) -> float:
	"""if elem < minimum:
		elem_m = minimum - elem
		return 1/(weight*elem_m+1) - 1
	elif elem > maximum:
		elem_m = elem - maximum
		return 1/(weight*elem_m+1) - 1
	if elem < ideal:
		return 1/(ideal-minimum) * (elem - minimum)
	else:
		return 1/(ideal-maximum) * (elem - maximum)"""
	if elem < minimum or elem > maximum:
		if elem < minimum:
			x = minimum - elem
		else:
			x = elem - maximum
		return score_error(x,weight,y_int)
	else:
		if elem < ideal:
			x_norm = 1 - normalize(elem, minimum, ideal)
		else:
			x_norm = normalize(elem, ideal, maximum)
		return score_good(x_norm,y_int)
		

def score_primers():
	scored_primers = df.copy()[df.hairpin == False]
	scores = {
	"tm":np.array([]),
	"gc":np.array([]),
	"end":scored_primers.end.apply(lambda elem: 1 if elem <= thermo_params["END_STABILITY"]["maximum"] else 0).to_numpy(),
	"entropy":scored_primers.entropy.apply(lambda elem : 1-(elem / scored_primers.entropy.max())).to_numpy(),
	"entropy_end":scored_primers.entropy_end.apply(lambda elem : 1-(elem / scored_primers.entropy_end.max())).to_numpy()
	}
	#Scoring: https://www.desmos.com/calculator/haitgwoprg
	for region in dict.fromkeys(scored_primers.region).keys():
		cur_tms = scored_primers[scored_primers.region == region].tm.apply(lambda elem : score_with_ideal_min_max(thermo_params["TM"][region]["minimum"],thermo_params["TM"][region]["maximum"],thermo_params["TM"][region]["ideal"],elem,weight=8.75,y_int=0.9)).to_numpy() #weight results in 0.1 score for 1 degree off from minimum or maximum
		cur_gcs = scored_primers[scored_primers.region == region].gc.apply(lambda elem : score_with_ideal_min_max(thermo_params["GC"][region]["minimum"],thermo_params["GC"][region]["maximum"],thermo_params["GC"][region]["ideal"],elem,weight=41,y_int=0.8)).to_numpy() #weight results in 70% decrease for 5% off from minimum or maximum
		scores["tm"] = np.append(scores["tm"],cur_tms)
		scores["gc"] = np.append(scores["gc"],cur_gcs)

	total_score = returned_total_score(scored_primers["region"].to_numpy(),{k:v for k,v in scores.items()})
	scored_primers["score"] = total_score
	out = pd.DataFrame(columns=scored_primers.columns)
	for region in ["F3","F2","F1c","B1c","B2","B3","LoopF","LoopB"]:
		cur_region = scored_primers[scored_primers.region==region]
		score_min = cur_region["score"].quantile(params["score_quantile"])
		out = out.append(cur_region[cur_region["score"] >= score_min])
	fig,ax = plt.subplots(nrows=2)
	sns.lineplot(x="order",y="score",data=scored_primers.sort_values(by="score",ascending=False).assign(order=range(1,scored_primers.shape[0]+1)),ax=ax[0],hue="region")
	sns.histplot(x="score",data=scored_primers,ax=ax[1],bins=10)
	sns.despine(ax=ax[0])
	sns.despine(ax=ax[1])
	plt.savefig(snakemake.output[1],bbox_inches="tight",dpi=600)
	#print(f'{scored_primers.sort_values(by="score",ascending=False).head(n=20)}')
	return out

def get_mean_method(mean_method):
	if mean_method == "arithmetic":
		def arithmetic_mean(scores,W):
			length = len(list(scores.values())[0])
			out = np.zeros(length)
			for k in scores.keys():
				out+=scores[k] * W[k]
			return out / sum(W.values())
		return arithmetic_mean

	elif mean_method == "geometric":
		def geometric_mean(scores,W):
			length = len(list(scores.values())[0])
			out = np.ones(length)
			for k in scores.keys():
				out*=scores[k]**W[k]
			return out ** (1 / sum(W.values()))
		return geometric_mean
	elif mean_method == "harmonic":
		def harmonic_mean(scores,W):
			length = len(list(scores.values())[0])
			out = np.zeros(length)
			for k in scores.keys():
				out+= W[k] / scores[k]
			return sum(W.values()) / out
		return harmonic_mean
	else:
		raise Exception("No mean method")

def returned_total_score(regions,scores):
	out = np.zeros(len(regions))
	mean = get_mean_method(score_params["score_mean"])
	regions_set = list(dict.fromkeys(regions).keys())
	end_starts = {region : [-1,-1] for region in regions_set}
	for i,k in enumerate(end_starts.keys()):
		end_starts[k][0] = regions.tolist().index(k)
		if i+1 >= len(regions_set):
			end_starts[k][1] = len(regions)
		else:
			end_starts[k][1] = regions.tolist().index(regions_set[i+1])
	for region in regions_set:
		W = score_params["weights_score"][region]
		cur_scores = {k:None for k in scores.keys()}
		for k in scores.keys():
			cur_scores[k] = scores[k][end_starts[region][0] : end_starts[region][1]]
		out[end_starts[region][0] : end_starts[region][1]] = mean(cur_scores,W)
	return out

df = pd.read_csv(snakemake.input[0])
thermo_params = json.load(open(snakemake.input[1],"r"))
params = json.load(open(snakemake.input[2],"r"))
score_params = yaml.safe_load(open(snakemake.params[0]))

output = score_primers()
output.to_csv(snakemake.output[0],index=False)
