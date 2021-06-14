import pandas as pd
import json
import numpy as np
import yaml

def score_with_ideal_min_max(minimum,maximum,ideal,elem,weight=1) -> float:
	if elem < minimum:
		elem_m = minimum - elem
		return 1/(weight*elem_m+1) - 1
	elif elem > maximum:
		elem_m = elem - maximum
		return 1/(weight*elem_m+1) - 1
	if elem < ideal:
		return 1/(ideal-minimum) * (elem - minimum)
	else:
		return 1/(ideal-maximum) * (elem - maximum)

def score_primers():
	scored_primers = df.copy()[df.hairpin == False]
	scores = {
	"tm":[],
	"gc":[],
	"end":[],
	"entropy":[],
	"entropy_end":[]
	}
	for region in dict.fromkeys(scored_primers.region).keys():
		scores["tm"]+=list(scored_primers[scored_primers.region == region].tm.apply(lambda elem : score_with_ideal_min_max(thermo_params["TM"][region]["minimum"],thermo_params["TM"][region]["maximum"],thermo_params["TM"][region]["ideal"],elem,2.3333)).to_numpy()) #weight results in -0.7 score for 1 degree off from minimum or maximum
		scores["gc"]+=list(scored_primers[scored_primers.region == region].gc.apply(lambda elem : score_with_ideal_min_max(thermo_params["GC"][region]["minimum"],thermo_params["GC"][region]["maximum"],thermo_params["GC"][region]["ideal"],elem,23.333)).to_numpy()) #weight results in -0.7 score for 10% off from minimum or maximum
	scores["end"]+=list(scored_primers.end.apply(lambda elem: 1 if elem <= thermo_params["END_STABILITY"]["maximum"] else 0).to_numpy())
	scores["entropy"]+=list(scored_primers.entropy.apply(lambda elem : 1-(elem / scored_primers.entropy.max())).to_numpy())
	#scored_entropy[region] = np.true_divide(scored_entropy[region],np.amax(scored_entropy[region]))
	scores["entropy_end"]+=list(scored_primers.entropy_end.apply(lambda elem : 1-(elem / scored_primers.entropy_end.max())).to_numpy())
	#total_score = (np.array(scored_tm) + np.array(scored_gc) + np.array(scored_end) + np.array(scored_entropy) + np.array(scored_entropy_end)) / 5
	total_score = returned_total_score(scored_primers["region"].to_numpy(),{k:np.array(v) for k,v in scores.items()})
	scored_primers["score"] = total_score
	out = pd.DataFrame(columns=scored_primers.columns)
	for region in ["F3","F2","F1c","B1c","B2","B3","LoopF","LoopB"]:
		cur_region = scored_primers[scored_primers.region==region]
		score_min = cur_region["score"].quantile(params["score_quantile"])
		out = out.append(cur_region[cur_region["score"] >= score_min])
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
