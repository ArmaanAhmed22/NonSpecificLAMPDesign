import pandas as pd
import json
import numpy as np

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
	scored_tm = []
	scored_gc = []
	scored_end = []
	scored_entropy = []
	scored_entropy_end = []
	for region in ["F3","F2","F1c","B1c","B2","B3"]:
		scored_tm+=list(scored_primers[scored_primers.region == region].tm.apply(lambda elem : score_with_ideal_min_max(thermo_params["TM"][region]["minimum"],thermo_params["TM"][region]["maximum"],thermo_params["TM"][region]["ideal"],elem,2.3333)).to_numpy()) #weight results in -0.7 score for 1 degree off from minimum or maximum
		scored_gc+=list(scored_primers[scored_primers.region == region].gc.apply(lambda elem : score_with_ideal_min_max(thermo_params["GC"][region]["minimum"],thermo_params["GC"][region]["maximum"],thermo_params["GC"][region]["ideal"],elem,23.333)).to_numpy()) #weight results in -0.7 score for 10% off from minimum or maximum
	scored_end+=list(scored_primers.end.apply(lambda elem: 1 if elem <= thermo_params["END_STABILITY"]["maximum"] else 0).to_numpy())
	scored_entropy+=list(scored_primers.entropy.apply(lambda elem : 1-(elem / scored_primers.entropy.max())).to_numpy())
	#scored_entropy[region] = np.true_divide(scored_entropy[region],np.amax(scored_entropy[region]))
	scored_entropy_end+=list(scored_primers.entropy_end.apply(lambda elem : 1-(elem / scored_primers.entropy_end.max())).to_numpy())
	total_score = np.array(scored_tm) + np.array(scored_gc) + np.array(scored_end) #+ np.array(scored_entropy)*1/max(scored_entropy) + np.array(scored_entropy_end) * 1/max(scored_entropy_end)

	scored_primers["score"] = total_score
	out = pd.DataFrame(columns=scored_primers.columns)
	for region in ["F3","F2","F1c","B1c","B2","B3"]:
		cur_region = scored_primers[scored_primers.region==region]
		score_min = cur_region["score"].quantile(params["score_quantile"])
		out.append(cur_region[cur_region["score"] >= score_min])
	return out










	"""
	scored_primers = primers.copy()[primers.hairpin == False]
	scored_tm_regions = []
	scored_gc_regions = []
	for region in ["F3","F2","F1c","B1c","B2","B3"]:
		scored_tm_regions.append(scored_primers[scored_primers.region == region].tm.apply(lambda elem : score_with_ideal_min_max(thermo_params["TM"][region]["minimum"],thermo_params["TM"][region]["maximum"],thermo_params["TM"][region]["ideal"],elem)).to_numpy())
		scored_gc_regions.append(scored_primers[scored_primers.region == region].gc.apply(lambda elem : score_with_ideal_min_max(thermo_params["GC"][region]["minimum"],thermo_params["GC"][region]["maximum"],thermo_params["GC"][region]["ideal"],elem)).to_numpy())

	scored_tm = (np.concatenate(tuple(scored_tm_regions)) + 1) / (np.concatenate(tuple(scored_tm_regions)) + 1).max()
	scored_gc = (np.concatenate(tuple(scored_gc_regions)) + 1) / (np.concatenate(tuple(scored_gc_regions)) + 1).max()
	scored_end = scored_primers.end.apply(lambda elem : 1 if elem <= thermo_params["END_STABILITY"]["maximum"] else 0).to_numpy()
	scored_entropy = scored_primers.entropy.apply(lambda elem : 1-(elem / scored_primers.entropy.max())).to_numpy()
	scored_entropy_end = scored_primers.entropy_end.apply(lambda elem : 1-(elem / scored_primers.entropy_end.max())).to_numpy()

	final_score = scored_entropy + scored_entropy_end + scored_tm + scored_gc + scored_end
	scored_primers["score"] = final_score
	true_output = pd.DataFrame(columns=scored_primers.columns)
	for region in ["F3","F2","F1c","B1c","B2","B3"]:
		region_scored_primers = scored_primers[scored_primers["region"] == region]
		score_thresh = region_scored_primers.quantile(params["score_quantile"])["score"]
		true_output = true_output.append(region_scored_primers[region_scored_primers["score"] >= score_thresh])
	return true_output"""


primers = pd.read_csv(snakemake.input[0])
thermo_params = json.load(open(snakemake.input[1],"r"))
params = json.load(open(snakemake.input[2],"r"))
output = score_primers()
output.to_csv(snakemake.output[0],index=False)