import pandas as pd
import json
import numpy as np

def score_with_ideal_min_max(minimum,maximum,ideal,elem):
	if elem < minimum or elem > maximum:
		return -1
	max_distance = max(abs(ideal-minimum),abs(ideal-maximum))
	distance = abs(ideal - elem)
	return max_distance - distance

def score_primers():
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
	return true_output


primers = pd.read_csv(snakemake.input[0])
thermo_params = json.load(open(snakemake.input[1],"r"))
params = json.load(open(snakemake.input[2],"r"))
output = score_primers()
output.to_csv(snakemake.output[0],index=False)