import json
import pandas as pd
import numpy as np

primer_sets = json.load(open(snakemake.input[0],"r"))
nuc_prevalence = pd.read_csv(snakemake.input[1])

num_degenerates = json.load(open(snakemake.input[2],"r"))["degenerates"]
ambiguous_to_nt = json.load(open(snakemake.input[3],"r"))

def calc_base_pair_shannon(nuc_prevalence):
	nuc_prevalence["total"] = nuc_prevalence["A"] + nuc_prevalence["G"] + nuc_prevalence["C"] + nuc_prevalence["T"]
	shannon = None
	for nuc in ["A","G","C","T"]:
		percent = nuc_prevalence[nuc] / nuc_prevalence["total"]
		if shannon is None:
			shannon = np.where(percent == 0,percent,-percent * np.log2(percent))
		else:
			shannon += np.where(percent == 0,percent,-percent * np.log2(percent))
	nuc_prevalence["shannon"] = shannon
def get_ambiguous(A,G,C,T,cutoff=0.1):
	total = A+G+C+T
	possible = []
	for nuc,prevalence in {"A":A,"G":G,"C":C,"T":T}.items():
		if prevalence / total >= cutoff:
			possible.append(nuc)
	if len(possible) == 0:
		#print("yes")
		return None
	if len(possible) == 1:
		return possible[0]
	for amb,nucs in ambiguous_to_nt.items():
		should_use_this_amb = True
		for cur_nuc in possible:
			if not (cur_nuc in nucs):
				should_use_this_amb = False
				break
		if not should_use_this_amb:
			continue
		return amb
	return None


calc_base_pair_shannon(nuc_prevalence)
degenerate_primer_sets_proto = [{f"{primer_num}_{item}":None for item in ["sequence","position"] for primer_num in range(6)} for _ in range(len(primer_sets))]

for primer_set_number,primer_set in enumerate(primer_sets.values()):
	for primer_num in range(6):
		least_conserved = [{"index":-1,"shannon":-1} for _ in range(num_degenerates)]
		most_conserved_index,most_conserved = (-1,{"index":-1,"shannon":-1})
		sequence = ""
		for npos in range(primer_set["set_position"][primer_num],primer_set["set_position"][primer_num]+primer_set["set_length"][primer_num]):
			##Get least conserved reginos for degenerate insertions
			if nuc_prevalence.iloc[npos]["shannon"] > most_conserved["shannon"]:
				least_conserved[most_conserved_index]["index"] = npos
				least_conserved[most_conserved_index]["shannon"] = nuc_prevalence.iloc[npos]["shannon"]

				most_conserved_index,most_conserved = min(enumerate(least_conserved),key=lambda x: x[1]["shannon"])
		for npos in range(primer_set["set_position"][primer_num],primer_set["set_position"][primer_num]+primer_set["set_length"][primer_num]):
			if npos in [nuc["index"] for nuc in least_conserved]:
				amb = get_ambiguous(nuc_prevalence.iloc[npos]["A"],nuc_prevalence.iloc[npos]["G"],nuc_prevalence.iloc[npos]["C"],nuc_prevalence.iloc[npos]["T"])
				sequence+=amb
			else:
				max_nuc,prevalence = max([(nuc,nuc_prevalence.iloc[npos][nuc]) for nuc in ["A","G","C","T"]],key=lambda x:x[1])
				sequence+=max_nuc
		degenerate_primer_sets_proto[primer_set_number][f"{primer_num}_sequence"] = sequence
		degenerate_primer_sets_proto[primer_set_number][f"{primer_num}_position"] = primer_set["set_position"][primer_num]

pd.DataFrame(degenerate_primer_sets_proto).to_csv(snakemake.output[0],index=False)