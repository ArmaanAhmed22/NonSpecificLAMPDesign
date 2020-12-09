from Bio import SeqIO
import json
from tqdm import tqdm
import pandas as pd
import numpy as np

#Data from https://primerexplorer.jp/e/v5_manual/
params = {
	"length": 
	{
		"F3":
		{
			"minimum":15,
			"maximum":20
		},
		"F2":
		{
			"minimum":15,
			"maximum":20
		},
		"F1c":
		{
			"minimum":15,
			"maximum":22
		},
		"B1c":
		{
			"minimum":15,
			"maximum":22
		},
		"B2":
		{
			"minimum":15,
			"maximum":20
		},
		"B3":
		{
			"minimum":15,
			"maximum":20
		}
	},
	"spacing":
	{
		"F3|3',F2|5'":
		{
			"minimum":0,
			"maximum":60
		},
		"F2|5',F1c|3'":
		{
			"minimum":40,
			"maximum":60
		},
		"B1c|3',B2|5'":
		{
			"minimum":40,
			"maximum":60
		},
		"F2|5',B2|5'":
		{
			"minimum":120,
			"maximum":160
		},
		"B2|5',B3|3'":
		{
			"minimum":0,
			"maximum":60
		}
	},
	"mismatch_tolerance":
	{
		"intolerance":
		{
			"F3":[1],
			"F2":[1],
			"F1c":[-1],
			"B1c":[-1],
			"B2":[1],
			"B3":[1]
		},
		"tolerance":1
	}

}

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

def get_primer_locs(start,primer_lengths,primer_spacing):
	primer_locs = {}

	primer_locs["F3"] = start
	primer_locs["F2"] = primer_locs["F3"] + primer_lengths["F3"] + primer_spacing["F3|3',F2|5'"]
	primer_locs["F1c"] = primer_locs["F2"] + primer_lengths["F2"] + (primer_spacing["F2|5',F1c|3'"] - primer_lengths["F2"])

	primer_locs["B2"] = primer_locs["F2"] + primer_lengths["F2"] + (primer_spacing["F2|5',B2|5'"] - primer_lengths["B2"] - primer_lengths["F2"])
	primer_locs["B1c"] = primer_locs["B2"] + primer_lengths["B2"] - (primer_spacing["B1c|3',B2|5'"] - primer_lengths["B1c"])
	primer_locs["B3"] = primer_locs["B2"]+primer_lengths["B2"]+primer_spacing["B2|5',B3|3'"]

	return primer_locs



def individual_primer_generate(reference,primer_lengths,primer_spacing,nuc_prevalence):
	#df = pd.DataFrame(columns=["region","pos","length","entropy"])

	def generate_region(region,downstream_requirements,upstream_requirements,depend_on_length=True):
		data = []
		for region_len in tqdm(range(primer_lengths[region]["minimum"],primer_lengths[region]["maximum"]+1),desc=region):
			region_pos = downstream_requirements
			while (region_pos + ((region_len if depend_on_length else 0)+ upstream_requirements) <= len(reference)):
				entropy = nuc_prevalence["shannon"][region_pos:region_pos+region_len].sum()
				#df = df.append(,ignore_index=True)
				data.append({"region":region,"pos":region_pos,"length":region_len,"entropy":entropy})
				region_pos+=1
		return data

	data_hold = []

	#Generate F3 primer data
	data_hold.extend(
	generate_region("F3",
		0,
		primer_spacing["F3|3',F2|5'"]["minimum"] + primer_spacing["F2|5',B2|5'"]["minimum"] + primer_spacing["B2|5',B3|3'"]["minimum"] + primer_lengths["B3"]["minimum"])
	)

	#Generate F2 primer data
	data_hold.extend(
	generate_region("F2",
		primer_lengths["F3"]["minimum"] + primer_spacing["F3|3',F2|5'"]["minimum"],
		primer_spacing["F2|5',B2|5'"]["minimum"] + primer_spacing["B2|5',B3|3'"]["minimum"] + primer_lengths["B3"]["minimum"],
		depend_on_length=False)
	)

	#Generate F1c primer data
	data_hold.extend(
	generate_region("F1c",
		primer_lengths["F3"]["minimum"]+primer_spacing["F3|3',F2|5'"]["minimum"]+primer_spacing["F2|5',F1c|3'"]["minimum"],
		-primer_spacing["F2|5',F1c|3'"]["maximum"] + primer_spacing["F2|5',B2|5'"]["minimum"] + primer_spacing["B2|5',B3|3'"]["minimum"] + primer_lengths["B3"]["minimum"],
		depend_on_length=False
		)
	)

	#Generate B1c primer data
	data_hold.extend(
	generate_region("B1c",
		primer_lengths["F3"]["minimum"] + primer_spacing["F3|3',F2|5'"]["minimum"] + primer_spacing["F2|5',B2|5'"]["minimum"] - primer_spacing["B1c|3',B2|5'"]["maximum"] - primer_lengths["B1c"]["maximum"],
		primer_spacing["B1c|3',B2|5'"]["minimum"]+primer_spacing["B2|5',B3|3'"]["minimum"] + primer_lengths["B3"]["minimum"]
		)
	)

	#Generate B2 primer data
	data_hold.extend(
	generate_region("B2",
		primer_lengths["F3"]["minimum"] + primer_spacing["F3|3',F2|5'"]["minimum"] + primer_spacing["F2|5',B2|5'"]["minimum"] - primer_lengths["B2"]["maximum"],
		primer_spacing["B2|5',B3|3'"]["minimum"] + primer_lengths["B3"]["minimum"]
		)
	)
	# Generate B3 primer data
	data_hold.extend(
	generate_region("B3",
		primer_lengths["F3"]["minimum"] + primer_spacing["F3|3',F2|5'"]["minimum"] + primer_spacing["F2|5',B2|5'"]["minimum"] + primer_spacing["B2|5',B3|3'"]["minimum"],
		0
		)
	)
	result = pd.DataFrame(data_hold)
	result["score"] = result["entropy"].apply(lambda x: x - result["entropy"].min()).apply(lambda x: -x / (result["entropy"].max() - result["entropy"].min()))
	print(result.head())
	return result


#open reference
reference = next(SeqIO.parse(open(snakemake.input[0],"r"),"fasta")).seq

#open nuc prevalence to calculate shannon entropy
nuc_prevalence = pd.read_csv(snakemake.input[1])

#Not working for some reason...
"""
#load parameters
with open(snakemake.input[2],"r") as handle:
	params = json.loads(str(handle.read()),encoding="UTF-8")
"""

calc_base_pair_shannon(nuc_prevalence)

df = individual_primer_generate(reference,params["length"],params["spacing"],nuc_prevalence)
df.to_csv(snakemake.output[0],index=False)