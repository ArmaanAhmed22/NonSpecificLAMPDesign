import pandas as pd
from Bio.Seq import Seq
import json






def main():
	parts_dict = {num_to_region[i] : primer_parts[f"{i}_sequence"] for i in range(6)}

	output_json = {}
	output_json["FILE_NAME"] = snakemake.input[0]
	output_json["F3_PRIMER"] = parts_dict["F3"]
	output_json["FIP"] = str(Seq(parts_dict["F1c"]).reverse_complement()) + parts_dict["F2"]
	output_json["BIP"] = parts_dict["B1c"] + str(Seq(parts_dict["B2"]).reverse_complement())
	output_json["B3_PRIMER"] = str(Seq(parts_dict["B3"]).reverse_complement())

	return output_json







num_to_region = {
	0:"F3",
	1:"F2",
	2:"F1c",
	3:"B1c",
	4:"B2",
	5:"B3"
}
primer_parts = pd.read_csv(snakemake.input[0]).iloc[0]

json.dump(main(),open(snakemake.output[0],"w"), indent=4)
