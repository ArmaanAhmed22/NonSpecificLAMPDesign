import pandas as pd
from Bio.Seq import Seq
import json

def main():
	parts_dict = {primer_region : primer_parts[f"{primer_region}_sequence"] for primer_region in ["F3","F2","F1c","B1c","B2","B3","LoopF","LoopB"]}

	output_json = {}
	output_json["FILE_NAME"] = snakemake.input[0]
	output_json["F3_PRIMER"] = parts_dict["F3"]
	output_json["FIP"] = str(Seq(parts_dict["F1c"]).reverse_complement()) + parts_dict["F2"]
	output_json["BIP"] = parts_dict["B1c"] + str(Seq(parts_dict["B2"]).reverse_complement())
	output_json["B3_PRIMER"] = str(Seq(parts_dict["B3"]).reverse_complement())
	output_json["LOOP_F"] = str(Seq(parts_dict["LoopF"]).reverse_complement())
	output_json["LOOP_B"] = str(Seq(parts_dict["LoopB"]))

	return output_json

primer_parts = pd.read_csv(snakemake.input[0]).iloc[0]

json.dump(main(),open(snakemake.output[0],"w"), indent=4)
