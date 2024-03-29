import pandas as pd
from tqdm import tqdm
import numpy as np
import json

params = json.load(open(snakemake.input[1],"r"))

def get_top_unit():
	return {"set_position":[-1,-1,-1,-1,-1,-1],"set_length":[-1,-1,-1,-1,-1,-1],"score":[-1,-1,-1,-1,-1,-1]}
top = [get_top_unit() for _ in range(params["top"])]
#top = {"set":[-1,-1,-1,-1,-1,-1],"score":[-1,-1,-1,-1,-1,-1]}

df_all_primers = pd.read_csv(snakemake.input[0])
df_all_primers["max_before_score"] = np.array([-100]*df_all_primers.shape[0])

sort_regions_by = lambda region,by : df_all_primers[df_all_primers["region"] == region].sort_values(by=by,ascending=False)
df_F3 = df_all_primers[df_all_primers["region"] == "F3"]
df_F2 = df_all_primers[df_all_primers["region"] == "F2"]
df_F1c = df_all_primers[df_all_primers["region"] == "F1c"]
df_B1c = df_all_primers[df_all_primers["region"] == "B1c"]
df_B2 = df_all_primers[df_all_primers["region"] == "B2"]
df_B3 = df_all_primers[df_all_primers["region"] == "B3"]

def get_satisfied_spacing(spacing,region):
	return not (spacing > params["spacing"][region]["maximum"] or spacing < params["spacing"][region]["minimum"])

for f3_data in tqdm(df_F3.itertuples(),total=df_F3.shape[0]):
	for f2_primer,f2_data in enumerate(df_F2.itertuples()):
		if not get_satisfied_spacing(f2_data.pos - (f3_data.pos + f3_data.length),"F3|3',F2|5'"):
			continue

		if f3_data.score < f2_data.max_before_score:
			continue

		df_F2.iat[f2_primer,-1] = f3_data.score
		for f1c_primer,f1c_data in enumerate(df_F1c.itertuples()):
			if not get_satisfied_spacing(f1c_data.pos - f2_data.pos,"F2|5',F1c|3'"):
				continue
			if f3_data.score + f2_data.score < f1c_data.max_before_score:
				continue
			df_F1c.iat[f1c_primer,-1] = f3_data.score + f2_data.score
			for b1c_primer,b1c_data in enumerate(df_B1c.itertuples()):
				if f3_data.score + f2_data.score + f1c_data.score < b1c_data.max_before_score:
					continue
				if f1c_data.pos+f1c_data.length >= b1c_data.pos:
					continue
				df_B1c.iat[b1c_primer,-1] = f3_data.score + f2_data.score + f1c_data.score
				for b2_primer,b2_data in enumerate(df_B2.itertuples()):
					#for B1c
					if not get_satisfied_spacing(b2_data.pos + b2_data.length - (b1c_data.pos + b1c_data.length),"B1c|3',B2|5'"):
						continue
					#for B2
					if not get_satisfied_spacing(b2_data.pos+b2_data.length - f2_data.pos,"F2|5',B2|5'"):
						continue
					if f3_data.score + f2_data.score + f1c_data.score + b1c_data.score < b2_data.max_before_score:
						continue
					df_B2.iat[b2_primer,-1] = f3_data.score + f2_data.score + f1c_data.score + b1c_data.score
					for b3_primer,b3_data in enumerate(df_B3.itertuples()):
						if not get_satisfied_spacing(b3_data.pos - (b2_data.pos + b2_data.length),"B2|5',B3|3'"):
							continue
						if f3_data.score + f2_data.score + f1c_data.score + b1c_data.score + b2_data.score < b3_data.max_before_score:
							continue
						df_B3.iat[b3_primer,-1] = f3_data.score + f2_data.score + f1c_data.score + b1c_data.score + b2_data.score
						total_score = f3_data.score + f2_data.score + f1c_data.score + b1c_data.score + b2_data.score + b3_data.score
						if total_score > sum(top[-1]["score"]):
							cur_index = len(top) - 1
							while (cur_index > 0 and total_score > sum(top[cur_index]["score"])):
								cur_index-=1
							#print(top)
							top = top[0:cur_index]+[get_top_unit()]+top[cur_index:-1]

							top[cur_index]["set_position"][0] = f3_data.pos
							top[cur_index]["set_position"][1] = f2_data.pos
							top[cur_index]["set_position"][2] = f1c_data.pos
							top[cur_index]["set_position"][3] = b1c_data.pos
							top[cur_index]["set_position"][4] = b2_data.pos
							top[cur_index]["set_position"][5] = b3_data.pos
							top[cur_index]["set_length"][0] = f3_data.length
							top[cur_index]["set_length"][1] = f2_data.length
							top[cur_index]["set_length"][2] = f1c_data.length
							top[cur_index]["set_length"][3] = b1c_data.length
							top[cur_index]["set_length"][4] = b2_data.length
							top[cur_index]["set_length"][5] = b3_data.length

							top[cur_index]["score"][0] = f3_data.score
							top[cur_index]["score"][1] = f2_data.score
							top[cur_index]["score"][2] = f1c_data.score
							top[cur_index]["score"][3] = b1c_data.score
							top[cur_index]["score"][4] = b2_data.score
							top[cur_index]["score"][5] = b3_data.score


#Output primer sets
json.dump({i:top[i] for i in range(len(top))},open(snakemake.output[0],"w"),indent=4)


