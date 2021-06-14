import pandas as pd
from tqdm import tqdm
import numpy as np
import json

params = json.load(open(snakemake.input[1],"r"))

def get_top_unit():
	get_dict = lambda: {k : [-1,-1,-1,-1,-1,-1,-1,-1] for k in ["F3","F2","F1c","B1c","B2","B3","LoopF","LoopB"]}
	return {"set_position":get_dict(),"set_length":get_dict(),"score":get_dict(),"thermo":{"tm":get_dict(),"end":get_dict(),"gc":get_dict()}}
top = [get_top_unit() for _ in range(params["top"])]
#top = {"set":[-1,-1,-1,-1,-1,-1],"score":[-1,-1,-1,-1,-1,-1]}

df_all_primers = pd.read_csv(snakemake.input[0])
df_all_primers["max_before_score"] = np.array([-100]*df_all_primers.shape[0])

sort_regions_by = lambda region,by : df_all_primers[df_all_primers["region"] == region].sort_values(by=by,ascending=False)
#df_F3 = df_all_primers[df_all_primers["region"] == "F3"].sort_values(by="score",ascending=False)
#df_F2 = df_all_primers[df_all_primers["region"] == "F2"].sort_values(by="score",ascending=False)
#df_F1c = df_all_primers[df_all_primers["region"] == "F1c"].sort_values(by="score",ascending=False)
#df_B1c = df_all_primers[df_all_primers["region"] == "B1c"].sort_values(by="score",ascending=False)
#df_B2 = df_all_primers[df_all_primers["region"] == "B2"].sort_values(by="score",ascending=False)
#df_B3 = df_all_primers[df_all_primers["region"] == "B3"].sort_values(by="score",ascending=False)

separated_dfs = {k : df_all_primers[df_all_primers["region"] == k].sort_values(by="score",ascending=False) for k in dict.fromkeys(df_all_primers.region).keys()}

max_before_score = {region:np.array([-100]*df_all_primers[df_all_primers["region"] == region].shape[0]) for region in dict.fromkeys(df_all_primers.region).keys()}

def get_satisfied_spacing(spacing,region):
	return not (spacing > params["spacing"][region]["maximum"] or spacing < params["spacing"][region]["minimum"])

for f3_primer,f3_data in tqdm(enumerate(separated_dfs["F3"].itertuples()),total=separated_dfs["F3"].shape[0],desc="F3"):
	iterate_f2 = tqdm(enumerate(separated_dfs["F2"].itertuples()),total=separated_dfs["F2"].shape[0],desc="F2") if f3_primer == 0 else enumerate(separated_dfs["F2"].itertuples())
	for f2_primer,f2_data in iterate_f2:
		if not get_satisfied_spacing(f2_data.pos - (f3_data.pos + f3_data.length),"F3|3',F2|5'"):
			continue

		if f3_data.score < max_before_score["F2"][f2_primer]:
			continue
		max_before_score["F2"][f2_primer] = f3_score

		for f1c_primer,f1c_data in enumerate(separated_dfs["F1c"].itertuples()):
			if not get_satisfied_spacing(f1c_data.pos - f2_data.pos,"F2|5',F1c|3'"):
				continue
			
			f2_score = f2_data.score
			if f3_data.score + f2_data.score < max_before_score["F1c"][f1c_primer]:
				continue
			max_before_score["F1c"][f1c_primer] = f3_data.score + f2_data.score

			for b1c_primer,b1c_data in enumerate(separated_dfs["B1c"].itertuples()):

				if f1c_data.pos+f1c_data.length >= b1c_data.pos:
					continue
				
				if f3_data.score + f2_data.score + f1c_data.score < max_before_score["B1c"][b1c_primer]:
					continue
				max_before_score["B1c"][b1c_primer] = f3_data.score + f2_data.score + f1c_data.score

				for b2_primer,b2_data in enumerate(separated_dfs["B2"].itertuples()):
					#for B1c
					if not get_satisfied_spacing(b2_data.pos + b2_data.length - (b1c_data.pos + b1c_data.length),"B1c|3',B2|5'"):
						continue
					#for B2
					if not get_satisfied_spacing(b2_data.pos+b2_data.length - f2_data.pos,"F2|5',B2|5'"):
						continue
					
					if f3_data.score + f2_data.score + f1c_data.score + b1c_data.score < max_before_score["B2"][b2_primer]:
						continue
					max_before_score["B2"][b2_primer] = f3_data.score + f2_data.score + f1c_data.score + b1c_data.score

					for b3_primer,b3_data in enumerate(separated_dfs["B3"].itertuples()):
						if not get_satisfied_spacing(b3_data.pos - (b2_data.pos + b2_data.length),"B2|5',B3|3'"):
							continue

						"""
						if f3_data.score + f2_data.score + f1c_data.score + b1c_data.score + b2_data.score < b3_data.max_before_score:
							continue
						df_B3.iat[b3_primer,-1] = f3_data.score + f2_data.score + f1c_data.score + b1c_data.score + b2_data.score
						total_score = f3_data.score + f2_data.score + f1c_data.score + b1c_data.score + b2_data.score + b3_data.score"""
						if f3_data.score + f2_data.score + f1c_data.score + b1c_data.score + b2_data.score < max_before_score["B3"][b3_primer]:
							continue
						max_before_score["B3"][b3_primer] = f3_data.score + f2_data.score + f1c_data.score + b1c_data.score + b2_data.score

						for loopF_primer,loopF_data in enumerate(separated_dfs["LoopF"].itertuples()):
							if not (loopF_data.pos >= f2_data.pos + f2_data.length and loopF_data.pos + loopF_data.length <= f1c_data.pos):
								continue
							if f3_data.score + f2_data.score + f1c_data.score + b1c_data.score + b2_data.score + b3_data.score < max_before_score["LoopF"][loopF_primer]:
								continue
							max_before_score["LoopF"][loopF_primer] = f3_data.score + f2_data.score + f1c_data.score + b1c_data.score + b2_data.score + b3_data.
							for loopB_primer,loopB_data in enumerate(separated_dfs["LoopB"].itertuples()):
								if not (loopB_data.pos >= b1c_data.pos + b1c_data.length and loopB_data.pos + loopB_data.length <= b2_data.pos):
									continue
								if f3_data.score + f2_data.score + f1c_data.score + b1c_data.score + b2_data.score + b3_data.score + loopF_data.score < max_before_score["LoopB"][loopB_primer]:
									continue
								max_before_score["LoopB"][loopB_primer] = f3_data.score + f2_data.score + f1c_data.score + b1c_data.score + b2_data.score + b3_data.score + loopF_data.score

								total_score = f3_data.score + f2_data.score + f1c_data.score + b1c_data.score + b2_data.score + b3_data.score + loopF_data.score + loopB_data.score

								if total_score > sum(top[-1]["score"].values()):
							
									cur_index = len(top) - 1
									while (cur_index > 0 and total_score > sum(top[cur_index]["score"].values())):
										cur_index-=1
									#print(top)
									top = top[0:cur_index]+[get_top_unit()]+top[cur_index:-1]

									top[cur_index]["set_position"]["F3"] = f3_data.pos
									top[cur_index]["set_position"]["F2"] = f2_data.pos
									top[cur_index]["set_position"]["F1c"] = f1c_data.pos
									top[cur_index]["set_position"]["B1c"] = b1c_data.pos
									top[cur_index]["set_position"]["B2"] = b2_data.pos
									top[cur_index]["set_position"]["B3"] = b3_data.pos
									top[cur_index]["set_position"]["LoopF"] = loopF_data.pos
									top[cur_index]["set_position"]["LoopB"] = loopB_data.pos
									top[cur_index]["set_length"]["F3"] = f3_data.length
									top[cur_index]["set_length"]["F2"] = f2_data.length
									top[cur_index]["set_length"]["F1c"] = f1c_data.length
									top[cur_index]["set_length"]["B1c"] = b1c_data.length
									top[cur_index]["set_length"]["B2"] = b2_data.length
									top[cur_index]["set_length"]["B3"] = b3_data.length
									top[cur_index]["set_length"]["LoopF"] = loopF_data.length
									top[cur_index]["set_length"]["LoopB"] = loopB_data.length

									top[cur_index]["score"]["F3"] = f3_data.score
									top[cur_index]["score"]["F2"] = f2_data.score
									top[cur_index]["score"]["F1c"] = f1c_score
									top[cur_index]["score"]["B1c"] = b1c_score
									top[cur_index]["score"]["B2"] = b2_score
									top[cur_index]["score"]["B3"] = b3_data.score
									top[cur_index]["score"]["LoopF"] = loopF_data.score
									top[cur_index]["score"]["LoopB"] = loopB_data.score

									top[cur_index]["thermo"]["tm"]["F3"] = f3_data.tm
									top[cur_index]["thermo"]["tm"]["F2"] = f2_data.tm
									top[cur_index]["thermo"]["tm"]["F1c"] = f1c_data.tm
									top[cur_index]["thermo"]["tm"]["B1c"] = b1c_data.tm
									top[cur_index]["thermo"]["tm"]["B2"] = b2_data.tm
									top[cur_index]["thermo"]["tm"]["B3"] = b3_data.tm
									top[cur_index]["thermo"]["tm"]["LoopF"] = loopF_data.tm
									top[cur_index]["thermo"]["tm"]["LoopB"] = loopB_data.tm

									top[cur_index]["thermo"]["end"]["F3"] = f3_data.end
									top[cur_index]["thermo"]["end"]["F2"] = f2_data.end
									top[cur_index]["thermo"]["end"]["F1c"] = f1c_data.end
									top[cur_index]["thermo"]["end"]["B1c"] = b1c_data.end
									top[cur_index]["thermo"]["end"]["B2"] = b2_data.end
									top[cur_index]["thermo"]["end"]["B3"] = b3_data.end
									top[cur_index]["thermo"]["end"]["LoopF"] = loopF_data.end
									top[cur_index]["thermo"]["end"]["LoopB"] = loopB_data.end

									top[cur_index]["thermo"]["gc"]["F3"] = f3_data.gc
									top[cur_index]["thermo"]["gc"]["F2"] = f2_data.gc
									top[cur_index]["thermo"]["gc"]["F1c"] = f1c_data.gc
									top[cur_index]["thermo"]["gc"]["B1c"] = b1c_data.gc
									top[cur_index]["thermo"]["gc"]["B2"] = b2_data.gc
									top[cur_index]["thermo"]["gc"]["B3"] = b3_data.gc
									top[cur_index]["thermo"]["gc"]["LoopF"] = loopF_data.gc
									top[cur_index]["thermo"]["gc"]["LoopB"] = loopB_data.gc


#Output primer sets
json.dump({i:top[i] for i in range(len(top))},open(snakemake.output[0],"w"),indent=4)


