with open(snakemake.input[0]) as f:
	r = f.read()
	with open(snakemake.output[0],"w") as f2:
		f2.write(r)
