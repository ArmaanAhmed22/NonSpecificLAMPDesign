rule compile:
	output:
		"test.class"
	conda:
		"env.yaml"
	script:
		 "javac test.java"

rule test:
	output:
		"text.txt"
	conda:
		"env.yaml"
	script:
		"java test"