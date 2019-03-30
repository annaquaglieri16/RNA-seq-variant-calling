import sys
from collections import defaultdict

samples=defaultdict(list)

with open(sys.argv[1], 'r') as infile:
	for line in infile:
		line = line.strip().split()
		samples[line[0]].append(line[1])

combine_sample = set()
with open(sys.argv[2], 'r') as infile:
	for line in infile:
		combine_sample.add(line.strip())

for sample in samples:
	if sample in combine_sample:
		with open(sample + "_combined_1.fastq", 'w') as outfile:
			for r1 in samples[sample]:
				with open(r1+"_1.fastq", 'r') as infile:
					for line in infile:
						outfile.write(line)
		with open(sample + "_combined_2.fastq", 'w') as outfile:
			for r1 in samples[sample]:
				with open(r1+"_2.fastq", 'r') as infile:
					for line in infile:
						outfile.write(line)
