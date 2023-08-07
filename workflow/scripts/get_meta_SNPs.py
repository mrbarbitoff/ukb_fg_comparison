# Get all SNP rsids available in meta analysis

import gzip

SNP_set = set()

for file in snakemake.input.meta_files:
    with gzip.open(file) as infile:
        infile.readline()  # skip header
        for line in infile:
            rsid = line.split()[0]
            SNP_set.add(rsid)

with open(snakemake.output.meta_rsids, "w") as outfile:
    for item in SNP_set:
        print(item, file=outfile)
