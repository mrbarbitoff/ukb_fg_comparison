# LDAK requires in summary:
# Predictor A1 A2 Direction Stat n
# rs2260895 A G 1 0.7656 424
# Predictor/SNP could be in chr:bp format
# Example input:
# Predictor A1 A2 Direction Stat n
# rs2260895 A G 1 0.7656 424

import numpy
import gzip

input = snakemake.input[0]
output = snakemake.output[0]

seen_SNPs = set()

with gzip.open(input, 'rt') as infile, open(output, "w") as outfile:
    outfile.write("SNP\tA1\tA2\tDirection\tP\tn\n")  # header
    infile.readline()  # skip header
    for line in infile:
        data = line.strip().split()
        pval = data[6]
        if pval == "NaN":  # remove NaNs
            continue
        direction = str(numpy.sign(float(data[7])))
        n = snakemake.params.n_samples
        # genotype = data[0].split(":")
        chr = data[-2].replace('chr', '')
        coord = data[-1]
        SNP = f'{chr}:{coord}'
        if SNP in seen_SNPs:  # remove repetitive SNPs
            continue
        else:
            seen_SNPs.add(SNP)
        A1 = data[2]
        A2 = data[3]
        if len(A2) > 1 or len(A1) > 1:  # remove indels (may be a bad thing)
            continue
        outfile.write(f"{SNP}\t{A1}\t{A2}\t{direction}\t{pval}\t{n}\n")
