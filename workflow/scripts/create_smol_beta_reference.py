import json
import gzip
from collections import defaultdict as dd
from datetime import datetime


finn_pheno = snakemake.input.finn_pheno
uk_pheno = snakemake.input.uk_pheno
SNPs = snakemake.input.SNPs

result = dd(dict)
SNP_list = list()
finn_list = list()
uk_list = list()

with open(finn_pheno) as finn_input_pheno, open(uk_pheno) as uk_input_pheno, open(SNPs) as SNPs_input:
    for line in finn_input_pheno:
        finn_list.append(line.strip())

    for line in uk_input_pheno:
        uk_list.append(line.strip())

    for line in SNPs_input:
        print(line)
        SNP_list.append(line.strip())

SNP_list = set(SNP_list)

print(f"Starting {datetime.now().strftime('%D:%H:%M:%S')}")

for finn_id in finn_list:
    print(f"Processing {finn_id}")
    with gzip.open(f"../resources/finngen/{finn_id}_hg19lifted.tsv.gz", "rt") as finn_data:
        finn_data.readline()
        for finn_line in finn_data:
            chr = finn_line.split()[19].replace("chr", "")
            coord = finn_line.split()[20]
            ref = finn_line.split()[2]
            alt = finn_line.split()[3]
            finn_varid = f'{chr}:{coord}:{ref}:{alt}'
            if finn_varid in SNP_list:
                result[finn_varid][finn_id] = finn_line.split()[8]
else:
    print(f"Finn processing complete {datetime.now().strftime('%D:%H:%M:%S')}")

for uk_id in uk_list:
    print(f"Processing {uk_id}")
    with gzip.open(f"../resources/UKB/{uk_id}.gwas.imputed_v3.both_sexes.tsv.bgz", "rt") as uk_data:
        uk_data.readline()
        for uk_line in uk_data:
            uk_varid = uk_line.split()[0]
            if uk_varid in SNP_list:
                result[uk_varid][uk_id] = uk_line.split()[8]
else:
    print("UK processing complete")

print(f"Finishing {datetime.now().strftime('%D:%H:%M:%S')}")

with open(snakemake.output.beta_reference, 'w') as beta_ref:
    json.dump(result, beta_ref)
