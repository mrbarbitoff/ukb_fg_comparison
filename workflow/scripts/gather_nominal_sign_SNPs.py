import pandas as pd
from collections import defaultdict as dd
import json


# dictionary dataset: [index for beta]
datasets = {"finn": 5, "uk": 11}

# load rsid lookup
with open(snakemake.input.rsid_lookup, 'r') as fp:
    rsid_lookup = json.load(fp)

# create genome-wide SNPs set
genome_wide_SNPs = set()

with open(snakemake.input.filtered_SNP_table) as filtered_SNP_table:
    filtered_SNP_table.readline()
    for line in filtered_SNP_table:
        varid = line.strip().split()[0]
        genome_wide_SNPs.add(varid)

class SNP:
    def __init__(self, pheno_id, beta):
        self.beta = beta
        self.pheno_id = pheno_id

    def __str__(self):
        return f"{self.pheno_id}:{self.beta}"

    def __eq__(self, obj):
        return isinstance(obj, SNP) and obj.pheno_id == self.pheno_id

    def __hash__(self):
        return hash(self.pheno_id)


# convert a set of SNP objects to string representations
def convert_set(inset):
    return set(map(str, inset))

# Create SNP dictionary
SNP_dict = dd(lambda: dd(set))

corr_table = pd.read_table(snakemake.input.corr_table)

for index, row in corr_table.iterrows():

    id = row["fg_phenotype"]

    for dataset in datasets:
        current_file=f"../results/filtered/nominal_{dataset}/{id}_filtered.tsv"
        with open(current_file) as snp_file:
            for line in snp_file:


                try:
                    varid = line.split()[4]
                    beta = line.split()[datasets[dataset]]
                except IndexError:
                    print(current_file)
                current_SNP = SNP(id, beta)
                if varid in genome_wide_SNPs:
                    SNP_dict[varid][f"phenonames_{dataset}"].add(current_SNP)

with open(snakemake.input.SNP_template) as SNP_template, \
    open(snakemake.output.nominal_sign_table, "w") as outfile:

    old_header = "varid\tFinn_MAF\tUKB_MAF\trsid"
    header_suffix = "nominal_phenonames_finn\tnominal_pleio_finn\tnominal_phenonames_uk\tnominal_pleio_uk\n"
    new_header = f'{old_header}\t{header_suffix}'
    outfile.write(new_header)


    for line in SNP_template:
        varid = line.strip().split()[0]

        finn_phenotypes = SNP_dict[varid]["phenonames_finn"]
        uk_phenotypes = SNP_dict[varid]["phenonames_uk"]

        phenonames_finn = ','.join(convert_set(finn_phenotypes))
        pleio_finn = len(finn_phenotypes)
        phenonames_uk = ','.join(convert_set(uk_phenotypes))
        pleio_uk = len(uk_phenotypes)

        outfile.write(f'{line.strip()}\t{phenonames_finn}\t{pleio_finn}\t{phenonames_uk}\t{pleio_uk}\n')
