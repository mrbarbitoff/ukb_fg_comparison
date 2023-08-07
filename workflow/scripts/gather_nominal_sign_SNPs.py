import pandas as pd
from collections import defaultdict as dd
import json
import gzip
import tracemalloc

tracemalloc.start()


# load rsid lookup

with open(snakemake.input.rsid_lookup, 'r') as fp:
    rsid_lookup = json.load(fp)


second_priority_set = set()
with open(snakemake.input.second_priority) as sp:
    for line in sp:
        finn_code, uk_code = line.split()
        meta_code = f"{finn_code}___{uk_code}"
        second_priority_set.update([finn_code, uk_code, meta_code])

meta_rsids_set = set()
with open(snakemake.input.meta_rsids) as meta_rsids:
    for line in meta_rsids:
        meta_rsids_set.add(line.strip().replace("b", "").replace("'", ""))

genome_wide_SNPs = set()
with open(snakemake.input.filtered_SNP_table) as filtered_SNP_table:
    filtered_SNP_table.readline()
    for line in filtered_SNP_table:
        varid = line.strip().split()[0]
        genome_wide_SNPs.add(varid)

print("Lookups loaded")

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

# create lookups for intersections
fu_lookup = dd(list)


for index, row in corr_table.iterrows():

    finn_id = row["Finn_code"]
    uk_id = row["UK_code"]
    print(f"Processing {finn_id} and {uk_id}")

    # trace memory usage in MB
    current, peak = tracemalloc.get_traced_memory()
    print(f"Memory usage: {current / (1024*1024)} and {peak / (1024*1024)}")
    # only work with second priority
    if not finn_id in second_priority_set or not uk_id in second_priority_set:
        continue

    # fill lookups
    fu_lookup[finn_id].append(uk_id)


    if not f"../results/filtered/nominal/{finn_id}_filtered.tsv" in snakemake.input.nominal_filtered_tables:
        print(f"Missing {finn_id}")
        continue
    elif not f"../results/filtered/nominal/{uk_id}_filtered.tsv" in snakemake.input.nominal_filtered_tables:
        print(f"Missing {uk_id}")
        continue

    with open(f"../results/filtered/nominal/{finn_id}_filtered.tsv") as finn_snp:
        for line in finn_snp:
            chr, coord, ref, alt = line.split()[0:4]
            finn_varid = f'{chr}:{coord}:{ref}:{alt}'
            beta = line.split()[5]
            finn_SNP = SNP(finn_id, beta)
            if finn_varid in genome_wide_SNPs:
                SNP_dict[finn_varid]["phenonames_finn"].add(finn_SNP)

    with open(f"../results/filtered/nominal/{uk_id}_filtered.tsv") as uk_snp:
        for line in uk_snp:
            uk_varid = line.split()[0]
            beta = line.split()[8]
            uk_SNP = SNP(uk_id, beta)
            if uk_varid in genome_wide_SNPs:
                SNP_dict[uk_varid]["phenonames_uk"].add(uk_SNP)


with open(snakemake.input.SNP_template) as SNP_template, \
    open(snakemake.output.nominal_sign_table, "w") as outfile:

    header_suffix = "meta_represented\tnominal_phenonames_finn\tnominal_pleio_finn\tnominal_phenonames_uk\tnominal_pleio_uk\n"
    new_header = f'{SNP_template.readline().strip()}\t{header_suffix}'


    for line in SNP_template:
        varid = line.strip().split()[0]
        finn_rsids = line.strip().split()[2].split(',')

        if any([x in meta_rsids_set for x in finn_rsids]):
            meta_represented = "yes"
        else:
            meta_represented = "no"

        finn_phenotypes = SNP_dict[varid]["phenonames_finn"]
        uk_phenotypes = SNP_dict[varid]["phenonames_uk"]


        phenonames_finn = ','.join(convert_set(finn_phenotypes))
        pleio_finn = len(finn_phenotypes)
        phenonames_uk = ','.join(convert_set(uk_phenotypes))
        pleio_uk = len(uk_phenotypes)

        outfile.write(f'{line.strip()}\t{meta_represented}\t{phenonames_finn}\t{pleio_finn}\t{phenonames_uk}\t{pleio_uk}\n')

tracemalloc.stop()
