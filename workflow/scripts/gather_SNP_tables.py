import pandas as pd
from collections import defaultdict as dd
import json
import gzip

# load rsid lookup

with open(snakemake.input.rsid_lookup, 'r') as fp:
    rsid_lookup = json.load(fp)

with open(snakemake.input.beta_reference, 'r') as beta_ref:
    beta_reference = json.load(beta_ref)


# load second_priority
second_priority_set = set()
with open(snakemake.input.second_priority) as sp:
    for line in sp:
        finn_code, uk_code = line.split()
        meta_code = f"{finn_code}___{uk_code}"
        second_priority_set.update([finn_code, uk_code, meta_code])

# LD score dictionary
ld_scores_dict = dd(lambda: "NA")
with open(snakemake.input.ld_scores) as ld_scores:
    for line in ld_scores:
        rsid = line.split()[1]
        LD = line.split()[5]
        try:
            varid = rsid_lookup[rsid]
        except KeyError:
            continue
        ld_scores_dict[varid] = LD

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


# find maximal effect size among several phenotypes
def find_max_effect_size(input_string):
    if not input_string:
        return "0"
    pairs = input_string.split(',')
    betas = [float(x.split(":")[1]) for x in pairs]
    return max(betas, key=abs)


# convert a set of SNP objects to string representations
def convert_set(inset):
    return set(map(str, inset))

def intersect_sets_with_lookup(set1, set2, lookup):
    # lookup: [item in set1] -> list of items in set2
    set1 = {x.split(":")[0] for x in set1}
    set2 = {x.split(":")[0] for x in set2}

    translated_set1 = set()

    for i in set1:
        for j in lookup[i]:
            translated_set1.add(j)
    intersection = set2 & translated_set1

    if translated_set1 == set2 and len(set2) > 0:  # check if sets are the same
        return [f"perfect:{str(len(intersection))}", intersection]
    elif len(intersection) > 0:
        return [f"partial:{str(len(intersection))}", intersection]
    else:
        return [f"none:{str(len(intersection))}", intersection]

# Create SNP dictionary
SNP_dict = dd(lambda: dd(set))

corr_table = pd.read_table(snakemake.input.corr_table)

# create lookups for intersections
fu_lookup = dd(list)
fm_lookup = dd(list)
um_lookup = dd(list)


index_dict = dd(lambda: "non_index")

for index, row in corr_table.iterrows():

    finn_id = row["Finn_code"]
    uk_id = row["UK_code"]

    # only work with second priority
    if not finn_id in second_priority_set or not uk_id in second_priority_set:
        continue

    meta_id = f'{finn_id}___{uk_id}'

    # fill lookups
    fu_lookup[finn_id].append(uk_id)
    fm_lookup[finn_id].append(meta_id)
    um_lookup[uk_id].append(meta_id)

    if not all(f"../results/filtered/{x}_filtered.tsv" in snakemake.input.filtered_tables for x in (finn_id, uk_id, meta_id)):
        continue

    if not all(f"../results/plink/{x}.clumped" in snakemake.input.plink_clumped for x in (finn_id, uk_id, meta_id)):
        continue

    with open(f"../results/filtered/{finn_id}_filtered.tsv") as finn_snp:
        for line in finn_snp:
            chr, coord, ref, alt = line.split()[0:4]
            finn_varid = f'{chr}:{coord}:{ref}:{alt}'
            beta = line.split()[5]
            finn_SNP = SNP(finn_id, beta)
            SNP_dict[finn_varid]["phenonames_finn"].add(finn_SNP)

    with open(f"../results/filtered/{uk_id}_filtered.tsv") as uk_snp:
        for line in uk_snp:
            uk_varid = line.split()[0]
            beta = line.split()[8]
            uk_SNP = SNP(uk_id, beta)
            SNP_dict[uk_varid]["phenonames_uk"].add(uk_SNP)

    with open(f"../results/filtered/{meta_id}_filtered.tsv") as meta_snp:
        for line in meta_snp:
            rsid = line.split()[0]
            meta_varid = rsid_lookup[rsid]
            direction = line.split()[10]
            meta_SNP = SNP(meta_id, direction)
            SNP_dict[meta_varid]["phenonames_meta"].add(meta_SNP)


# check clumps
    for current_id, current_set_name in {finn_id: "clump_finn", uk_id: "clump_uk", meta_id: "clump_meta"}.items():
        with open(f"../results/plink/{current_id}.clumped") as clump_file:
            clump_file.readline()
            for line in clump_file:
                try:
                    rsid = line.split()[2]
                    varid = rsid_lookup[rsid]

                    index_dict[varid] = "index"
                except IndexError:
                    break
                varid = rsid_lookup[rsid]
                SNP_dict[varid][current_set_name].add(current_id)

with open(snakemake.input.SNP_template) as SNP_template, \
    open(snakemake.output.phenonames_table, "w") as outfile, \
    open(snakemake.output.filtered_SNP_table, "w") as filtered_SNP_table, \
    open("meta_specific_SNP.txt", "w") as meta_specific_SNP:

    orphan_finn_set = set()
    orphan_uk_set = set()

    header_suffix = "beta_finn\tbeta_UK\tLD_score\tphenonames_finn\tpleio_finn\tclump_finn\tphenonames_uk\tpleio_uk\tclump_uk\tphenonames_meta\tpleio_meta\tclump_meta\tinter_fu\trepr_fu\tinter_fm\trepr_fm\tinter_um\trepr_um\tis_index\n"
    new_header = f'{SNP_template.readline().strip()}\t{header_suffix}'

    outfile.write(new_header)  # skip header and write outfile header
    filtered_SNP_table.write(new_header)


    for line in SNP_template:
        varid = line.strip().split()[0]

        # check if index in a clump
        is_index = index_dict[varid]

        # lists of phenotypes from SNP dictionary
        finn_phenotypes = SNP_dict[varid]["phenonames_finn"]
        uk_phenotypes = SNP_dict[varid]["phenonames_uk"]
        meta_phenotypes = SNP_dict[varid]["phenonames_meta"]

        phenonames_finn = ','.join(convert_set(finn_phenotypes))
        pleio_finn = len(finn_phenotypes)
        phenonames_uk = ','.join(convert_set(uk_phenotypes))
        pleio_uk = len(uk_phenotypes)
        phenonames_meta = ','.join(convert_set(meta_phenotypes))
        pleio_meta = len(meta_phenotypes)

        # effect sizes
        if pleio_uk == 0 and pleio_finn == 0 and pleio_meta > 0:

            pairs = phenonames_meta.split(',')
            meta_phenotypes_uk = set([x.split(":")[0].split("___")[1] for x in pairs])

            meta_phenotypes_finn = set([x.split(":")[0].split("___")[0] for x in pairs])

            try:
                betas_finn = {float(beta_reference[varid][x]) for x in meta_phenotypes_finn}
            except KeyError:
                beta_finn = {0, }
            betas_uk = {float(beta_reference[varid][x]) for x in meta_phenotypes_uk}

            beta_finn = max(betas_finn, key=abs)
            beta_UK = max(betas_uk, key=abs)

            orphan_finn_set.update(meta_phenotypes_finn)
            orphan_uk_set.update(meta_phenotypes_uk)

        else:
            beta_UK = find_max_effect_size(phenonames_uk)
            beta_finn = find_max_effect_size(phenonames_finn)

        # clumps
        clump_finn = ','.join(SNP_dict[varid]["clump_finn"])
        clump_uk = ','.join(SNP_dict[varid]["clump_uk"])
        clump_meta = ','.join(SNP_dict[varid]["clump_meta"])

        # LD scores
        LD_score = ld_scores_dict[varid]


        # intersections
        intersection_finn_uk = intersect_sets_with_lookup(convert_set(finn_phenotypes), convert_set(uk_phenotypes), fu_lookup)
        intersection_finn_meta = intersect_sets_with_lookup(convert_set(finn_phenotypes), convert_set(meta_phenotypes), fm_lookup)
        intersection_uk_meta = intersect_sets_with_lookup(convert_set(uk_phenotypes), convert_set(meta_phenotypes), um_lookup)

        # are sets perfectly intersected?
        repr_fu = intersection_finn_uk[0]
        repr_fm = intersection_finn_meta[0]
        repr_um = intersection_uk_meta[0]

        # string with intersecting phenonames
        ifu_str =  ",".join(intersection_finn_uk[1])
        ifm_str = ",".join(intersection_finn_meta[1])
        ium_str = ",".join(intersection_uk_meta[1])

        # check replace empty entries with NA
        added_text_columns = [beta_UK, beta_finn, LD_score, \
        phenonames_finn, phenonames_uk, phenonames_meta, \
        clump_uk, clump_finn, clump_meta, \
        intersection_finn_uk, intersection_finn_meta, intersection_uk_meta, \
        repr_fu, repr_fm, repr_um, \
        ifu_str, ifm_str, ium_str]

        for i in range(len(added_text_columns)):
            if not added_text_columns[i]:
                added_text_columns[i] = "NA"

        if pleio_uk > 0 or pleio_finn > 0 or pleio_meta > 0:
            filtered_SNP_table.write(f'{line.strip()}\t{beta_finn}\t{beta_UK}\t{LD_score}\t{phenonames_finn}\t{pleio_finn}\t{clump_finn}\t{phenonames_uk}\t{pleio_uk}\t{clump_uk}\t{phenonames_meta}\t{pleio_meta}\t{clump_meta}\t{ifu_str}\t{repr_fu}\t{ifm_str}\t{repr_fm}\t{ium_str}\t{repr_um}\t{is_index}\n')
        outfile.write(f'{line.strip()}\t{beta_finn}\t{beta_UK}\t{LD_score}\t{phenonames_finn}\t{pleio_finn}\t{clump_finn}\t{phenonames_uk}\t{pleio_uk}\t{clump_uk}\t{phenonames_meta}\t{pleio_meta}\t{clump_meta}\t{ifu_str}\t{repr_fu}\t{ifm_str}\t{repr_fm}\t{ium_str}\t{repr_um}\t{is_index}\n')

    orphan_pheno_finn = open("orphan_pheno_finn.txt", "a")
    orphan_pheno_uk = open("orphan_pheno_uk.txt", "a")

    for opf in orphan_finn_set:
        print(f"{opf}", file=orphan_pheno_finn)
    for opu in orphan_uk_set:
        print(f"{opu}", file=orphan_pheno_uk)
