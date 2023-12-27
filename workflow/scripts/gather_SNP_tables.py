import pandas as pd
from collections import defaultdict as dd
import json

# load rsid lookup

# dictionary dataset: [index for beta]
datasets = {"finn": 5, "uk": 11, "meta": 16}

# create beta_reference
beta_reference = dd(lambda: dd(lambda: dd(str)))


with open(snakemake.input.rsid_lookup, 'r') as fp:
    rsid_lookup = json.load(fp)

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

def intersect_sets(set1, set2):
    set1 = {x.split(":")[0] for x in set1}
    set2 = {x.split(":")[0] for x in set2}

    intersection = set1 & set2

    if set1 == set2 and len(set2) > 0:  # check if sets are the same
        return [f"perfect:{str(len(intersection))}", intersection]
    elif len(intersection) > 0:
        return [f"partial:{str(len(intersection))}", intersection]
    else:
        return [f"none:{str(len(intersection))}", intersection]

# Create SNP dictionary
SNP_dict = dd(lambda: dd(set))

corr_table = pd.read_table(snakemake.input.corr_table)

index_dict = dd(set)

for index, row in corr_table.iterrows():

    id = row["fg_phenotype"]

    for dataset in datasets:
        with open(f"../results/filtered/{dataset}/{id}_filtered.tsv") as snp_file:
            for line in snp_file:
                varid = line.split()[4]
                beta = line.split()[datasets[dataset]]

                # get betas for non-GW associations
                for inner_dataset in datasets:
                    inner_beta = line.split()[datasets[inner_dataset]]
                    beta_reference[varid][id][inner_dataset] = inner_beta

                current_SNP = SNP(id, beta)
                SNP_dict[varid][f"phenonames_{dataset}"].add(current_SNP)

# check clumps
    for dataset in ["finn", "ukbb", "meta"]:
        clump_file_name = f"../results/plink/{dataset}/{id}.clumped"

        with open(clump_file_name) as clump_file:
            clump_file.readline()
            for line in clump_file:
                try:
                    rsid = line.split()[2]
                except IndexError:
                    break
                varid = rsid_lookup[rsid]
                index_dict[varid].add(f"index_{dataset}")
                SNP_dict[varid][f"clump_{dataset}"].add(id)

with open(snakemake.input.SNP_template) as SNP_template, \
    open(snakemake.output.phenonames_table, "w") as outfile, \
    open(snakemake.output.filtered_SNP_table, "w") as filtered_SNP_table:

    old_header = "varid\tFinn_MAF\tUKB_MAF\trsid"
    header_suffix = "beta_finn\tbeta_UK\tLD_score\tphenonames_finn\tpleio_finn\tclump_finn\tphenonames_uk\tpleio_uk\tclump_uk\tphenonames_meta\tpleio_meta\tclump_meta\tinter_fu\trepr_fu\tinter_fm\trepr_fm\tinter_um\trepr_um\tis_index\n"
    new_header = f'{old_header}\t{header_suffix}'

    outfile.write(new_header)  # skip header and write outfile header
    filtered_SNP_table.write(new_header)


    for line in SNP_template:
        varid = line.split()[0]

        # check if index in a clump
        is_index = ":".join(index_dict[varid])

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

        # intersections
        intersection_finn_uk = intersect_sets(convert_set(finn_phenotypes), convert_set(uk_phenotypes))
        intersection_finn_meta = intersect_sets(convert_set(finn_phenotypes), convert_set(meta_phenotypes))
        intersection_uk_meta = intersect_sets(convert_set(uk_phenotypes), convert_set(meta_phenotypes))


        if pleio_uk == 0 and pleio_finn == 0 and pleio_meta > 0:

            meta_ps = set([x.split(":")[0] for x in phenonames_meta.split(',')])


            betas_finn = {float(beta_reference[varid][x]["finn"]) for x in meta_ps}
            betas_uk = {float(beta_reference[varid][x]["uk"]) for x in meta_ps}

            beta_finn = max(betas_finn, key=abs)
            beta_UK = max(betas_uk, key=abs)
        else:
            # Workaround for betas if no GW-significant associations found in one BB
            beta_UK = find_max_effect_size(phenonames_uk)
            if beta_UK == "0" and pleio_finn > 0:

                finn_ps = set([x.split(":")[0] for x in phenonames_finn.split(',')])
                betas_uk = {float(beta_reference[varid][x]["uk"]) if not beta_reference[varid][x]["uk"] == "NA" else 0 for x in finn_ps}
                beta_UK = max(betas_uk, key=abs)

            beta_finn = find_max_effect_size(phenonames_finn)
            if beta_finn == "0" and pleio_uk > 0:

                uk_ps = set([x.split(":")[0] for x in phenonames_uk.split(',')])
                betas_finn = {float(beta_reference[varid][x]["finn"]) if not beta_reference[varid][x]["finn"] == "NA" else 0 for x in uk_ps}
                beta_finn = max(betas_finn, key=abs)

        # clumps
        clump_finn = ','.join(SNP_dict[varid]["clump_finn"])
        clump_uk = ','.join(SNP_dict[varid]["clump_uk"])
        clump_meta = ','.join(SNP_dict[varid]["clump_meta"])

        # LD scores
        LD_score = ld_scores_dict[varid]


        # are sets perfectly intersected?
        repr_fu = intersection_finn_uk[0]
        repr_fm = intersection_finn_meta[0]
        repr_um = intersection_uk_meta[0]

        # string with intersecting phenonames
        ifu_str =  ",".join(intersection_finn_uk[1])
        ifm_str = ",".join(intersection_finn_meta[1])
        ium_str = ",".join(intersection_uk_meta[1])


        new_data_line = f'{line.strip()}\t{beta_finn}\t{beta_UK}\t{LD_score}\t{phenonames_finn}\t{pleio_finn}\t{clump_finn}\t{phenonames_uk}\t{pleio_uk}\t{clump_uk}\t{phenonames_meta}\t{pleio_meta}\t{clump_meta}\t{ifu_str}\t{repr_fu}\t{ifm_str}\t{repr_fm}\t{ium_str}\t{repr_um}\t{is_index}\n'
        if pleio_uk > 0 or pleio_finn > 0 or pleio_meta > 0:
            filtered_SNP_table.write(new_data_line)
        outfile.write(new_data_line)
