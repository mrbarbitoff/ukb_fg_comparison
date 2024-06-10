import pandas as pd
from scipy.stats import norm
from collections import defaultdict as dd
import gzip

alpha = snakemake.params.threshold
inverted_threshold = norm.ppf(alpha / 2)

dataset_thresholds = dict()
dataset_thresholds["finn"] = norm.ppf(snakemake.params.finn_repr_threshold / 2)
dataset_thresholds["uk"] = norm.ppf(snakemake.params.uk_repr_threshold / 2)

def repr(beta, se, inverted_threshold):
    first = -abs(beta)/se
    return norm.cdf(first + inverted_threshold) + 1 - norm.cdf(first - inverted_threshold)

datasets = {"finn": 5, "uk": 11}
corr_table = pd.read_table(snakemake.input.corr_table)

outfile = open(snakemake.output.repr_table, "w")
outfile.write("ID\tFinn_SNP_N\tFinn_index_N\tUKB_SNP_N\tUKB_index_N\n")

total_lead_SNPs = dd(lambda: dd(float))

for index, row in corr_table.iterrows():

    id = row["fg_phenotype"]

    outline = f"{id}\t"

    lead_SNPs = {"finn": set(), "ukbb": set()}

# check clumps
    for dataset in ["finn", "ukbb"]:
        clump_file_name = f"../results/plink/{dataset}/{id}.clumped"

        with open(clump_file_name) as clump_file:
            clump_file.readline()
            for line in clump_file:
                try:
                    rsid = line.split()[2]
                    lead_SNPs[dataset].add(rsid)

                except IndexError:
                    break


    for dataset in ["finn", "uk"]:
        index_counter = 0
        prob_sum = 0
        with open(f"../results/filtered/{dataset}/{id}_filtered.tsv") as snp_file:
            for line in snp_file:
                rsid = line.split()[21]

                # revert datasets
                if dataset == "uk":
                    key = "ukbb"
                    repr_dataset = "finn"
                else:
                    key = "finn"
                    repr_dataset = "uk"

                if not rsid in lead_SNPs[key]:
                    continue
                index_counter += 1
                beta = line.split()[datasets[repr_dataset]]

                if beta == "NA":
                    beta = 0
                    se = 1
                else:
                    beta = float(beta)
                    se = float(line.split()[datasets[repr_dataset] + 1])
                inv_thresh = dataset_thresholds[repr_dataset]
                repr_prob = repr(beta, se, inv_thresh)

                # update total lead SNPs
                if total_lead_SNPs[dataset][rsid] < repr_prob:
                    total_lead_SNPs[dataset][rsid] = repr_prob
                    total_lead_SNPs["Total"][rsid] = repr_prob

                prob_sum += repr_prob
        outline += f"{str(prob_sum)}\t{str(index_counter)}\t"
    outline += "\n"
    outfile.write(outline)


total_line = "Total\t"
for key in ["finn", "uk"]:
    key_sum = 0
    key_counter = 0
    for inner_key in total_lead_SNPs[key]:
        key_counter += 1
        key_sum += total_lead_SNPs[key][inner_key]
    total_line += f"{str(key_sum)}\t{str(key_counter)}\t"
total_line += "\n"
outfile.write(total_line)
outfile.close()
