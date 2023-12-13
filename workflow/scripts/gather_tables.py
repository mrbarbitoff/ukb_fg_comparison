import pandas as pd
import subprocess
import json

corr_table = pd.read_table(snakemake.input.corr_table)

datasets = ["finn", "uk", "meta"]

finn_SNP_n = list()
uk_SNP_n = list()
meta_SNP_n = list()

shared_SNP_fu_n = list()
shared_SNP_fm_n = list()
shared_SNP_um_n = list()
shared_SNP_fum_n = list()

uk_clump_n = list()
finn_clump_n = list()
meta_clump_n = list()

intersect_clump_fu_n = list()
intersect_clump_fm_n = list()
intersect_clump_um_n = list()
intersect_clump_fum_n = list()

subprocess.call("mkdir ../results/plink/intersections/", shell=True)

for index, row in corr_table.iterrows():
    id = row["fg_phenotype"]

    SNP_dict = {dataset: set() for dataset in datasets}

    for dataset in datasets:
        with open(f"../results/filtered/{dataset}/{id}_filtered.tsv") as filtered_file:
            for line in filtered_file:
                varid = line.split()[4]
                SNP_dict[dataset].add(varid)

    finn_SNP_n.append(len(SNP_dict["finn"]))
    uk_SNP_n.append(len(SNP_dict["uk"]))
    meta_SNP_n.append(len(SNP_dict["meta"]))

    shared_SNP_fu_n.append(len(SNP_dict["finn"] & SNP_dict["uk"]))
    shared_SNP_fm_n.append(len(SNP_dict["finn"] & SNP_dict["meta"]))
    shared_SNP_um_n.append(len(SNP_dict["meta"] & SNP_dict["uk"]))
    shared_SNP_fum_n.append(len(SNP_dict["meta"] & SNP_dict["uk"] & SNP_dict["finn"]))

    # Part about shared clumps

    # 1. Line number of beds
    uk_clump_file = f"../results/plink/ukbb/{id}_signSNP_clump.bed"
    uk_interval_number = len(open(uk_clump_file).readlines())
    uk_clump_n.append(uk_interval_number)


    # finn_clump_file = f"../results/plink/finn/{id}_signSNP_clump.bed"
    finn_clump_file = f"../results/plink/finn/{id}_signSNP_clump.bed"
    finn_interval_number = len(open(finn_clump_file).readlines())
    finn_clump_n.append(finn_interval_number)

    meta_clump_file = f"../results/plink/meta/{id}_signSNP_clump.bed"
    meta_interval_number = len(open(meta_clump_file).readlines())
    meta_clump_n.append(meta_interval_number)

    # 2. Run bedtools intersect

    intersected_fu_filename = f"../results/plink/intersections/uk_finn_{id}.intersect.bed"
    subprocess.call(f"bedtools intersect -a {uk_clump_file} -b {finn_clump_file} | sort -u > {intersected_fu_filename}", shell=True)

    intersected_fm_filename = f"../results/plink/intersections/finn_meta_{id}.intersect.bed"
    subprocess.call(f"bedtools intersect -a {finn_clump_file} -b {meta_clump_file} | sort -u > {intersected_fm_filename}", shell=True)

    intersected_um_filename = f"../results/plink/intersections/uk_meta_{id}.intersect.bed"
    subprocess.call(f"bedtools intersect -a {uk_clump_file} -b {meta_clump_file} | sort -u > {intersected_um_filename}", shell=True)

    intersected_fum_filename = f"../results/plink/intersections/uk_finn_meta_{id}.intersect.bed"
    subprocess.call(f"bedtools intersect -a {uk_clump_file} -b {finn_clump_file} {meta_clump_file} | sort -u > {intersected_fum_filename}", shell=True)



    # 3. Line number of intersect
    intersect_clump_fu_n.append(len(open(intersected_fu_filename).readlines()))
    intersect_clump_fm_n.append(len(open(intersected_fm_filename).readlines()))
    intersect_clump_um_n.append(len(open(intersected_um_filename).readlines()))
    intersect_clump_fum_n.append(len(open(intersected_fum_filename).readlines()))

# add new columns to df

corr_table['finn_SNP_n'] = finn_SNP_n
corr_table['uk_SNP_n'] = uk_SNP_n
corr_table['meta_SNP_n'] = meta_SNP_n

corr_table['shared_SNP_fu_n'] = shared_SNP_fu_n
corr_table['shared_SNP_fm_n'] = shared_SNP_fm_n
corr_table['shared_SNP_um_n'] = shared_SNP_um_n
corr_table['shared_SNP_fum_n'] = shared_SNP_fum_n


corr_table['uk_clump_n'] = uk_clump_n
corr_table['finn_clump_n'] = finn_clump_n
corr_table['meta_clump_n'] = meta_clump_n


corr_table['intersect_clump_fu_n'] = intersect_clump_fu_n
corr_table['intersect_clump_fm_n'] = intersect_clump_fm_n
corr_table['intersect_clump_um_n'] = intersect_clump_um_n
corr_table['intersect_clump_fum_n'] = intersect_clump_fum_n

subprocess.call("rm -r ../results/plink/intersections/", shell=True)

corr_table.to_csv(snakemake.output.updated_table, sep="\t", index=False)
