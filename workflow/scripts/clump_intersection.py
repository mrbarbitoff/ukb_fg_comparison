import pandas as pd
import subprocess
import os
from collections import defaultdict as dd
from pathlib import Path



corr_table = pd.read_table(snakemake.input.corr_table)

intersect_clump_fu_n = list()
intersect_clump_fm_n = list()
intersect_clump_um_n = list()
intersect_clump_fum_n = list()


for index, row in corr_table.iterrows():
    id = row["fg_phenotype"]

    for dataset in ["finn", "ukbb", "meta"]:

        clump_file_name = f"../results/plink/{dataset}/{id}_signSNP_clump.bed"

        with open(clump_file_name) as clump_file:
            if not os.path.isfile(clump_file_name):
                print(f"{id} file is missing")
                Path(clump_file_name).touch()



    with open(f"../results/plink/{id}_intersection.unmerged.bed", "w") as outfile, \
        open(f"../results/plink/ukbb/{id}_signSNP_clump.bed") as uk_bed, \
        open(f"../results/plink/finn/{id}_signSNP_clump.bed") as finn_bed, \
        open(f"../results/plink/meta/{id}_signSNP_clump.bed") as meta_bed:

        for line in uk_bed:
            outfile.write(f"{line.strip()}\t1\n")

        for line in finn_bed:
            outfile.write(f"{line.strip()}\t0\n")

        for line in meta_bed:
            outfile.write(f"{line.strip()}\t2\n")


    subprocess.call(f"sort -k1,1 -k2,2n ../results/plink/{id}_intersection.unmerged.bed > ../results/plink/{id}_intersection.sorted.bed", shell=True)
    subprocess.call(f"bedtools merge -i ../results/plink/{id}_intersection.sorted.bed -c 5 -o collapse -delim '|' > ../results/plink/{id}_intersection.merged.bed", shell=True)

    intersection_counter = dd(int)

    with open(f"../results/plink/{id}_intersection.merged.bed") as merged:
        for line in merged:
            values = [int(i) for i in line.split()[3].split("|")]
            if len(values) == 3:  # triple intersection
                intersection_counter["fum"] += 1
            elif len(values) == 1:
                continue

            # All double intersections
            if sum(values) == 1:
                intersection_counter["fu"] += 1
            elif sum(values) == 2:
                intersection_counter["fm"] += 1
            else:
                intersection_counter["um"] += 1

    intersect_clump_fu_n.append(intersection_counter['fu'])
    intersect_clump_fm_n.append(intersection_counter['fm'])
    intersect_clump_um_n.append(intersection_counter['um'])
    intersect_clump_fum_n.append(intersection_counter['fum'])


    os.remove(f"../results/plink/{id}_intersection.unmerged.bed")
    os.remove(f"../results/plink/{id}_intersection.sorted.bed")
    os.remove(f"../results/plink/{id}_intersection.merged.bed")


corr_table['intersect_clump_fu_n'] = intersect_clump_fu_n
corr_table['intersect_clump_fm_n'] = intersect_clump_fm_n
corr_table['intersect_clump_um_n'] = intersect_clump_um_n
corr_table['intersect_clump_fum_n'] = intersect_clump_fum_n


corr_table.to_csv(snakemake.output.clump_intersections, sep="\t", index=False)
