import pandas as pd
import subprocess
import json

# load rsid lookup

with open(snakemake.input.rsid_lookup, 'r') as fp:
    rsid_lookup = json.load(fp)

corr_table = pd.read_table(snakemake.input.corr_table)

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

all_added_columns = [finn_SNP_n, uk_SNP_n, meta_SNP_n, shared_SNP_fu_n, shared_SNP_fm_n, shared_SNP_um_n, shared_SNP_fum_n, uk_clump_n, finn_clump_n, meta_clump_n, intersect_clump_fu_n, intersect_clump_fm_n, intersect_clump_um_n, intersect_clump_fum_n]

for index, row in corr_table.iterrows():
    finn_id = row["Finn_code"]
    uk_id = row["UK_code"]
    meta_id = f'{finn_id}___{uk_id}'


    finn_SNPs = set()
    uk_SNPs = set()
    meta_SNPs = set()

    with open(f"../results/filtered/{finn_id}_filtered.tsv") as finn_snp:
        for line in finn_snp:
            chr, coord, ref, alt = line.split()[0:4]
            finn_varid = f'{chr}:{coord}:{ref}:{alt}'
            finn_SNPs.add(finn_varid)

    with open(f"../results/filtered/{uk_id}_filtered.tsv") as uk_snp:
        for line in uk_snp:
            uk_varid = line.split()[0]
            uk_SNPs.add(uk_varid)

    with open(f"../results/filtered/{meta_id}_filtered.tsv") as meta_snp:
        for line in meta_snp:
            rsid = line.split()[0]
            meta_varid = rsid_lookup[rsid]
            meta_SNPs.add(meta_varid)



    finn_SNP_n.append(len(finn_SNPs))
    uk_SNP_n.append(len(uk_SNPs))
    meta_SNP_n.append(len(meta_SNPs))

    shared_SNP_fu_n.append(len(finn_SNPs & uk_SNPs))
    shared_SNP_fm_n.append(len(finn_SNPs & meta_SNPs))
    shared_SNP_um_n.append(len(meta_SNPs & uk_SNPs))
    shared_SNP_fum_n.append(len(meta_SNPs & uk_SNPs & finn_SNPs))

    # Part about shared clumps

    # 1. Line number of beds
    uk_interval_number = len(open(f"../results/plink/{uk_id}_signSNP_clump.bed").readlines())
    uk_clump_n.append(uk_interval_number)

    finn_interval_number = len(open(f"../results/plink/{finn_id}_signSNP_clump.bed").readlines())
    finn_clump_n.append(finn_interval_number)

    meta_interval_number = len(open(f"../results/plink/{meta_id}_signSNP_clump.bed").readlines())
    meta_clump_n.append(meta_interval_number)

    # 2. Run bedtools intersect
    intersected_fu_filename = f"../results/plink/{uk_id}_{finn_id}.intersect.bed"
    subprocess.call(f"bedtools intersect -a ../results/plink/{uk_id}_signSNP_clump.bed -b ../results/plink/{finn_id}_signSNP_clump.bed | sort -u > {intersected_fu_filename}", shell=True)

    intersected_fm_filename = f"../results/plink/{finn_id}_{meta_id}.intersect.bed"
    subprocess.call(f"bedtools intersect -a ../results/plink/{finn_id}_signSNP_clump.bed -b ../results/plink/{meta_id}_signSNP_clump.bed | sort -u > {intersected_fm_filename}", shell=True)

    intersected_um_filename = f"../results/plink/{uk_id}_{meta_id}.intersect.bed"
    subprocess.call(f"bedtools intersect -a ../results/plink/{uk_id}_signSNP_clump.bed -b ../results/plink/{meta_id}_signSNP_clump.bed | sort -u > {intersected_um_filename}", shell=True)

    intersected_fum_filename = f"../results/plink/{uk_id}_{finn_id}_{meta_id}.intersect.bed"
    subprocess.call(f"bedtools intersect -a ../results/plink/{uk_id}_signSNP_clump.bed -b ../results/plink/{finn_id}_signSNP_clump.bed ../results/plink/{meta_id}_signSNP_clump.bed | sort -u > {intersected_fum_filename}", shell=True)



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


corr_table.to_csv(snakemake.output.updated_table, sep="\t", index=False)
