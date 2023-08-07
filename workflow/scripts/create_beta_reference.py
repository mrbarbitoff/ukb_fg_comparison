import pandas as pd
from collections import defaultdict as dd
import json
import gzip
import h5py


# create index lookup for every SNP varid
varid_index_dict = dict()

with open(snakemake.input.SNP_template) as SNP_template:
    SNP_template.readline()
    for index, line in enumerate(SNP_template):
        varid = line.strip().split()[0]
        varid_index_dict[varid] = index

SNP_N = len(varid_index_dict)

orphanSNPs = open("orphan_SNP.txt", "a")

# open hdf5 file for writing and create groups
with h5py.File(snakemake.output.beta_reference, 'w') as beta_ref_h5:

    finn_group = beta_ref_h5.create_group("Finngen")
    ukb_group = beta_ref_h5.create_group("UKB")


    seen = set()
    # read phenotable and create datasets
    corr_table = pd.read_table(snakemake.input.corr_table)

    for index, row in corr_table.iterrows():

        finn_id = row["Finn_code"]
        uk_id = row["UK_code"]

        if finn_id in seen:
            pass
        else:
            finn_dataset = finn_group.create_dataset(f"{finn_id}", (SNP_N,), fillvalue=0)
            seen.add(finn_id)
            try:
                with gzip.open(f"../resources/finngen/{finn_id}_hg19lifted.tsv.gz", "rt") as finn_data:
                    finn_data.readline()
                    for finn_line in finn_data:
                        chr = finn_line.split()[19].replace("chr", "")
                        coord = finn_line.split()[20]
                        ref = finn_line.split()[2]
                        alt = finn_line.split()[3]
                        finn_varid = f'{chr}:{coord}:{ref}:{alt}'
                        try:
                            index = varid_index_dict[finn_varid]
                            beta_finn = float(finn_line.split()[8])
                            finn_dataset[index] = beta_finn
                        except KeyError:
                            print(f"{finn_varid}", file=orphanSNPs)
            except FileNotFoundError:
                pass

        if uk_id in seen:
            pass
        else:
            ukb_dataset = ukb_group.create_dataset(f"{uk_id}", (SNP_N,), fillvalue=0)
            seen.add(uk_id)
            try:
                with gzip.open(f"../resources/UKB/{uk_id}.gwas.imputed_v3.both_sexes.tsv.bgz", "rt") as uk_data:
                    uk_data.readline()
                    for uk_line in uk_data:
                        uk_varid = uk_line.split()[0]
                        index = varid_index_dict[uk_varid]
                        beta_UK = float(uk_line.split()[7])
                        ukb_dataset[index] = beta_UK
            except FileNotFoundError:
                pass
