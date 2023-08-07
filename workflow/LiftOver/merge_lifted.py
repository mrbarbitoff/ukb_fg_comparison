import pandas as pd
import gzip

table = pd.read_table(snakemake.input.table, compression="gzip")

lifted_bed = pd.read_table(snakemake.input.lifted_bed,
    names=["new_chr", "new_coord", "coord_2", "ID"],
    dtype='str')

table['ID'] = table['#chrom'].astype('str') + ":" + table['pos'].astype('str')

hg19lifted = table.set_index('ID').join(lifted_bed.set_index('ID'))

hg19lifted.pop('coord_2')  # remove dummy coordinate

hg19lifted.to_csv(snakemake.output.hg19lifted, sep="\t",
        index=False,
        compression="gzip",
        na_rep="NA")
