import gzip


finn_repr_threshold = snakemake.params.finn_repr_threshold
uk_repr_threshold = snakemake.params.uk_repr_threshold

with gzip.open(snakemake.input.sum_table, 'rt') as infile, \
    open(snakemake.output.fg_repr, "w") as fg_repr_file, \
    open(snakemake.output.uk_repr, "w") as uk_repr_file:

    infile.readline()  # skip header
    for line in infile:
        fields = line.split()
        FG_pval = fields[7]
        UKBB_pval = fields[13]
        meta_pval = fields[18]

        if FG_pval != "NA" and float(FG_pval) < finn_repr_threshold:
            fg_repr_file.write(line)
        if UKBB_pval != "NA" and float(UKBB_pval) < uk_repr_threshold:
            uk_repr_file.write(line)
