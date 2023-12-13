import gzip

gw_threshold = snakemake.params.threshold
nominal_threshold = snakemake.params.nominal_threshold

with gzip.open(snakemake.input.sum_table, 'rt') as infile, \
    open(snakemake.output.fg_gw, "w") as fg_gw_file, \
    open(snakemake.output.uk_gw, "w") as uk_gw_file, \
    open(snakemake.output.mt_gw, "w") as mt_gw_file, \
    open(snakemake.output.fg_nom, "w") as fg_nom_file, \
    open(snakemake.output.uk_nom, "w") as uk_nom_file:

    infile.readline()  # skip header
    for line in infile:
        fields = line.split()
        SNP = fields[4]
        FG_beta = fields[5]
        FG_pval = fields[7]
        FG_af_alt = fields[8]
        UKBB_beta = fields[11]
        UKBB_pval = fields[13]
        UKBB_af_alt = fields[14]
        meta_beta = fields[16]
        meta_pval = fields[18]
        rsid = fields[21]
        if FG_pval != "NA" and float(FG_pval) < gw_threshold:
            fg_gw_file.write(line.strip())
        if UKBB_pval != "NA" and float(UKBB_pval) < gw_threshold:
            uk_gw_file.write(line.strip())
        if meta_pval != "NA" and float(meta_pval) < gw_threshold:
            mt_gw_file.write(line.strip())
        if FG_pval != "NA" and float(FG_pval) < nominal_threshold:
            fg_nom_file.write(line.strip())
        if UKBB_pval != "NA" and float(UKBB_pval) < nominal_threshold:
            uk_nom_file.write(line.strip())
