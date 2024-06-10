rule generic_filter:
    input:
        sum_table="../resources/meta/{sample}_meta_out.tsv.gz"
    output:
        fg_gw="../results/filtered/finn/{sample}_filtered.tsv",
        uk_gw="../results/filtered/uk/{sample}_filtered.tsv",
        mt_gw="../results/filtered/meta/{sample}_filtered.tsv",
        fg_nom="../results/filtered/nominal_finn/{sample}_filtered.tsv",
        uk_nom="../results/filtered/nominal_uk/{sample}_filtered.tsv",
    params:
        threshold=THRESHOLD,
        nominal_threshold=0.05,
    script:
        "../scripts/generic_filter.py"

# requires numbers of index variants for FG and UKB first
rule repr_filter:
    input:
        sum_table="../resources/meta/{sample}_meta_out.tsv.gz",
        filtered_SNP_table="../results/filtered_SNP_table.tsv",
    output:
        fg_repr="../results/filtered/repr_finn/{sample}_filtered.tsv",
        uk_repr="../results/filtered/repr_uk/{sample}_filtered.tsv",
    params:
        finn_repr_threshold=snakemake.ukbb_index_counter,
        uk_repr_threshold=snakemake.finn_index_counter,
    script:
        "../scripts/repr_filter.py"
