rule generic_filter:
    input:
        sum_table="../resources/meta/{sample}_meta_out.tsv.gz"
    output:
        fg_gw="../results/filtered/finn/{sample}_filtered.tsv",
        uk_gw="../results/filtered/uk/{sample}_filtered.tsv",
        mt_gw="../results/filtered/meta/{sample}_filtered.tsv",
        fg_nom="../results/filtered/nominal_finn/{sample}_filtered.tsv",
        uk_nom="../results/filtered/nominal_uk/{sample}_filtered.tsv",
        fg_repr="../results/filtered/repr_finn/{sample}_filtered.tsv",
        uk_repr="../results/filtered/repr_uk/{sample}_filtered.tsv",
    params:
        threshold=THRESHOLD,
        nominal_threshold=0.05,
        finn_repr_threshold=0.05 / 8207,
        uk_repr_threshold=0.05 / 19514,
    script:
        "../scripts/generic_filter.py"
