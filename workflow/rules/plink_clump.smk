pval_col_lookup = {"finn": "FINNGEN_pval", "ukbb": "UKBB_pval", "meta": "all_inv_var_meta_p"}

rule plink_total_clump:
    input:
        sum_table="../resources/meta/{sample}_meta_out.tsv.gz",
    params:
        threshold=THRESHOLD,
        snp_field="rsid",
        pval_field=lambda wildcards: pval_col_lookup[wildcards.dataset],
        output_prefix="../results/plink/{dataset}/{sample}",
    conda:
        "../envs/plink.yaml"
    wildcard_constraints:
        dataset="finn|ukbb|meta"
    output:
        "../results/plink/{dataset}/{sample}.clumped.ranges",
        "../results/plink/{dataset}/{sample}.clumped",
    shell:
        "plink --bfile {PLINK_GENOTYPE} "
        "--clump {input} --clump-kb 250 "
        "--clump-p1 {params.threshold} "
        "--clump-snp-field {params.snp_field} "
        "--clump-field {params.pval_field} "
        "--allow-extra-chr "
        "--out {params.output_prefix} "
        "--clump-range {GLIST};"
        "test -f {output[0]} || touch {output[0]};"
        "test -f {output[1]} || touch {output[1]};"

rule range2bed:
    input:
        range="../results/plink/{dataset}/{sample}.clumped.ranges",
        rsid_lookup="../results/rsid_lookup.json",
    output:
        bed=temp("../results/plink/{dataset}/{sample}_signSNP_clump.unmerged.bed")
    params:
        current_dataset=lambda wc: wc.get("dataset"),
        flank=50000
    script:
        "../scripts/range2bed.py"

rule sort_bed:
    input:
        bed="../results/plink/{dataset}/{sample}_signSNP_clump.unmerged.bed",
    output:
        bed=temp("../results/plink/{dataset}/{sample}_signSNP_clump.sorted.bed"),
    shell:
        "sort -k1,1 -k2,2n {input.bed} > {output.bed}"

rule merge_bed:
    input:
        bed="../results/plink/{dataset}/{sample}_signSNP_clump.sorted.bed",
    output:
        bed="../results/plink/{dataset}/{sample}_signSNP_clump.bed",
    shell:
        "bedtools merge -i {input.bed} -c 4 -o distinct > {output.bed}"

rule total_bed:
    input:
        plink_bed=expand("../results/plink/{dataset}/{sample}_signSNP_clump.bed", sample=samples, dataset=["finn", "ukbb", "meta"]),
    output:
        total_bed="../results/merged.total.bed",
    shell:
        "cat ../results/plink/*/*.bed | sort -k1,1 -k2,2n | "
        "bedtools merge -i stdin -c 4 -o -collapse -delim '|' > {output.total_bed};"
