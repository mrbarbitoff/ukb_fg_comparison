rule prep_uk_clump:
    input:
        sum_table="../resources/UKB/{uk_sample}.gwas.imputed_v3.both_sexes.tsv.bgz",
        variants=VARIANTS
    output:
        temp("../results/{uk_sample}.annot.tsv")
    shell:
        "paste -d, <(zcat {input.variants}) <(zcat {input.sum_table}) | grep -E '(false|low_confidence)' > {output}"

rule plink_clump_uk:
    input:
        "../results/{uk_sample}.annot.tsv"
    output:
        "../results/plink/{uk_sample}.clumped.ranges",
        "../results/plink/{uk_sample}.clumped"
    params:
        threshold=THRESHOLD,
        snp_field="rsid",
        output_prefix="../results/plink/{uk_sample}"
    conda:
        "../envs/plink.yaml"
    shell:
        "plink --bfile {PLINK_GENOTYPE} --clump {input} "
        "--clump-p1 {params.threshold} "
        "--clump-kb 250 "
        "--clump-snp-field {params.snp_field} "
        "--clump-field pval --allow-extra-chr "
        "--out {params.output_prefix} --clump-range {GLIST};"
        "test -f {output[0]} || touch {output[0]};"
        "test -f {output[1]} || touch {output[1]};"  # empty plink output workaround

rule plink_clump_finn:
    input:
        "../resources/finngen/{finn_sample}_hg19lifted.tsv.gz"
    output:
        "../results/plink/{finn_sample}.clumped.ranges",
        "../results/plink/{finn_sample}.clumped"
    params:
        threshold=THRESHOLD,
        snp_field="rsids",
        output_prefix="../results/plink/{finn_sample}"
    conda:
        "../envs/plink.yaml"
    shell:
        "plink --bfile {PLINK_GENOTYPE} --clump {input} "
        "--clump-p1 {params.threshold} "
        "--clump-kb 250 "
        "--clump-snp-field {params.snp_field} "
        "--clump-field pval --allow-extra-chr "
        "--out {params.output_prefix} --clump-range {GLIST};"
        "test -f {output[0]} || touch {output[0]};"
        "test -f {output[1]} || touch {output[1]};"  # empty plink output workaround

rule plink_clump_meta:
    input:
        "../resources/meta/{meta_sample}1.TBL.gz"
    output:
        "../results/plink/{meta_sample}.clumped.ranges",
        "../results/plink/{meta_sample}.clumped"
    params:
        threshold=THRESHOLD,
        snp_field="MarkerName",
        output_prefix="../results/plink/{meta_sample}"
    conda:
        "../envs/plink.yaml"
    shell:
        "plink --bfile {PLINK_GENOTYPE} --clump {input} "
        "--clump-p1 {params.threshold} "
        "--clump-kb 250 "
        "--clump-snp-field {params.snp_field} "
        "--clump-field P-value --allow-extra-chr "
        "--out {params.output_prefix} --clump-range {GLIST};"
        "test -f {output[0]} || touch {output[0]};"
        "test -f {output[1]} || touch {output[1]};"  # empty plink output workaround

rule range2bed:
    input:
        range="../results/plink/{sample}.clumped.ranges",
        rsid_lookup="../results/rsid_lookup.json",
    output:
        bed=temp("../results/plink/{sample}_signSNP_clump.unmerged.bed")
    params:
        flank=50000
    script:
        "../scripts/range2bed.py"

rule sort_bed:
    input:
        bed="../results/plink/{sample}_signSNP_clump.unmerged.bed",
    output:
        bed=temp("../results/plink/{sample}_signSNP_clump.sorted.bed"),
    shell:
        "sort -k1,1 -k2,2n {input.bed} > {output.bed}"

rule merge_bed:
    input:
        bed="../results/plink/{sample}_signSNP_clump.sorted.bed",
    output:
        bed="../results/plink/{sample}_signSNP_clump.bed",
    shell:
        "bedtools merge -i {input.bed} > {output.bed}"
