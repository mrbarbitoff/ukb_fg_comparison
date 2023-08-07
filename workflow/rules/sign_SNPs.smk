rule filter_significant_finn:
    input:
        "../resources/finngen/{finn_sample}_hg19lifted.tsv.gz"
    output:
        "../results/filtered/{finn_sample}_filtered.tsv"
    params:
        threshold=THRESHOLD
    shell:
        "zcat {input} | cut -f 3,4,7,9,20,21 | sed 's/chr//g' | awk '{{if($3 <= {params.threshold}) {{print $5, $6, $1, $2, $3, $4}}}}' > {output}"

rule filter_significant_uk:
    input:
        "../resources/UKB/{uk_sample}.gwas.imputed_v3.both_sexes.tsv.bgz"
    output:
        "../results/filtered/{uk_sample}_filtered.tsv"
    params:
        threshold=THRESHOLD,
    shell:
        "zcat {input} | grep \"false\" | awk '{{if($12 <= {params.threshold}) {{print}}}}' > {output}"

rule nominal_filter_significant_finn:
    input:
        "../resources/finngen/{finn_sample}_hg19lifted.tsv.gz"
    output:
        "../results/filtered/nominal/{finn_sample}_filtered.tsv"
    params:
        threshold=5e-02
    shell:
        "zcat {input} | cut -f 3,4,7,9,20,21 | sed 's/chr//g' | awk '{{if($3 <= {params.threshold}) {{print $5, $6, $1, $2, $3, $4}}}}' > {output}"

rule nominal_filter_significant_uk:
    input:
        "../resources/UKB/{uk_sample}.gwas.imputed_v3.both_sexes.tsv.bgz"
    output:
        "../results/filtered/nominal/{uk_sample}_filtered.tsv"
    params:
        threshold=5e-02,
    shell:
        "zcat {input} | grep \"false\" | awk '{{if($12 <= {params.threshold}) {{print}}}}' > {output}"

rule filter_significant_meta:
    input:
        "../resources/meta/{meta_sample}1.TBL.gz"
    output:
        "../results/filtered/{meta_sample}_filtered.tsv"
    params:
        threshold=THRESHOLD
    shell:
        "zcat {input} | awk '{{if($10 <= {params.threshold}) {{print}}}}' > {output}"
