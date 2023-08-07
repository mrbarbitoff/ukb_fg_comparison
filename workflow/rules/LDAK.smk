rule convert_ukbb:
    input:
        "../resources/UKB/{uk_sample}.gwas.imputed_v3.both_sexes.tsv.bgz"
    output:
        temp("../results/LDAK/{uk_sample}.sum")
    script:
        "../scripts/convert_UKBB_to_LDAK_summary.py"

rule convert_finn:
    input:
        "../resources/finngen/{finn_sample}_hg19lifted.tsv.gz",
    output:
        temp("../results/LDAK/{finn_sample}.sum")
    params:
        n_samples = lambda wildcards: finn_manifest_dict[wildcards.finn_sample]
    script:
        "../scripts/convert_finngen_to_LDAK_summary.py"

rule LDAK:
    input:
        summary="../results/LDAK/{sample}.sum",
        taggings=TAGGINGS,
    output:
        "../results/LDAK/{sample}.hers"
    params:
        output_prefix = "../results/LDAK/{sample}"
    threads:
        THREADS
    run:
        cmd = f"./{LDAK} --max-threads {threads} --sum-hers {params.output_prefix} --summary {input.summary} --tagfile {input.taggings} --check-sums NO"
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()

        if p.returncode != 0:  # if some error happens
            print("LDAK error occured")
            print(stderr.decode())
            with open(output[0], "w") as outfile:
                outfile.write("Her_All\tNaN")
