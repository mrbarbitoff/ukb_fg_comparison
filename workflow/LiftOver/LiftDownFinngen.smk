import pandas as pd

chain_file = config["CHAIN"]  # hg38ToHg19.over.chain.gz

finn_ids = config["FINN_IDS"]  # ../LDAK/finn/finn_ids_uniq.txt
# original_coords = config["ORIGINAL_COORDS"]
finn_sample_table = pd.read_table(finn_ids, header=None, names=("Samples", ))
finn_samples = list(finn_sample_table.Samples.unique())

onsuccess:
    print("Workflow successful")

wildcard_constraints:
    finn_sample = "|".join(finn_samples)

rule all:
    input:
        expand("../finngen/R6/{finn_sample}_hg19lifted.tsv.gz", finn_sample=finn_samples)
    output:
        touch(".liftdown_done")

rule sum2bed:
    input:
        "../finngen/R6/finngen_R6_{finn_sample}.gz"
    output:
        temp("{finn_sample}.bed")
    shell:
        "zcat {input} | awk 'NR>1 {{print \"chr\"$1, $2, $2+1, $1\":\"$2}}' | sed 's/chr23/chrX/g'> {output}"
        # unzip the file, add 'chr' to chromosome number, create single nucleotide intervals
        # use chr:coord notation as ID to track unlifted coordinates
        # ad hoc replacement of chr23 to chrX

rule liftover:
    input:
        inbed="{finn_sample}.bed",
        chain=chain_file,
    output:
        lifted=temp("{finn_sample}.lifted.bed"),
        unlifted="{finn_sample}.unlifted.bed",
    shell:
        "./liftOver {input.inbed} {input.chain} {output.lifted} {output.unlifted}"

rule merge_lifted:
    input:
        lifted_bed="{finn_sample}.lifted.bed",
        table="../finngen/R6/finngen_R6_{finn_sample}.gz",
    output:
        hg19lifted="../finngen/R6/{finn_sample}_hg19lifted.tsv.gz"
    script:
        "merge_lifted.py"
#
# rule paste_new_coords:
#     input:
#         hg19lifted="original2lifted_lookup.tsv",
#         table="../finngen/finngen_R5_{finn_sample}.gz"
#     output:
#         "{finn_sample}_hg19lifted.tsv.gz"
#     shell:
#         "paste -d '\\t' {input.hg19lifted} <(zcat {input.table}) | gzip > {output}"
