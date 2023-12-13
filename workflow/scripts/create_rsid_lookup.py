import json
from collections import defaultdict as dd

rsid_lookup = dd(str)
rsid_lookup["NA"] = "NA"

with open(snakemake.input.SNP_template) as SNP_template:
    for line in SNP_template:
        varid, _, _, rsid = line.strip().split()[0:4]
        if rsid == "NA":
            rsid_lookup[varid] = varid
        else:
            rsid_lookup[rsid] = varid

with open(snakemake.output.rsid_lookup, 'w') as fp:
    json.dump(rsid_lookup, fp)
