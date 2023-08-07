import json
from collections import defaultdict as dd

rsid_lookup = dd(str)
rsid_lookup["NA"] = "NA"

with open(snakemake.input.SNP_template) as SNP_template:
    SNP_template.readline()  # skip header
    for line in SNP_template:
        varid, _, finn_rsids, uk_rsids = line.strip().split()[0:4]
        if finn_rsids == "NA" and uk_rsids == "NA":
            rsid_lookup[varid] = varid
        elif finn_rsids == "NA":
            for rsid in uk_rsids.split(","):
                rsid_lookup[rsid] = varid
        elif uk_rsids == "NA":
            for rsid in finn_rsids.split(","):
                rsid_lookup[rsid] = varid
        else:
            for rsid in uk_rsids.split(","):
                rsid_lookup[rsid] = varid
            for rsid in finn_rsids.split(","):
                rsid_lookup[rsid] = varid


with open(snakemake.output.rsid_lookup, 'w') as fp:
    json.dump(rsid_lookup, fp)
