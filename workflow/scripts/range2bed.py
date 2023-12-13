import json

# load rsid lookup

with open(snakemake.input.rsid_lookup, 'r') as fp:
    rsid_lookup = json.load(fp)

with open(snakemake.input.range) as range_file, open(snakemake.output.bed, "w") as bed_file:
    range_file.readline()  # skip header
    for line in range_file:

        # disregard initial clump size by default. If you want actual ranges use this:
        # range = line.split()[4]
        # chrom, coords = range.split(":")
        # start, stop = (int(x) for x in coords.split(".."))

        name = line.split()[1]
        index_SNP = rsid_lookup[name]
        index_coord = int(index_SNP.split(":")[1])
        start = index_coord - snakemake.params.flank
        stop = index_coord + snakemake.params.flank

        start = max(0, start)
        print('\t'.join((chrom, str(start), str(stop), snakemake.params.current_dataset)), file=bed_file)
