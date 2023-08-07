library(tidyverse)

finn_data <- read_tsv(snakemake@input[['finn_example_dataset']]) %>%
  transmute(varid = paste(str_replace(new_chr, "chr", ''), new_coord, ref, alt, sep = ":"),
            af_alt,
            af_ctrl = 1 - af_alt,
            Finn_rsid = rsids) %>%
  transmute(varid, Finn_MAF = pmin(af_ctrl, af_alt), Finn_rsid) %>%
  na.omit()

ukb_data <- read_tsv(snakemake@input[['uk_variants']]) %>%
  transmute(varid = variant,
            UKB_rsid = rsid,
            UKB_MAF = minor_AF)

joined <- finn_data %>%
  full_join(ukb_data, by = 'varid') %>%
  distinct(varid, .keep_all = T)

joined %>%
  write_tsv(snakemake@output[['SNP_template']])
