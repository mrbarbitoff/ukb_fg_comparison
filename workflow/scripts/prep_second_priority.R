library(tidyverse)

read_tsv(snakemake@input[['corr_table']]) %>%
  filter(priority == 2) %>%
  transmute(Finn_code, UK_code) %>%
  write_tsv(file = snakemake@output[['second_priority']], col_names = F)
