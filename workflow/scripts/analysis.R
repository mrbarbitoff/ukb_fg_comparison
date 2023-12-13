library(tidyverse)
library(reshape2)
library(vioplot)
library(viridis)
library(egg)
library(ggpubr)
library(ggupset)
library(tidytext)
library(ggvenn)
library(RColorBrewer)

# default text size
text_size = 20

SNP_data <- read_tsv(snakemake@input[["filtered_SNP_table"]], col_types = "cnncnnncicciccicccccccc") %>%
  mutate(meta_index = str_detect(is_index, "index_meta"),
         ukbb_index = str_detect(is_index, "index_ukbb"),
         finn_index = str_detect(is_index, "index_finn")) %>%
  mutate(Finn_MAF = pmin(Finn_MAF, 1 - Finn_MAF), UKB_MAF = pmin(UKB_MAF, 1 - UKB_MAF))

nominal_data <- read_tsv(snakemake@input[["nominal_sign_table"]], col_types = "cnnccici", col_names = c("varid", "Finn_MAF", "UKB_MAF", "rsid", "nominal_phenonames_finn", "nom_pleio_finn", "nominal_phenonames_uk", "nom_pleio_uk"))


SNP_data <- SNP_data %>%
  left_join(nominal_data, by = c("varid", "Finn_MAF", "UKB_MAF", "rsid"))

index_SNPs <- SNP_data %>%
  filter(!is.na(is_index)) %>%
  filter(meta_index | ukbb_index | finn_index)

total_bed <- read_tsv(snakemake@input[["updated_table"]], col_names = c("chr", "start", "stop", "index_field")) %>%
  mutate(clump_id = paste(chr, start, stop, sep = ":"),
         meta_index = str_detect(index_field, "meta"),
         ukbb_index = str_detect(index_field, "ukbb"),
         finn_index = str_detect(index_field, "finn"))

trait_data <- read_tsv(snakemake@input[["updated_table"]])

intersection_data <- read_tsv(snakemake@input[["clump_intersections"]])

new_clump_data <- trait_data %>%
  select(fg_phenotype, finn_clump_n, uk_clump_n, meta_clump_n) %>%
  left_join(intersection_data, by = "fg_phenotype") %>%
  select(fg_phenotype, contains("clump"))


# Function for fig2 and fig3

plot_battery <- function(df1, df2, col_palette, labels)
{

# labels = c(deparse(substitute(df1)), deparse(substitute(df2)))

df1$source <- labels[1]
df2$source <- labels[2]

df <- rbind(df1, df2)

UK_MAF_test <- df %>%
  melt(measure = c("Finn_MAF", "UKB_MAF"), variable = "SNP_source") %>%
  wilcox.test(value ~ source, data = ., subset = SNP_source == "UKB_MAF")

Finn_MAF_test <- df %>%
  melt(measure = c("Finn_MAF", "UKB_MAF"), variable = "SNP_source") %>%
  wilcox.test(value ~ source, data = ., subset = SNP_source == "Finn_MAF")

MAF_finn_p <- paste("Wilcox test p-value", format(Finn_MAF_test$p.value, digits = 3))
MAF_uk_p <- paste("Wilcox test p-value", format(UK_MAF_test$p.value, digits = 3))

MAF_plot <- df %>%
  melt(measure = c("Finn_MAF", "UKB_MAF"), variable = "SNP_source") %>%
  ggplot(aes(x = source, y = value, fill = source))+
  geom_violin(width=1)+
  geom_boxplot(width=0.1, color="grey", fill="white", outlier.shape = NA)+
  scale_fill_manual(values = col_palette, name = element_blank())+
  theme_bw()+
  theme(text = element_text(size = 20), axis.text.x=element_blank(), strip.background = element_blank())+
  labs(y = "MAF", x = "")+
  facet_wrap("SNP_source", labeller = labeller("SNP_source" = c("Finn_MAF" = paste("FinnGen MAF", MAF_finn_p, sep = "\n"),
                                        "UKB_MAF" = paste("UKB MAF", MAF_uk_p, sep = "\n"))))

MAF_diff_test <- df %>%
  mutate(maf_diff = abs(UKB_MAF - Finn_MAF)/pmax(UKB_MAF, Finn_MAF)) %>%
  wilcox.test(maf_diff ~ source, data = .)

MAF_diff_plot <- df %>%
  mutate(maf_diff = abs(UKB_MAF - Finn_MAF)/pmax(UKB_MAF, Finn_MAF)) %>%
  ggplot(aes(x = source, y = maf_diff, fill = source))+
  geom_violin(width=1)+
  geom_boxplot(width=0.1, color="grey", fill="white", outlier.shape = NA)+
  scale_fill_manual(values = col_palette, name = element_blank()) +
  labs(caption = paste("Wilcox test p-value", format(MAF_diff_test$p.value, digits = 3))) +
  theme_bw()+
  theme(text = element_text(size = 20), axis.text.x=element_blank())+
  labs(y = "Difference in MAF", x = "")

pleio_plot <- df %>%
  melt(measure = c("pleio_finn", "pleio_uk", "pleio_meta"), variable = "SNP_source") %>%
  ggplot(aes(x = source, y = value, fill = source))+
  geom_violin(width=1, outlier.shape = NA)+
  geom_boxplot(width=0.1, color="grey", fill="white", outlier.shape = NA)+
  scale_fill_manual(values = col_palette, name = element_blank())+
  theme_bw()+
  theme(text = element_text(size = 20), axis.text.x=element_blank(), strip.background = element_blank())+
  ylim(0, 30)+
  labs(x = "", y = "Pleiotropy")+
  scale_x_discrete(labels = c("FinnGen", "UKB", "Meta-analysis"))+
  facet_wrap("SNP_source", labeller = labeller("SNP_source" = c("pleio_finn" = "FinnGen",
                                        "pleio_uk" = "UKB",
                                        "pleio_meta" = "Meta-analysis")))

beta_test_uk <- df %>%
  melt(id = c("varid", "source"), measure = c("beta_finn", "beta_UK"), variable = "SNP_source") %>%
  wilcox.test(value ~ source, data = ., subset = SNP_source == "beta_UK")

beta_test_finn <- df %>%
  melt(id = c("varid", "source"), measure = c("beta_finn", "beta_UK"), variable = "SNP_source") %>%
  wilcox.test(value ~ source, data = ., subset = SNP_source == "beta_finn")

beta_uk_p <- paste("Wilcox test p-value", format(beta_test_uk$p.value, digits = 3))
beta_finn_p <- paste("Wilcox test p-value", format(beta_test_finn$p.value, digits = 3))

beta_plot <- df %>%
  melt(id = c("varid", "source"), measure = c("beta_finn", "beta_UK"), variable = "SNP_source") %>%
  ggplot(aes(x = source, y = value, fill = source))+
  geom_violin(width=1)+
  geom_boxplot(width=0.1, color="grey", fill="white", outlier.shape = NA)+
  scale_fill_manual(values = col_palette, name = element_blank())+
  theme_bw()+
  theme(text = element_text(size = 20),
        axis.text.x=element_blank(),
        strip.background = element_blank())+
  labs(x = "", y = "Effect size")+
  ylim(0, 1.1)+
  facet_wrap("SNP_source", labeller = labeller("SNP_source" = c("beta_finn" = paste("FinnGen effect size", beta_finn_p, sep = "\n"),
                                        "beta_UK" = paste("UKB effect size", beta_uk_p, sep = "\n"))))

LD_test <- df %>%
  wilcox.test(LD_score ~ source, data = .)

LD_plot <- df %>%
  ggplot(aes(x = source, y = LD_score, fill = source))+
  geom_violin(width=1)+
  geom_boxplot(width=0.1, color="grey", fill="white", outlier.shape = NA)+
  scale_fill_manual(values = col_palette, name = element_blank())+
  labs(x = "", caption = paste("Wilcox test p-value", format(LD_test$p.value, digits = 3)), y = "LD score") +
  theme_bw()+
  theme(text = element_text(size = 20), axis.text.x=element_blank())

return(list(MAF_plot, MAF_diff_plot, pleio_plot, beta_plot, LD_plot))
}


# shared clump boxplot
unnested_data <- index_SNPs %>%
  unnest_tokens(clump_finn, clump_finn, token = "regex", pattern = ",") %>%
  unnest_tokens(clump_uk, clump_uk, token = "regex", pattern = ",") %>%
  unnest_tokens(clump_meta, clump_meta, token = "regex", pattern = ",")

area1 = unnested_data %>% filter(finn_index) %>% pull(varid)
area2 = unnested_data %>% filter(ukbb_index) %>% pull(varid)
area3 = unnested_data %>% filter(meta_index) %>% pull(varid)


fig1_1 <- ggvenn(list(FinnGen = area1, UKB = area2, 'Meta-analysis' = area3))+
  theme(text = element_text(size = 10))

fig1_2 <- new_clump_data %>%
  na.omit() %>%
  mutate(across(c(intersect_clump_fm_n, intersect_clump_fu_n, intersect_clump_um_n), ~ . + intersect_clump_fum_n)) %>%
  mutate(`Meta-analysis` = meta_clump_n, UKB = uk_clump_n, FinnGen = finn_clump_n) %>%
  mutate(FinnGen_UKB = intersect_clump_fu_n, `FinnGen_Meta-analysis` = intersect_clump_fm_n, `UKB_Meta-analysis` = intersect_clump_um_n, `FinnGen_UKB_Meta-analysis` = intersect_clump_fum_n) %>%
  melt(measure = c("FinnGen_UKB", "FinnGen_Meta-analysis", "UKB_Meta-analysis", "FinnGen_UKB_Meta-analysis", "FinnGen", "UKB", "Meta-analysis"), id = "fg_phenotype", variable = "SNP_source") %>%
  filter(value > 0) %>%
  ggplot(aes(SNP_source, value))+
  geom_violin(width=1, outlier.shape = NA)+
  geom_boxplot(width=0.1, color="grey", fill="white", outlier.shape = NA)+
  axis_combmatrix(sep = "_")+
  theme_bw()+
  ylim(0, 50)+
  theme(text = element_text(size = 20))+
  labs(title="Number associated genomic regions by association source",
        x ="", y = "Number of associated regions")

# Fig 1
pdf(snakemake@output[[1]], 20, 10)

ggarrange(fig1_1, fig1_2)

dev.off()

# Meta-reproducible
finn_meta_repr_SNPs <- index_SNPs %>%
  filter(str_detect(repr_fm, "perfect|partial")) %>%
  filter(finn_index)

uk_meta_repr_SNPs <- index_SNPs %>%
  filter(str_detect(repr_um, "perfect|partial")) %>%
  filter(ukbb_index)

uk_non_repr <- index_SNPs %>% filter(pleio_meta == 0, ukbb_index)
finn_non_repr <- index_SNPs %>% filter(pleio_meta == 0, finn_index)

meta_repr_SNPs <- rbind(finn_meta_repr_SNPs, uk_meta_repr_SNPs) %>%
  distinct()

meta_non_repr_SNPs <- rbind(uk_non_repr, finn_non_repr) %>%
  distinct()

meta_repr_plots <- plot_battery(meta_repr_SNPs, meta_non_repr_SNPs, c("magenta", "cyan"), labels = c("Reproducing variants", "Non-reproducing variants"))

# Fig 2
pdf(snakemake@output[[2]], 20, 10)

ggarrange(meta_repr_plots[[1]], meta_repr_plots[[2]], meta_repr_plots[[4]], meta_repr_plots[[5]], common.legend = T, legend = "bottom")

dev.off()

meta_SNPs <- index_SNPs %>%
  filter(pleio_finn == 0) %>%
  filter(pleio_meta > 0) %>%
  filter(pleio_uk == 0) %>%
  filter(meta_index)

non_meta_SNPs <- index_SNPs %>%
  filter(pleio_finn > 0 | pleio_uk > 0) %>%
  anti_join(meta_SNPs, by = "varid")

meta_plots <- plot_battery(meta_SNPs, non_meta_SNPs, col_palette = c("yellow", "green"), labels = c("Meta-analysis specific SNPs", "Nonmeta-analysis specific SNPs"))

# Fig 3
pdf(snakemake@output[[3]], 20, 10)

ggarrange(meta_plots[[1]], meta_plots[[2]], meta_plots[[4]], meta_plots[[5]], common.legend = T, legend = "bottom")

dev.off()

# Supp1

pdf(snakemake@output[[4]], 20, 10)

trait_data %>%

  ggplot(aes(fg_n_cases, ukbb_n_cases))+
  geom_point()+
  geom_abline(aes(intercept=0, slope=1))+
  coord_fixed(clip = "off")+
  theme_bw()+
  theme(text = element_text(size = 30))+
  labs(y = "Cases in UKB", x = "Cases in FinnGen")+
  ggtitle("Cases across biobanks")

dev.off()

# Supp2

trait_data_for_plots <- trait_data %>%
  na.omit() %>%
  filter(meta_clump_n > 0, uk_clump_n > 0, finn_clump_n > 0)

p1 <- trait_data_for_plots %>%
  ggplot(aes(finn_clump_n, uk_clump_n))+
  geom_hex()+
  theme_bw()+
  theme(text = element_text(size = text_size))+
  labs(x = "FinnGen associated regions", y = "UKB associated regions")

p2 <- trait_data_for_plots %>%
  ggplot(aes(meta_clump_n, uk_clump_n))+
  geom_hex()+
  theme_bw()+
  theme(text = element_text(size = text_size))+
  labs(x = "Meta-analysis associated regions", y = "UKB associated regions")

p3 <- trait_data_for_plots %>%
  ggplot(aes(meta_clump_n, finn_clump_n))+
  geom_hex()+
  theme_bw()+
  theme(text = element_text(size = text_size))+
  labs(x = "Meta-analysis associated regions", y = "FinnGen associated regions")


Supp2 <- ggarrange(p1, p2, p3)

pdf(snakemake@output[[5]], 20, 15)

Supp2

dev.off()

#Supp3

pleio_plot <- function(df1, df2, col_palette, labels)
{

df1$source <- labels[1]
df2$source <- labels[2]

df <- rbind(df1, df2)

pleio_plot <- df %>%
  melt(measure = c("pleio_finn", "pleio_uk", "pleio_meta"), variable = "SNP_source") %>%
  mutate(is_pleio = value > 1) %>%
  ggplot(aes(x = SNP_source, fill = is_pleio))+
  geom_bar(position = "fill", stat = "count")+
  scale_fill_brewer(guide = guide_legend(reverse=TRUE), palette = col_palette, name = "Is variant pleiotropic?", labels = c("no", "yes"))+
  theme_bw()+
  theme(text = element_text(size = 25),
        axis.title.x=element_blank())+
  labs(y = "Proportion of SNPs")+
  scale_x_discrete(labels = c("FinnGen", "UKB", "Meta-analysis"))+

  facet_wrap("source")

return(pleio_plot)
}


pdf(snakemake@output[[6]], 30, 10)

pleio1 <- pleio_plot(meta_repr_SNPs, meta_non_repr_SNPs, "Set1", labels = c("Reproducing variants", "Non-reproducing variants"))
pleio2 <- pleio_plot(meta_SNPs, non_meta_SNPs, "Set3", labels = c("Meta-analysis specific SNPs", "Nonmeta-analysis specific SNPs"))

ggarrange(pleio1, pleio2)

# Supp4

# Variants w/o GW significance in Finngen
UKB_GW <- SNP_data %>%
  filter(pleio_finn == 0) %>%
  filter(pleio_uk > 0)

# Nominal significance and listed in our pairs
UKB_GW_paired <- UKB_GW %>%
  transmute(varid, phenonames_uk, nominal_phenonames_finn) %>%
  unnest_tokens(phenonames_uk, phenonames_uk, token = "regex", pattern = ",") %>%
  unnest_tokens(nominal_phenonames_finn, nominal_phenonames_finn, token = "regex", pattern = ",") %>%
  mutate(nominal_phenonames_finn = str_replace(nominal_phenonames_finn, pattern = ":.*", replacement = "") %>% toupper(),
         phenonames_uk = str_replace(phenonames_uk, pattern = ":.*", replacement = "") %>% toupper()) %>%
  distinct(varid, phenonames_uk, nominal_phenonames_finn, .keep_all = T) %>%
  filter(nominal_phenonames_finn == phenonames_uk)



# Variants w/o GW significance in UKB
Finn_GW <- SNP_data %>%
  filter(pleio_finn > 0) %>%
  filter(pleio_uk == 0)

# Nominal significance and listed in our pairs
Finn_GW_paired <- Finn_GW %>%
  transmute(varid, phenonames_finn, nominal_phenonames_uk) %>%
  unnest_tokens(phenonames_finn, phenonames_finn, token = "regex", pattern = ",") %>%
  unnest_tokens(nominal_phenonames_uk, nominal_phenonames_uk, token = "regex", pattern = ",") %>%
  mutate(nominal_phenonames_uk = str_replace(nominal_phenonames_uk, pattern = ":.*", replacement = "") %>% toupper(),
         phenonames_finn = str_replace(phenonames_finn, pattern = ":.*", replacement = "") %>% toupper()) %>%
  distinct(varid, phenonames_finn, nominal_phenonames_uk, .keep_all = T) %>%
  filter(phenonames_finn == nominal_phenonames_uk)


# Nominal group means that associations are not supported by nominal from other BB
stacked_finn_data <- Finn_GW %>%
  anti_join(Finn_GW_paired, by = 'varid') %>%
  mutate(group = "Nominal") %>%
  transmute(varid, group) %>%
  rbind(Finn_GW_paired %>%
          mutate(group = "Paired") %>%
          transmute(varid, group)) %>%
  mutate(source = "Finngen")

stacked_uk_data <- UKB_GW %>%
  anti_join(UKB_GW_paired, by = 'varid') %>%
  mutate(group = "Nominal") %>%
  transmute(varid, group) %>%
  rbind(UKB_GW_paired %>%
          mutate(group = "Paired") %>%
          transmute(varid, group)) %>%
  mutate(source = "UKB")

pdf(snakemake@output[[7]], 30, 10)

rbind(stacked_uk_data, stacked_finn_data) %>%
  ggplot(aes(x = source, fill = group))+
  scale_fill_brewer(palette = "Set1", name = "Group of SNPs", labels = c("Nominal association in another biobank", "No nominal association in another biobank"))+
  geom_bar(position = 'fill')+
  theme_bw()+
  theme(text = element_text(size = 30),
        axis.title.x=element_blank())+
  scale_x_discrete(labels = c("FinnGen", "UKB"))+
  ylab("Proportion of SNPs")

dev.off()

# Supp5

pdf(snakemake@output[[8]], 30, 10)
clump_area1 = total_bed %>% filter(finn_index) %>% pull(clump_id)
clump_area2 = total_bed %>% filter(ukbb_index) %>% pull(clump_id)
clump_area3 = total_bed %>% filter(meta_index) %>% pull(clump_id)

clump_venn <- ggvenn(list(FinnGen = clump_area1, UKB = clump_area2, 'Meta-analysis' = clump_area3))+
  theme(text = element_text(size = 10))

clump_venn
dev.off()
