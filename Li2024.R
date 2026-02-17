##########################################################################
### Author: Hannah L. Dixon
### Dec 2025
##########################################################################
### Script for HES1 expression across transcriptomic datasets
########################################################################## 
### Paper: Li et al 2024 


library(tidyverse)
library(tximport)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(ggsignif)

setwd("Li_2024/cancer_files_raw")

# Read sample information
samples <- read_csv("/cancer_files_raw/mapping_ids_patients.csv")   # path to your CSV file
head(samples)

# Create named vector of abundance file paths by Run ID, each folder name = Run column values
files <- file.path("1_raw", samples$Run, "abundance.tsv")
names(files) <- samples$Run
head(files)

# Import quantification using tximport
txi <- tximport(files, type = "kallisto", txOut = TRUE)

# TPM matrix
tpm <- txi$abundance

# Set columns to run IDs
colnames(tpm) <- samples$Run   # Must match metadata column

# Convert to long format
tpm_long <- tpm %>%
  as.data.frame() %>%
  rownames_to_column("transcript") %>%
  pivot_longer(-transcript, names_to = "Run", values_to = "TPM") %>%
  left_join(samples, by = "Run")  # merges subtype info

# library(DESeq2)
# 
# dds <- DESeqDataSetFromTximport(
#   txi,
#   colData = samples,
#   design = ~ subtype
# )
# 
# dds <- DESeq(dds)
# res <- results(dds)

gene_to_plot <- "ENST00000232424.4" # Canonical HES1 isoform 

# gene_to_plot <- "ENST00000909760.1"

# isoforms <- c("ENST00000232424.4", "ENST00000476918.1")
# df_iso <- tpm_long %>%
#   filter(transcript %in% isoforms)
# 

tpm_gene <- txi$abundance[gene_to_plot, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("transcript")

colnames(tpm_gene)[2:ncol(tpm_gene)] <- samples$Run


tpm_gene_long <- tpm_gene %>%
  pivot_longer(-transcript, names_to = "Run", values_to = "TPM") %>%
  left_join(samples, by = "Run")

tpm_gene_long <- tpm_gene_long %>%
  mutate(pair = parse_number(`Sample Name`))

tpm_gene_long <- tpm_gene_long %>%
  group_by(pair) %>% 
  mutate(pairID = cur_group_id()) %>%
  ungroup()


tpm_gene_long$subtype <- factor(
  tpm_gene_long$subtype,
  levels = c("Normal", "Luminal A")
)

ggplot(tpm_gene_long, aes(x = subtype, y = TPM, color = subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, linewidth = 1) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Normal" = "#1f78b4", "Luminal A" = "#e31a1c")) +
  theme_bw(base_size = 14) +
  ggtitle(paste("TPM expression of", gene_to_plot)) +
  xlab("") +
  ylab("TPM")


df_gene <- tpm_long %>%
  filter(transcript == gene_to_plot) %>%
  filter(subtype %in% c("Normal", "Luminal A")) %>%
  mutate(subtype = factor(subtype, levels = c("Normal", "Luminal A")))

wilcox.test(df_gene$TPM, df_gene$subtype, paired = TRUE)
df_gene %>% arrange(Patient, subtype)
wilcox.test(TPM ~ subtype, data = df_gene, paired = TRUE)


df_wide <- df_gene %>%
  select(Patient, subtype, TPM) %>%
  pivot_wider(names_from = subtype, values_from = TPM)
wilcox.test(df_wide$`Luminal A`, df_wide$Normal, paired = TRUE)


ggplot(df_gene, aes(x = subtype, y = TPM, fill = subtype)) +
  geom_boxplot(width = 0.5, show.legend = FALSE, alpha = 1, outlier.shape = NA) +
  geom_point(aes(alpha = 0.1), size = 2, show.legend = FALSE) +
  
  # Dotted lines connecting groupes samples
  geom_line(aes(group = Patient), color = "black", alpha = 0.7, linetype = "dotted") +
  
  # Significance test 
  geom_signif(
    comparisons = list(c("Normal", "Luminal A")),
    map_signif_level = TRUE,
    test = "wilcox.test",
    y_position = max(df_gene$TPM) * 1.05,
    tip_length = 0.02
  ) +
  
  labs(
    y = "TPM",
    x = NULL,
    title = paste("TPM expression of", gene_to_plot)
  ) +
  
  scale_fill_manual(values = c(
    "Normal" = "#D55E00",
    "Luminal A" = "#0072B2"
  )) +
  
  scale_y_continuous() +
  theme(
    axis.title      = element_text(size = 16),
    axis.text       = element_text(size = 16),
    strip.text      = element_text(size = 16),
    axis.text.x     = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title      = element_text(size = 18, hjust = 0.5)
  )


