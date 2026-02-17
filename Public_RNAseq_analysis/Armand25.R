##########################################################################
### Author: Hannah L. Dixon
### Dec 2025
##########################################################################
### Script for HES1 expression across transcriptomic datasets
########################################################################## 
### Paper 2: Armand et al., 2025


library(tximport)
library(dplyr)
library(tidyr)
library(tibble)
library(ggsignif)
library(dplyr)
library(ggplot2)

setwd("/Ollie")
## Create a table
samples <- tibble::tribble(
  ~Run,           ~condition,
  "SRR30931915",  "Palbocilib + INX-315 + fulverstrant",
  "SRR30931916",  "Palbocilib + INX-315 + fulverstrant",
  "SRR30931917",  "INX-315 + fulverstrant",
  "SRR30931918",  "INX-315 + fulverstrant",
  "SRR30931919",  "Palbocilib + fulverstrant",
  "SRR30931920",  "Palbocilib + fulverstrant",
  "SRR30931921",  "DMSO",
  "SRR30931922",  "DMSO"
)

# Be in your working directory and have a folder called tpm_data where all the abundance files (created by kallisto) are in 
files <- file.path("tpm_data", paste0(samples$Run, "_abundance.tsv"))
names(files) <- samples$Run
txi <- tximport(files, type = "kallisto", txOut = TRUE)

tpm_long <- txi$abundance %>%
  as.data.frame() %>%
  rownames_to_column("transcript") %>%
  pivot_longer(-transcript, names_to = "Run", values_to = "TPM") %>%
  left_join(samples, by = "Run")

# Check both HES1 isoforms
gene_to_plot <- "ENST00000232424.4" # Canonical isoform
gene_to_plot <- "ENST00000476918.1"

df_gene <- tpm_long %>% filter(transcript == gene_to_plot)

# mean + SD per condition
df_sum <- df_gene %>%
  group_by(condition) %>%
  summarise(
    mean_TPM = mean(TPM),
    sd_TPM   = sd(TPM),
    sem_TPM  = sd(TPM) / sqrt(n()),
    .groups  = "drop"
  )

# define order 
df_sum$condition <- factor(df_sum$condition,
                           levels = c("Palbocilib + fulverstrant",
                                      "DMSO",
                                      "Palbocilib + INX-315 + fulverstrant",
                                      "INX-315 + fulvestrant")
)

df_gene$condition <- factor(df_gene$condition,
                            levels = levels(df_sum$condition)
)

# Calculate SEM as there's only 2 samples
ggplot() +
  geom_col(
    data = df_sum,
    aes(x = condition, y = mean_TPM, fill = condition),
    width = 0.6, alpha = 1, show.legend = FALSE
  ) +
  geom_point(
    data = df_gene,
    aes(x = condition, y = TPM),
    size = 2, color = "black", alpha = 0.5
  ) +
  geom_errorbar(
    data = df_sum,
    aes(x = condition,
        ymin = mean_TPM - sem_TPM,
        ymax = mean_TPM + sem_TPM),
    width = 0.15, linewidth = 0.8
  ) +
  labs(
    x = NULL,
    y = "TPM",
    title = "ENST00000232424.4 expression"
  ) +
    scale_fill_manual(values = c(
    "DMSO"                                  = "#0072B2",
    "Palbocilib + fulverstrant"             = "#D55E00",  
    "Palbocilib + INX-315 + fulverstrant"   = "#D55E00"   
  )) +
  
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.title.y = element_text(size = 20)
  ) +
  
  expand_limits(y = max(df_sum$mean_TPM + df_sum$sd_TPM) * 1.1)



### Computing fold change for plot 

df_sum <- df_gene %>%
  group_by(condition) %>%
  summarise(
    mean_TPM = mean(TPM),
    sd_TPM = sd(TPM),
    n = n(),
    sem_TPM = sd_TPM / sqrt(n)
  )

fc <- df_sum %>%
  mutate(
    fold_change = mean_TPM / mean_TPM[condition == "DMSO"],
    fc_sd = fold_change * sqrt( (sd_TPM/mean_TPM)^2 +
                                  (sd_TPM[condition == "DMSO"] / mean_TPM[condition == "DMSO"])^2 )
  )

ggplot(fc, aes(x = condition, y = fold_change, fill = condition)) +
  geom_col(width = 0.5, alpha = 1, color = "black") +
  geom_errorbar(aes(ymin = fold_change - fc_sd,
                    ymax = fold_change + fc_sd),
                width = 0.2, size = 0.7) +
  geom_hline(yintercept = 1, linewidth = 0.9, linetype = "dotted") +
  labs(
    x = NULL,
    y = "Fold Change vs DMSO",
    title = "HES1 gene-level fold change (combined isoforms)"
  ) +
  scale_fill_manual(values = c(
    "Palbocilib + fulverstrant"             = "#D55E00",
    "DMSO"                                  = "#009E73",
    "Palbocilib + INX-315 + fulverstrant"   = "#D55E00",
    "INX-315 + fulverstrant"                = "#009E73"
  )) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size=16)
  ) +
  expand_limits(y = max(fc$fold_change + fc$fc_sd) * 1.1)+
  scale_x_discrete(limits = c(
    "Palbocilib + fulverstrant",
    "DMSO",
    "Palbocilib + INX-315 + fulverstrant",
    "INX-315 + fulverstrant"
  ))

