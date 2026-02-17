##########################################################################
### Author: Hannah L. Dixon
### Dec 2025
##########################################################################
### Script for HES1 expression across transcriptomic datasets
########################################################################## 
### Paper 1: Schmidt et al 2024 (MCF7 paper treated with Palbo vs non treated cells)

library(tidyverse)
library(tximport)
library(ggpubr)
library(dplyr)

# Assuming files have this structure: 
# /Paper_1/kallisto_output/SRR24978566_trimmed/abundance.tsv
# Set wd to the folder that contains kallisto_output

setwd("/wd")

# List all abundance.tsv files 
files <- list.files(path = ".", pattern = "abundance.tsv$", recursive = TRUE, full.names = TRUE)
print(files)

# Extract SRR IDs automatically from the folder names
sample_names <- basename(dirname(files))
names(files) <- sample_names

# Define conditions
condition <- c(rep("Palbociclib", 3), rep("Control", 3))

# Create a metadata dataframe
samples <- data.frame(
  sample = sample_names,
  condition = condition,
  file = files
)

print(samples)

# Import TPM values using tximport 
# tximport automatically extracts TPMs from kallisto outputs
txi <- tximport(
  files = samples$file,
  type = "kallisto",
  txOut = TRUE
)

# Extract TPM matrix
tpm <- txi$abundance
colnames(tpm) <- samples$sample

# Combine TPMs with metadata 
tpm_long <- as.data.frame(tpm) %>%
  rownames_to_column("transcript") %>%
  pivot_longer(-transcript, names_to = "sample", values_to = "TPM") %>%
  left_join(samples, by = "sample")

# Example plot
# Compare TPM distributions across conditions for all genes
ggplot(tpm_long, aes(x = condition, y = TPM, fill = condition)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_log10() +
  theme_bw() +
  labs(title = "TPM Distribution by Condition", y = "TPM (log10)", x = "")


# Define the gene or transcript of interest
gene_to_plot <- "ENST00000232424.4" #HES1 201 1575bp 280 aa 

# gene_to_plot <- "ENST00000476918.1"# HES1 1009 no protein retained intron


# Filter your data for HES1 canonical isofrom
gene_data <- tpm_long %>%
  filter(transcript == gene_to_plot)

condition_cols <- c("Control" = "#0072B2", "Palbociclib" = "#D55E00")

# Run t-test
t_test_result <- t.test(TPM ~ condition, data = gene_data)

## boxplot of raw TPMs 
ggplot(gene_data, aes(x = condition, y = TPM, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, colour = "black") +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.5, colour = "black") +
  scale_fill_manual(values = condition_cols) +
  stat_compare_means(
    comparisons = list(c("Control", "Palbociclib")),
    method = "t.test",
    label = "p.signif",
    size = 6,
    colour = "black"
  ) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 16),
    strip.text = element_text(size = 16),
    axis.title.x = element_blank()
  ) +
  labs(
    title = paste("Expression of", gene_to_plot),
    subtitle = paste0("t = ", round(t_test_result$statistic, 2),
                      ", p = ", signif(t_test_result$p.value, 3)),
    y = "mRNA expression (TPM)"
  )



# Calculate the median of the control TPM values to get fold change 
control_median <- gene_data %>%
  filter(condition == "Control") %>%
  summarise(median_TPM = median(TPM)) %>%
  pull(median_TPM)

# Add a new column with fold change values
gene_fc <- gene_data %>%
  mutate(fold_change = TPM / control_median)

# Boxplot of the fold change wihtout the stats 
ggplot(gene_fc, aes(x = condition, y = fold_change, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, colour = "black") +
  geom_jitter(width = 0.12, size = 2.2, alpha = 0.9, colour = "black") +
  scale_fill_manual(values = condition_cols) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    title = paste("Relative expression of", gene_to_plot),
    y = "Fold change vs Control",
    x = ""
  ) 

#+
#geom_hline(yintercept = 1, linetype = "dashed")  # baseline


# Calculate t-test and add stars to plot
fc_ttest <- t.test(fold_change ~ condition, data = gene_fc)
p_val <- fc_ttest$p.value

stars <- ifelse(p_val < 0.0001, "****",
                ifelse(p_val < 0.001, "***",
                       ifelse(p_val < 0.01, "**",
                              ifelse(p_val < 0.05, "*", "ns"))))

# Define y-position for line
y_max <- max(gene_fc_summary$mean_fc) * 1.15

## Barplot of the fold change with added stats
ggplot(gene_fc_summary, aes(x = condition, y = mean_fc, fill = condition)) +
  geom_bar(stat = "identity", colour = "black", width = 0.5) +
  geom_errorbar(aes(ymin = mean_fc - se, ymax = mean_fc + se),
                width = 0.15, size = 0.8) +
  scale_fill_manual(values = condition_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + 
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_blank()
  ) +
  labs(
    title = paste("Fold change of", gene_to_plot),
    y = "Fold change vs Control"
  ) +
  # Significance bar 
  geom_segment(
    aes(x = 1, xend = 2, y = y_max, yend = y_max),
    linewidth = 0.7
  ) +
  geom_segment(
    aes(x = 1, xend = 1, y = y_max, yend = y_max - y_max*0.03),
    linewidth = 0.7
  ) +
  geom_segment(
    aes(x = 2, xend = 2, y = y_max, yend = y_max - y_max*0.03),
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = 1.5,
    y = y_max + y_max*0.05,
    label = stars,
    size = 8
  )
