##########################################################################
### Author: Hannah L. Dixon
### Dec 2025
##########################################################################
### Script for HES1 expression across transcriptomic datasets
########################################################################## 
### Paper: Frost et al 2021 (Cancer vs Normal expression of HES1)


library(tidyverse)
library(tidyr)
library(readr)
library(dplyr)

norm_df <- read.delim("~/TissueMeanGeneExp.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cancer_df <- read.delim("~/CancerMeanGeneExp.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# View first few rows to check
colnames(norm_df)
head(cancer_df)

gene_id <- "ENSG00000114315"
norm_gene   <- norm_df[rownames(norm_df) == gene_id, ]
cancer_gene <- cancer_df[rownames(cancer_df) == gene_id, ]

# extract the values for breast samples only 
norm_value   <- as.numeric(norm_gene["breast"])  
cancer_value <- as.numeric(cancer_gene["BRCA"])  

# join data for plotting 
expr_values <- c(Normal = norm_value, Cancer = cancer_value)

# simple barplot
barplot(expr_values,
        col = c("#0072B2", "#D55E00"),
        main = "Expression of ENSG00000114315",
        ylab = "Mean expression (FPKM)")


df_expr <- data.frame(
  group = names(expr_values),
  value = as.numeric(expr_values)
)

df_expr$group <- factor(df_expr$group, levels = c("Normal", "Cancer"))

ggplot(df_expr, aes(x = group, y = value, fill = group)) +
  geom_col(width = 0.5,              # narrow bars
           color = "black",          # black border
           alpha = 0.9) +            # solid fill
  labs(
    title = "Expression of ENSG00000114315",
    x = NULL,
    y = "Mean expression (FPKM)"
  ) +
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#D55E00"   # orange
  )) +        # grey panel background
  theme(  # remove grid lines
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    legend.position = "none"         # no legend needed for 2 bars
  ) +
  expand_limits(y = max(df_expr$value) * 1.2)  # extra space above bars
