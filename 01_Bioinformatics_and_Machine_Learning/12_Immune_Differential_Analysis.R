# ==============================================================================
# Script: 11_Immune_Differential_Analysis.R
# Purpose: Compare 22 immune cell infiltration fractions between Control and 
#          Treat (Pulmonary Fibrosis) groups using Wilcoxon rank-sum tests.
# Relevant Figures: Figure 5A, Figure 5B
# ==============================================================================

suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(broom)
  library(ggsci)
})

inputFile <- "CIBERSORT_Results.csv"

# 1. Data Ingestion & Group Label Assignment -----------------------------------
rt <- read.csv(inputFile, header = TRUE, row.names = 1, check.names = FALSE)
cell_types <- colnames(rt)[1:22]
rt_cells <- rt[, cell_types, drop = FALSE]

# Extract group status from sample naming
group_labels <- ifelse(grepl("Control|_con$", rownames(rt_cells), ignore.case = TRUE), "Control", "Treat")
rt_cells$Group <- factor(group_labels, levels = c("Control", "Treat"))
rt_cells$Sample <- rownames(rt_cells)

# Convert to long-format dataframe
data_long <- melt(rt_cells, id.vars = c("Sample", "Group"), 
                  variable.name = "ImmuneCell", value.name = "Fraction")
data_long <- data_long[!is.na(data_long$Group), ]

# 2. Statistical Testing across Immune Populations -----------------------------
pvalue_table <- data_long %>%
  group_by(ImmuneCell) %>%
  do(tidy(wilcox.test(Fraction ~ Group, data = .))) %>%
  select(ImmuneCell, p.value) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value)

write.csv(pvalue_table, file = "immune_differential_statistics.csv", row.names = FALSE)

# Summary table (Mean, Median, SD)
summary_table <- data_long %>%
  group_by(ImmuneCell, Group) %>%
  summarise(
    Mean = mean(Fraction, na.rm = TRUE),
    Median = median(Fraction, na.rm = TRUE),
    SD = sd(Fraction, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(summary_table, file = "immune_fraction_summary.csv", row.names = FALSE)

# 3. Comprehensive Differential Boxplot (Figure 5A) ----------------------------
myColors <- c("Control" = "#4A90E2", "Treat" = "#E74C3C")

pdf("Figure5A_Immune_Differential_Boxplot.pdf", width = 11, height = 6.5)
p_box <- ggplot(data_long, aes(x = ImmuneCell, y = Fraction, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.8) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
             size = 0.8, alpha = 0.4, color = "black") +
  stat_compare_means(
    aes(group = Group),
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    )
  ) +
  scale_fill_manual(values = myColors) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    axis.title.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(title = "Infiltration Distribution of 22 Immune Cell Populations", y = "Relative Fraction")

print(p_box)
dev.off()

# 4. Individual Cell Boxplots for Key Disrupted Populations (Figure 5B) --------
output_dir <- "individual_immune_boxplots"
if (!dir.exists(output_dir)) dir.create(output_dir)

for (cell in unique(data_long$ImmuneCell)) {
  sub_df <- subset(data_long, ImmuneCell == cell)
  y_limit <- max(sub_df$Fraction, na.rm = TRUE) * 1.2
  
  p_single <- ggplot(sub_df, aes(x = Group, y = Fraction, fill = Group)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(shape = 21, color = "black", alpha = 0.6, width = 0.15, size = 2) +
    scale_fill_manual(values = myColors) +
    stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.4, label.y = y_limit * 0.9) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(title = cell, x = NULL, y = "Infiltration Proportion")
  
  ggsave(file.path(output_dir, paste0(gsub("[ /]", "_", cell), "_boxplot.pdf")),
         p_single, width = 4, height = 4.5)
}