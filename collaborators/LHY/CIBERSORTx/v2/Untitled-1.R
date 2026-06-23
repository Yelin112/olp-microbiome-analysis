# =============================================================================
# CIBERSORTx Overlap and Similarity Analysis - English Version (Enhanced)
# Includes: Cell Composition Stacked Bar Chart + All Previous Analyses
# Output Format: PDF
# =============================================================================

# Load required packages
required_packages <- c(
  "tidyverse",
  "pheatmap",
  "ggdendro",
  "igraph",
  "ggrepel",
  "circlize",
  "RColorBrewer",
  "gridExtra",
  "vegan"
)

new_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]
if (length(new_packages)) {
  install.packages(new_packages)
}

library(tidyverse)
library(pheatmap)
library(ggdendro)
library(igraph)
library(ggrepel)
library(circlize)
library(RColorBrewer)
library(gridExtra)
library(vegan)

# =============================================================================
# Settings
# =============================================================================

output_dir <- "overlap_similarity_analysis_english_2"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

cat(
  "================================================================================\n"
)
cat("     CIBERSORTx Overlap & Similarity Analysis (with Statistics)\n")
cat(
  "================================================================================\n\n"
)

# =============================================================================
# 1. Data Loading and Preprocessing
# =============================================================================
cat("📊 Loading data...\n")

job2 <- read_csv("CIBERSORTx_Job2_Results.csv", show_col_types = FALSE)
job4 <- read_csv("CIBERSORTx_Job4_Results.csv", show_col_types = FALSE)

# UPDATED: Remove Myocytes (contamination)
cat("  ⚠️  Removing Myocytes (muscle contamination)...\n")

# Get all cell types except Myocytes
all_cell_types <- colnames(job2)[2:17]
cell_types <- all_cell_types[all_cell_types != "Myocytes"]

cat(sprintf("  Original cell types: %d\n", length(all_cell_types)))
cat(sprintf("  After removal: %d\n", length(cell_types)))
cat(sprintf("  Removed: %s\n\n", setdiff(all_cell_types, cell_types)))

# Remove Myocytes from data
job2 <- job2 |>
  select(-Myocytes) |>
  mutate(Group = str_extract(Mixture, "^[A-Za-z]+"))

job4 <- job4 |>
  select(-Myocytes) |>
  mutate(Group = str_extract(Mixture, "^[A-Za-z]+"))

# Calculate group means
job2_means <- job2 |>
  group_by(Group) |>
  summarise(across(all_of(cell_types), mean, na.rm = TRUE))

job4_means <- job4 |>
  group_by(Group) |>
  summarise(across(all_of(cell_types), mean, na.rm = TRUE))

# Calculate changes
osf_change <- (job2_means |>
  filter(Group == "OSF") |>
  select(-Group) |>
  as.numeric()) -
  (job2_means |> filter(Group == "Normal") |> select(-Group) |> as.numeric())
names(osf_change) <- cell_types

pbs_mean <- job4_means |>
  filter(Group == "PBS") |>
  select(-Group) |>
  as.numeric()

job4_changes <- job4_means |>
  filter(Group != "PBS") |>
  column_to_rownames("Group") |>
  sweep(2, pbs_mean, "-")

# Merge all changes
all_changes <- bind_rows(
  data.frame(Group = "OSF", t(osf_change)),
  job4_changes |> rownames_to_column("Group")
) |>
  column_to_rownames("Group")

cat("✓ Data preparation complete\n")
cat("  Number of groups:", nrow(all_changes), "\n")
cat("  Number of cell types (after removal):", ncol(all_changes), "\n\n")
# =============================================================================
# NEW: 0. Cell Composition Stacked Bar Chart
# =============================================================================

cat("📊 NEW Analysis: Cell Composition Stacked Bar Chart...\n")

# Prepare composition data (absolute abundance, not changes)
composition_data <- bind_rows(
  job2_means |> filter(Group %in% c("Normal", "OSF")),
  job4_means
) |>
  pivot_longer(
    cols = all_of(cell_types),
    names_to = "Cell_Type",
    values_to = "Abundance"
  ) |>
  group_by(Group) |>
  mutate(Percentage = Abundance / sum(Abundance) * 100) |>
  ungroup()

# Reorder groups for better visualization
group_order <- c("Normal", "OSF", "PBS", "Am", "Pm", "Av", "Ai", "AB", "B")
composition_data <- composition_data |>
  mutate(Group = factor(Group, levels = group_order))

# Calculate mean percentage for ordering cell types
cell_type_order <- composition_data |>
  group_by(Cell_Type) |>
  summarise(mean_pct = mean(Percentage)) |>
  arrange(desc(mean_pct)) |>
  pull(Cell_Type)

composition_data <- composition_data |>
  mutate(Cell_Type = factor(Cell_Type, levels = cell_type_order))

# Create color palette (16 distinct colors)
color_palette <- c(
  "#E41A1C",
  "#377EB8",
  "#4DAF4A",
  "#984EA3",
  "#FF7F00",
  "#FFFF33",
  "#A65628",
  "#F781BF",
  "#999999",
  "#66C2A5",
  "#FC8D62",
  "#8DA0CB",
  "#E78AC3",
  "#A6D854",
  "#FFD92F",
  "#E5C494"
)
names(color_palette) <- cell_type_order

# Custom theme (adapted from user's alluvial theme)
theme_composition <- theme_bw() +
  theme(
    # Legend on the right
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12, color = "black"),
    legend.text = element_text(face = "bold", size = 10, color = "black"),
    legend.key.size = unit(0.5, "cm"),

    # Remove grid lines
    panel.grid = element_blank(),

    # Panel spacing
    panel.spacing.x = unit(0, units = "cm"),

    # Strip (facet label) styling
    strip.background = element_rect(
      color = "black",
      fill = "white",
      linewidth = 0.8
    ),
    strip.placement = "outside",
    strip.text.x = element_text(size = 14, face = "bold"),

    # Axis lines
    axis.line.y.left = element_line(color = "black", linewidth = 0.8),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.8),

    # Text styling - all bold
    axis.text = element_text(face = "bold", size = 12, color = "black"),
    axis.title = element_text(face = "bold", size = 14, colour = "black"),

    # Ticks
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.3),

    # Plot margins
    plot.margin = margin(10, 10, 10, 10)
  )

# Create stacked bar chart
p_composition <- ggplot(
  composition_data,
  aes(x = Group, y = Percentage, fill = Cell_Type)
) +
  geom_bar(
    stat = "identity",
    position = "stack",
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = color_palette, name = "Cell Type") +
  labs(
    title = "Cell Composition Across Groups",
    x = "Groups",
    y = "Percentage (%)"
  ) +
  theme_composition +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))

# Save the plot
ggsave(
  file.path(output_dir, "0_cell_composition_stacked_bar.pdf"),
  p_composition,
  width = 12,
  height = 8,
  device = "pdf"
)

cat("✓ Cell composition chart saved\n\n")

# Save composition data
composition_summary <- composition_data |>
  select(Group, Cell_Type, Percentage) |>
  pivot_wider(names_from = Cell_Type, values_from = Percentage)

write_csv(
  composition_summary,
  file.path(output_dir, "cell_composition_percentages.csv")
)

# =============================================================================
# 2. Hierarchical Clustering + Statistics
# =============================================================================

cat("🌳 Analysis 1: Hierarchical Clustering + Statistics...\n")

dist_matrix <- dist(all_changes, method = "euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")

# Cophenetic correlation
coph_corr <- cor(dist_matrix, cophenetic(hc))
cat(sprintf("  Cophenetic correlation: %.3f\n", coph_corr))
cat("  (>0.75 = Excellent, >0.60 = Good, <0.60 = Poor)\n\n")

# Plot dendrogram
pdf(
  file.path(output_dir, "1_hierarchical_clustering_dendrogram.pdf"),
  width = 12,
  height = 6
)

par(mar = c(5, 5, 4, 2))
plot(
  hc,
  main = sprintf(
    "Hierarchical Clustering Dendrogram (Cophenetic r = %.3f)",
    coph_corr
  ),
  xlab = "Groups",
  ylab = "Distance (Ward's method)",
  cex.main = 1.5,
  cex.lab = 1.2,
  cex.axis = 1.1,
  hang = -1
)

rect.hclust(hc, k = 3, border = "red")
text(
  x = par("usr")[2] * 0.02,
  y = par("usr")[4] * 0.95,
  labels = sprintf("Clustering quality: %.3f", coph_corr),
  adj = c(0, 1),
  cex = 1.1,
  col = "blue"
)

dev.off()

cat("✓ Dendrogram saved\n\n")

# =============================================================================
# 3. Principal Component Analysis (PCA) + Statistics
# =============================================================================

cat("📉 Analysis 2: PCA + Statistics...\n")

# Data cleaning
cat("  Cleaning data...\n")
valid_cols <- sapply(all_changes, function(x) {
  !any(is.na(x)) && !any(is.infinite(x)) && var(x, na.rm = TRUE) > 1e-10
})

all_changes_pca <- all_changes[, valid_cols]
cat(sprintf(
  "  Using %d/%d cell types for PCA\n",
  sum(valid_cols),
  length(valid_cols)
))

# PCA
pca_result <- prcomp(all_changes_pca, center = TRUE, scale. = FALSE)
pca_scores <- as.data.frame(pca_result$x[, 1:2])
pca_scores$Group <- rownames(pca_scores)

explained_var <- summary(pca_result)$importance[2, ]

cat(sprintf("  PC1 variance explained: %.2f%%\n", explained_var[1] * 100))
cat(sprintf("  PC2 variance explained: %.2f%%\n", explained_var[2] * 100))
cat(sprintf("  Cumulative variance: %.2f%%\n", sum(explained_var[1:2]) * 100))

# Loadings analysis
loadings <- pca_result$rotation[, 1:2]
top_pc1 <- names(sort(abs(loadings[, 1]), decreasing = TRUE))[
  1:min(3, nrow(loadings))
]
top_pc2 <- names(sort(abs(loadings[, 2]), decreasing = TRUE))[
  1:min(3, nrow(loadings))
]

cat("\n  Top contributors to PC1:\n")
for (i in seq_along(top_pc1)) {
  cat(sprintf(
    "    %d. %s (loading = %.3f)\n",
    i,
    top_pc1[i],
    loadings[top_pc1[i], 1]
  ))
}

cat("\n  Top contributors to PC2:\n")
for (i in seq_along(top_pc2)) {
  cat(sprintf(
    "    %d. %s (loading = %.3f)\n",
    i,
    top_pc2[i],
    loadings[top_pc2[i], 2]
  ))
}

# PERMANOVA test
cat("\n  PERMANOVA test (between-group differences):\n")
perm_result <- tryCatch(
  {
    adonis2(
      all_changes_pca ~ factor(rownames(all_changes_pca)),
      permutations = 999,
      method = "euclidean"
    )
  },
  error = function(e) {
    cat("  ⚠️  PERMANOVA test failed\n")
    NULL
  }
)

if (!is.null(perm_result)) {
  perm_r2 <- perm_result$R2[1]
  perm_p <- perm_result$`Pr(>F)`[1]
  cat(sprintf("  R² = %.3f\n", perm_r2))
  cat(sprintf("  p-value = %.4f %s\n", perm_p, ifelse(perm_p < 0.05, "*", "")))
} else {
  perm_r2 <- NA
  perm_p <- NA
}

# Save statistics
pca_stats <- data.frame(
  Component = paste0("PC", 1:length(explained_var)),
  Variance_Explained = explained_var * 100,
  Cumulative_Variance = cumsum(explained_var) * 100
)
write_csv(pca_stats, file.path(output_dir, "pca_variance_explained.csv"))

loadings_df <- as.data.frame(loadings) |>
  rownames_to_column("Cell_Type") |>
  arrange(desc(abs(PC1)))
write_csv(loadings_df, file.path(output_dir, "pca_loadings.csv"))

# Plot PCA
p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point(
    aes(
      size = ifelse(Group == "OSF", 8, 5),
      color = Group,
      shape = ifelse(Group == "OSF", 18, 16)
    ),
    alpha = 0.8
  ) +
  geom_label_repel(
    aes(label = Group, fontface = ifelse(Group == "OSF", "bold", "plain")),
    size = 4,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = 20
  ) +
  scale_shape_identity() +
  scale_size_identity() +
  scale_color_manual(
    values = c(
      "OSF" = "#E06681",
      "Am" = "#8087E2",
      "Pm" = "#90EE90",
      "Av" = "#FFB6C1",
      "Ai" = "#87CEEB",
      "AB" = "#DDA0DD",
      "B" = "#F0E68C"
    )
  ) +
  labs(
    title = if (!is.na(perm_p)) {
      sprintf(
        "Principal Component Analysis\n(PERMANOVA: R² = %.3f, p = %.3f)",
        perm_r2,
        perm_p
      )
    } else {
      "Principal Component Analysis"
    },
    x = sprintf("PC1 (%.1f%%)", explained_var[1] * 100),
    y = sprintf("PC2 (%.1f%%)", explained_var[2] * 100)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggsave(
  file.path(output_dir, "2_pca_analysis.pdf"),
  p_pca,
  width = 10,
  height = 8,
  device = "pdf"
)

cat("\n✓ PCA analysis complete\n\n")

# =============================================================================
# 4. Correlation Matrix + Significance Test
# =============================================================================

cat("🔥 Analysis 3: Correlation Matrix + Significance Test...\n")

n_groups <- nrow(all_changes)
corr_matrix <- matrix(NA, n_groups, n_groups)
pval_matrix <- matrix(NA, n_groups, n_groups)
rownames(corr_matrix) <- colnames(corr_matrix) <- rownames(all_changes)
rownames(pval_matrix) <- colnames(pval_matrix) <- rownames(all_changes)

for (i in 1:n_groups) {
  for (j in 1:n_groups) {
    test_result <- cor.test(
      as.numeric(all_changes[i, ]),
      as.numeric(all_changes[j, ]),
      method = "pearson"
    )
    corr_matrix[i, j] <- test_result$estimate
    pval_matrix[i, j] <- test_result$p.value
  }
}

sig_matrix <- pval_matrix < 0.05

cat(sprintf(
  "  Significant pairs (p < 0.05): %d / %d\n",
  sum(sig_matrix[upper.tri(sig_matrix)]),
  sum(upper.tri(sig_matrix))
))

# Save correlation and p-values
corr_df <- as.data.frame(corr_matrix) |> rownames_to_column("Group")
write_csv(corr_df, file.path(output_dir, "correlation_matrix.csv"))

pval_df <- as.data.frame(pval_matrix) |> rownames_to_column("Group")
write_csv(pval_df, file.path(output_dir, "correlation_pvalues.csv"))

# Plot heatmap with significance markers
pdf(
  file.path(output_dir, "3_correlation_clustered_heatmap.pdf"),
  width = 11,
  height = 10
)

sig_labels <- matrix("", n_groups, n_groups)
sig_labels[pval_matrix < 0.001] <- "***"
sig_labels[pval_matrix >= 0.001 & pval_matrix < 0.01] <- "**"
sig_labels[pval_matrix >= 0.01 & pval_matrix < 0.05] <- "*"

pheatmap(
  corr_matrix,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  breaks = seq(-1, 1, length.out = 101),
  display_numbers = matrix(
    sprintf("%.2f%s", corr_matrix, sig_labels),
    nrow = nrow(corr_matrix)
  ),
  fontsize_number = 9,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  main = "Clustered Correlation Heatmap\n(* p<0.05, ** p<0.01, *** p<0.001)",
  fontsize = 11,
  fontsize_row = 11,
  fontsize_col = 11,
  border_color = "white",
  cellwidth = 55,
  cellheight = 55
)

dev.off()

cat("✓ Correlation heatmap saved\n\n")

# =============================================================================
# 5. Correlation Network
# =============================================================================

cat("🕸️  Analysis 4: Correlation Network...\n")

threshold <- 0.3
edges <- data.frame()

for (i in 1:(nrow(corr_matrix) - 1)) {
  for (j in (i + 1):nrow(corr_matrix)) {
    corr_val <- corr_matrix[i, j]
    pval <- pval_matrix[i, j]
    if (abs(corr_val) > threshold) {
      edges <- rbind(
        edges,
        data.frame(
          from = rownames(corr_matrix)[i],
          to = rownames(corr_matrix)[j],
          weight = abs(corr_val),
          correlation = corr_val,
          pvalue = pval,
          significant = pval < 0.05
        )
      )
    }
  }
}

cat(sprintf(
  "  Network edges: %d (threshold |r| > %.2f)\n",
  nrow(edges),
  threshold
))
cat(sprintf("  Significant edges: %d\n", sum(edges$significant)))

write_csv(edges, file.path(output_dir, "network_edges.csv"))

# Create network
g <- graph_from_data_frame(
  edges,
  directed = FALSE,
  vertices = rownames(corr_matrix)
)

# Plot network
pdf(file.path(output_dir, "4_correlation_network.pdf"), width = 12, height = 10)

set.seed(42)
layout <- layout_with_fr(g, niter = 1000)

node_colors <- ifelse(V(g)$name == "OSF", "#FFD700", "#90EE90")
node_sizes <- ifelse(V(g)$name == "OSF", 30, 20)

edge_colors <- ifelse(E(g)$correlation > 0, "#E06681", "#8087E2")
edge_lty <- ifelse(E(g)$significant, 1, 2)

par(mar = c(1, 1, 3, 1))
plot(
  g,
  layout = layout,
  vertex.color = node_colors,
  vertex.size = node_sizes,
  vertex.frame.color = "black",
  vertex.frame.width = 2,
  vertex.label.color = "black",
  vertex.label.font = 2,
  vertex.label.cex = 1.2,
  edge.color = edge_colors,
  edge.width = E(g)$weight * 5,
  edge.lty = edge_lty,
  edge.curved = 0.2,
  main = sprintf(
    "Correlation Network (threshold |r| > %.1f, solid=p<0.05)",
    threshold
  )
)

legend(
  "topright",
  legend = c("Positive", "Negative", "Significant", "Non-sig", "OSF", "Models"),
  col = c("#E06681", "#8087E2", "black", "black", "#FFD700", "#90EE90"),
  lty = c(1, 1, 1, 2, NA, NA),
  pch = c(NA, NA, NA, NA, 21, 21),
  pt.bg = c(NA, NA, NA, NA, "#FFD700", "#90EE90"),
  pt.cex = 2,
  lwd = c(3, 3, 2, 2, NA, NA),
  bty = "n",
  cex = 1
)

dev.off()

cat("✓ Network plot saved\n\n")

# =============================================================================
# 6. Overlap Analysis + Fisher's Exact Test
# =============================================================================

cat("📊 Analysis 5: Overlap Analysis + Fisher's Test...\n")

sig_threshold <- 0.01

# Identify significant cells
significant_cells <- list()
for (group in rownames(all_changes)) {
  sig_cells <- colnames(all_changes)[abs(all_changes[group, ]) > sig_threshold]
  significant_cells[[group]] <- sig_cells
}

cat("\nSignificantly changed cell types (|change| > 0.01):\n")
cat(strrep("-", 70), "\n")
for (group in names(significant_cells)) {
  cells <- significant_cells[[group]]
  cat(sprintf(
    "%-5s: %2d types - %s\n",
    group,
    length(cells),
    paste(sort(cells), collapse = ", ")
  ))
}
cat(strrep("-", 70), "\n\n")

# Overlap with OSF + Fisher's test
osf_cells <- significant_cells[["OSF"]]
overlap_data <- data.frame()

cat("Fisher's Exact Test (overlap significance):\n")
cat(strrep("-", 70), "\n")

for (group in names(significant_cells)) {
  if (group != "OSF") {
    group_cells <- significant_cells[[group]]
    overlap <- intersect(osf_cells, group_cells)

    # Fisher's exact test
    a <- length(overlap)
    b <- length(setdiff(osf_cells, group_cells))
    c <- length(setdiff(group_cells, osf_cells))
    d <- length(cell_types) - a - b - c

    fisher_test <- fisher.test(matrix(c(a, b, c, d), nrow = 2))

    cat(sprintf(
      "  %s: OR = %.2f, p = %.4f %s\n",
      group,
      fisher_test$estimate,
      fisher_test$p.value,
      ifelse(fisher_test$p.value < 0.05, "*", "")
    ))

    overlap_data <- rbind(
      overlap_data,
      data.frame(
        Group = group,
        OSF_Only = b,
        Overlap = a,
        Group_Only = c,
        Overlap_Rate_OSF = a / length(osf_cells),
        Overlap_Rate_Group = a / length(group_cells),
        Fisher_OR = as.numeric(fisher_test$estimate),
        Fisher_pvalue = fisher_test$p.value
      )
    )
  }
}
cat(strrep("-", 70), "\n\n")

write_csv(
  overlap_data,
  file.path(output_dir, "overlap_analysis_with_stats.csv")
)

# Plot overlap analysis
pdf(file.path(output_dir, "5_overlap_analysis.pdf"), width = 12, height = 10)

par(mfrow = c(2, 1), mar = c(5, 5, 4, 2))

# Plot 1: Stacked barplot
overlap_matrix <- as.matrix(overlap_data[, c(
  "OSF_Only",
  "Overlap",
  "Group_Only"
)])
rownames(overlap_matrix) <- overlap_data$Group

barplot(
  t(overlap_matrix),
  beside = FALSE,
  col = c("#E06681", "#9467BD", "#8087E2"),
  border = "black",
  main = "Overlap of Significantly Changed Cell Types: OSF vs Models",
  xlab = "Model Groups",
  ylab = "Number of Cell Types",
  cex.main = 1.3,
  cex.lab = 1.2,
  cex.names = 1.1,
  ylim = c(0, max(rowSums(overlap_matrix)) * 1.2)
)

# Add significance stars
for (i in 1:nrow(overlap_data)) {
  if (overlap_data$Fisher_pvalue[i] < 0.05) {
    text(
      x = i,
      y = sum(overlap_matrix[i, ]) + 0.5,
      labels = "*",
      cex = 2,
      col = "red"
    )
  }
}

legend(
  "topright",
  legend = c("OSF only", "Overlap", "Model only", "* p<0.05"),
  fill = c("#E06681", "#9467BD", "#8087E2", NA),
  border = c("black", "black", "black", NA),
  pch = c(NA, NA, NA, 42),
  col = c(NA, NA, NA, "red"),
  bty = "n",
  cex = 1.1
)

# Plot 2: Overlap rate + p-values
overlap_data_sorted <- overlap_data |> arrange(desc(Overlap_Rate_OSF))
bp <- barplot(
  overlap_data_sorted$Overlap_Rate_OSF * 100,
  names.arg = overlap_data_sorted$Group,
  col = ifelse(overlap_data_sorted$Fisher_pvalue < 0.05, "#9467BD", "#CCCCCC"),
  border = "black",
  main = "Overlap Rate with OSF (* Fisher's p<0.05)",
  xlab = "Model Groups",
  ylab = "Overlap Rate (%)",
  cex.main = 1.2,
  cex.lab = 1.2,
  cex.names = 1.1,
  horiz = TRUE,
  las = 1,
  ylim = c(0, max(overlap_data_sorted$Overlap_Rate_OSF * 100) * 1.2)
)

for (i in 1:nrow(overlap_data_sorted)) {
  text(
    x = overlap_data_sorted$Overlap_Rate_OSF[i] * 100 + 2,
    y = bp[i],
    labels = sprintf(
      "%.1f%% (p=%.3f)%s",
      overlap_data_sorted$Overlap_Rate_OSF[i] * 100,
      overlap_data_sorted$Fisher_pvalue[i],
      ifelse(overlap_data_sorted$Fisher_pvalue[i] < 0.05, "*", "")
    ),
    pos = 4,
    cex = 0.9
  )
}

dev.off()

cat("✓ Overlap analysis saved\n\n")

# =============================================================================
# 7. Chord Diagram
# =============================================================================

cat("🎵 Analysis 6: Chord Diagram...\n")

chord_threshold <- 0.4
chord_data <- data.frame()

for (i in 1:(nrow(corr_matrix) - 1)) {
  for (j in (i + 1):nrow(corr_matrix)) {
    corr_val <- corr_matrix[i, j]
    if (abs(corr_val) > chord_threshold) {
      chord_data <- rbind(
        chord_data,
        data.frame(
          from = rownames(corr_matrix)[i],
          to = rownames(corr_matrix)[j],
          value = abs(corr_val)
        )
      )
    }
  }
}

if (nrow(chord_data) > 0) {
  pdf(file.path(output_dir, "6_chord_diagram.pdf"), width = 12, height = 12)

  groups <- unique(c(chord_data$from, chord_data$to))
  grid_colors <- c(
    "OSF" = "#E06681",
    "Am" = "#8087E2",
    "Pm" = "#90EE90",
    "Av" = "#FFB6C1",
    "Ai" = "#87CEEB",
    "AB" = "#DDA0DD",
    "B" = "#F0E68C"
  )[groups]

  circos.clear()
  circos.par(start.degree = 90, gap.degree = 4)

  chordDiagram(
    chord_data,
    grid.col = grid_colors,
    transparency = 0.5,
    directional = 0,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.1)
  )

  circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      sector.name = get.cell.meta.data("sector.index")
      circos.text(
        mean(xlim),
        ylim[1],
        sector.name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 1.2,
        font = ifelse(sector.name == "OSF", 2, 1)
      )
    },
    bg.border = NA
  )

  title(
    sprintf(
      "Chord Diagram: Inter-group Correlations (|r| > %.1f)",
      chord_threshold
    ),
    cex.main = 1.5
  )

  dev.off()
  circos.clear()

  cat("✓ Chord diagram saved\n\n")
} else {
  cat("⚠️  No correlations above threshold, skipping chord diagram\n\n")
}

# =============================================================================
# 8. Generate Statistical Summary Report
# =============================================================================

cat("📝 Generating statistical summary report...\n\n")

# Safe retrieval of variables
pc1_contrib_text <- if (exists("top_pc1") && length(top_pc1) > 0) {
  paste(
    sprintf(
      "    %d. %s (loading = %.3f)",
      1:length(top_pc1),
      top_pc1,
      loadings[top_pc1, 1]
    ),
    collapse = "\n"
  )
} else {
  "    Data unavailable"
}

pc2_contrib_text <- if (exists("top_pc2") && length(top_pc2) > 0) {
  paste(
    sprintf(
      "    %d. %s (loading = %.3f)",
      1:length(top_pc2),
      top_pc2,
      loadings[top_pc2, 2]
    ),
    collapse = "\n"
  )
} else {
  "    Data unavailable"
}

permanova_text <- if (exists("perm_result") && !is.null(perm_result)) {
  sprintf(
    "    R² = %.3f\n    F = %.2f\n    p-value = %.4f %s",
    perm_result$R2[1],
    perm_result$F[1],
    perm_result$`Pr(>F)`[1],
    ifelse(
      perm_result$`Pr(>F)`[1] < 0.001,
      "***",
      ifelse(
        perm_result$`Pr(>F)`[1] < 0.01,
        "**",
        ifelse(perm_result$`Pr(>F)`[1] < 0.05, "*", "")
      )
    )
  )
} else {
  "    Not executed or failed"
}

fisher_text <- if (exists("overlap_data") && nrow(overlap_data) > 0) {
  paste(
    sprintf(
      "    %s: overlap=%.1f%%, OR=%.2f, p=%.4f %s",
      overlap_data$Group,
      overlap_data$Overlap_Rate_OSF * 100,
      overlap_data$Fisher_OR,
      overlap_data$Fisher_pvalue,
      ifelse(
        overlap_data$Fisher_pvalue < 0.001,
        "***",
        ifelse(
          overlap_data$Fisher_pvalue < 0.01,
          "**",
          ifelse(overlap_data$Fisher_pvalue < 0.05, "*", "")
        )
      )
    ),
    collapse = "\n"
  )
} else {
  "    Data unavailable"
}

summary_text <- sprintf(
  "
================================================================================
           Overlap & Similarity Analysis - Statistical Summary
================================================================================

【0. Cell Composition Overview】
  Total groups analyzed: %d (Normal, OSF, PBS, Am, Pm, Av, Ai, AB, B)
  Cell types profiled: %d
  Composition data saved to: cell_composition_percentages.csv
  Visualization: 0_cell_composition_stacked_bar.pdf

【1. Hierarchical Clustering Statistics】
  Method: Ward's D2
  Distance: Euclidean
  Cophenetic correlation: %.3f
  Interpretation: %s
  
【2. Principal Component Analysis (PCA)】
  PC1 variance explained: %.2f%%
  PC2 variance explained: %.2f%%
  Cumulative variance: %.2f%%
  
  Top PC1 contributors:
%s
  
  Top PC2 contributors:
%s
  
  PERMANOVA test:
%s
    
【3. Correlation Analysis】
  Total pairs: %d
  Significant pairs (p<0.05): %d (%.1f%%)
  Positive correlations: %d
  Negative correlations: %d
  
【4. Overlap Analysis (Fisher's Exact Test)】
  OSF significant cells: %d
  
  Overlap with each group:
%s

【5. Network Analysis】
  Threshold: |r| > %.2f
  Total edges: %d
  Significant edges: %d (p<0.05)
  Mean correlation: %.3f
  
【Statistical Interpretation】
  * p < 0.05: Significant
  ** p < 0.01: Highly significant
  *** p < 0.001: Extremely significant
  
  Cophenetic correlation:
    > 0.75: Excellent clustering
    0.60-0.75: Good clustering
    < 0.60: Poor clustering
    
  PERMANOVA:
    Tests if between-group differences are significant
    R²: Proportion of variance explained by groups
    
  Fisher's Exact Test:
    Tests if overlap exceeds random expectation
    OR > 1: Overlap higher than expected
    OR < 1: Overlap lower than expected

【Output Files】
  Visualizations (PDF format):
    • 0_cell_composition_stacked_bar.pdf - NEW!
    • 1_hierarchical_clustering_dendrogram.pdf
    • 2_pca_analysis.pdf
    • 3_correlation_clustered_heatmap.pdf
    • 4_correlation_network.pdf
    • 5_overlap_analysis.pdf
    • 6_chord_diagram.pdf
  
  Data tables:
    • cell_composition_percentages.csv - NEW!
    • correlation_matrix.csv
    • correlation_pvalues.csv
    • pca_variance_explained.csv
    • pca_loadings.csv
    • overlap_analysis_with_stats.csv
    • network_edges.csv

================================================================================
✅ Statistical analysis complete!
================================================================================
",
  # Cell composition
  length(group_order),
  length(cell_types),

  # Hierarchical clustering
  coph_corr,
  ifelse(
    coph_corr > 0.75,
    "Excellent",
    ifelse(coph_corr > 0.60, "Good", "Poor")
  ),

  # PCA
  explained_var[1] * 100,
  explained_var[2] * 100,
  sum(explained_var[1:2]) * 100,

  pc1_contrib_text,
  pc2_contrib_text,

  permanova_text,

  # Correlation
  sum(upper.tri(corr_matrix)),
  sum(pval_matrix[upper.tri(pval_matrix)] < 0.05),
  mean(pval_matrix[upper.tri(pval_matrix)] < 0.05) * 100,
  sum(corr_matrix[upper.tri(corr_matrix)] > 0),
  sum(corr_matrix[upper.tri(corr_matrix)] < 0),

  # Overlap
  length(osf_cells),
  fisher_text,

  # Network
  threshold,
  nrow(edges),
  sum(edges$significant),
  mean(edges$weight)
)

cat(summary_text)

writeLines(
  summary_text,
  file.path(output_dir, "STATISTICAL_ANALYSIS_SUMMARY.txt")
)

cat("\n✓ Statistical summary saved\n")
cat(
  "\n================================================================================\n"
)
cat("              🎉 All analyses complete (with statistics)!\n")
cat(
  "================================================================================\n"
)
