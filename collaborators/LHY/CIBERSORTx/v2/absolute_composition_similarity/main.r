# =============================================================================
# CIBERSORTx Absolute Composition Similarity Analysis
# Direct comparison of cell composition between OSF and modeling groups
# Methods: Bray-Curtis, Aitchison, Cosine, PCoA, NMDS
# =============================================================================

library(tidyverse)
library(vegan) # Bray-Curtis, PERMANOVA, NMDS
library(ape) # PCoA
library(compositions) # Aitchison distance (CLR transformation)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(officer)
library(rvg)

# =============================================================================
# Settings
# =============================================================================

# setwd("~/help_others/LHY/CIBERSORTx")

output_dir <- "absolute_composition_similarity"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

group_colors <- c(
  "OSF" = "#E06681",
  "Am" = "#fdd456",
  "Pm" = "#6358cb",
  "Av" = "#4fb984",
  "Ai" = "#ed9548",
  "AB" = "#ee6161",
  "B" = "#55b7ec"
)

save_pptx_gg <- function(gg_obj, file_path) {
  doc <- officer::read_pptx() |>
    officer::add_slide(layout = "Blank", master = "Office Theme") |>
    officer::ph_with(
      value = rvg::dml(ggobj = gg_obj),
      location = officer::ph_location_fullsize()
    )
  print(doc, target = file_path)
  invisible(file_path)
}

save_pptx_base <- function(plot_expr, file_path) {
  doc <- officer::read_pptx() |>
    officer::add_slide(layout = "Blank", master = "Office Theme") |>
    officer::ph_with(
      value = rvg::dml(code = plot_expr),
      location = officer::ph_location_fullsize()
    )
  print(doc, target = file_path)
  invisible(file_path)
}

cat(
  "================================================================================\n"
)
cat("     Absolute Cell Composition Similarity Analysis\n")
cat("     (Direct comparison without controls)\n")
cat(
  "================================================================================\n\n"
)

# =============================================================================
# 1. Data Loading (Remove Myocytes)
# =============================================================================

cat("📊 Loading data...\n")

job2 <- read_csv("CIBERSORTx_Job2_Results.csv", show_col_types = FALSE)
job4 <- read_csv("CIBERSORTx_Job4_Results.csv", show_col_types = FALSE)

# Remove Myocytes (contamination)
cat("  Removing Myocytes...\n")

all_cell_types <- colnames(job2)[2:17]
cell_types <- all_cell_types[all_cell_types != "Myocytes"]

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

# Extract OSF and modeling groups (exclude Normal and PBS)
osf_data <- job2_means |> filter(Group == "OSF")
modeling_data <- job4_means |> filter(Group != "PBS")

# Combine for analysis
analysis_data <- bind_rows(osf_data, modeling_data) |>
  column_to_rownames("Group")

cat(sprintf(
  "✓ Data loaded: %d groups, %d cell types\n",
  nrow(analysis_data),
  ncol(analysis_data)
))
cat(sprintf(
  "  Groups: %s\n\n",
  paste(rownames(analysis_data), collapse = ", ")
))

# =============================================================================
# 2. Distance/Similarity Calculations
# =============================================================================

cat("📏 Calculating distance/similarity metrics...\n\n")

# ---------------------------------------------------------------------------
# 2.1 Bray-Curtis Dissimilarity (Gold standard for composition data)
# ---------------------------------------------------------------------------
cat("  1. Bray-Curtis Dissimilarity...\n")

bc_dist <- vegdist(analysis_data, method = "bray")
bc_matrix <- as.matrix(bc_dist)

# Convert to similarity (1 - dissimilarity)
bc_similarity <- 1 - bc_matrix

# Extract OSF vs each modeling group
bc_osf_vs_models <- data.frame(
  Group = rownames(analysis_data)[-1], # Exclude OSF itself
  Bray_Curtis_Similarity = bc_similarity["OSF", -1],
  Bray_Curtis_Dissimilarity = bc_matrix["OSF", -1]
)

cat(sprintf("    OSF vs modeling groups (Bray-Curtis Similarity):\n"))
for (i in 1:nrow(bc_osf_vs_models)) {
  cat(sprintf(
    "      %s: %.3f\n",
    bc_osf_vs_models$Group[i],
    bc_osf_vs_models$Bray_Curtis_Similarity[i]
  ))
}
cat("\n")

# ---------------------------------------------------------------------------
# 2.2 Aitchison Distance (Compositional data specific)
# ---------------------------------------------------------------------------
cat("  2. Aitchison Distance (CLR-based)...\n")

# CLR transformation (Centered Log-Ratio)
# Add small constant to avoid log(0)
analysis_data_clr <- analysis_data + 1e-6
clr_data <- as.data.frame(clr(acomp(analysis_data_clr)))

# Euclidean distance on CLR-transformed data = Aitchison distance
aitchison_dist <- dist(clr_data, method = "euclidean")
aitchison_matrix <- as.matrix(aitchison_dist)

# Convert to similarity (inverse distance)
max_ait <- max(aitchison_matrix)
aitchison_similarity <- 1 - (aitchison_matrix / max_ait)

# Extract OSF vs each modeling group
aitchison_osf_vs_models <- data.frame(
  Group = rownames(analysis_data)[-1],
  Aitchison_Similarity = aitchison_similarity["OSF", -1],
  Aitchison_Distance = aitchison_matrix["OSF", -1]
)

cat(sprintf("    OSF vs modeling groups (Aitchison Similarity):\n"))
for (i in 1:nrow(aitchison_osf_vs_models)) {
  cat(sprintf(
    "      %s: %.3f\n",
    aitchison_osf_vs_models$Group[i],
    aitchison_osf_vs_models$Aitchison_Similarity[i]
  ))
}
cat("\n")

# ---------------------------------------------------------------------------
# 2.3 Cosine Similarity (Pattern similarity)
# ---------------------------------------------------------------------------
cat("  3. Cosine Similarity...\n")

# Cosine similarity function
cosine_sim <- function(x, y) {
  sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
}

# Calculate cosine similarity matrix
n_groups <- nrow(analysis_data)
cosine_matrix <- matrix(NA, n_groups, n_groups)
rownames(cosine_matrix) <- colnames(cosine_matrix) <- rownames(analysis_data)

for (i in 1:n_groups) {
  for (j in 1:n_groups) {
    cosine_matrix[i, j] <- cosine_sim(
      as.numeric(analysis_data[i, ]),
      as.numeric(analysis_data[j, ])
    )
  }
}

# Extract OSF vs each modeling group
cosine_osf_vs_models <- data.frame(
  Group = rownames(analysis_data)[-1],
  Cosine_Similarity = cosine_matrix["OSF", -1]
)

cat(sprintf("    OSF vs modeling groups (Cosine Similarity):\n"))
for (i in 1:nrow(cosine_osf_vs_models)) {
  cat(sprintf(
    "      %s: %.3f\n",
    cosine_osf_vs_models$Group[i],
    cosine_osf_vs_models$Cosine_Similarity[i]
  ))
}
cat("\n")

# ---------------------------------------------------------------------------
# 2.4 Euclidean Distance (on percentages)
# ---------------------------------------------------------------------------
cat("  4. Euclidean Distance (on percentages)...\n")

# Convert to percentages
analysis_data_pct <- analysis_data / rowSums(analysis_data) * 100

euclidean_dist <- dist(analysis_data_pct, method = "euclidean")
euclidean_matrix <- as.matrix(euclidean_dist)

# Convert to similarity
max_euc <- max(euclidean_matrix)
euclidean_similarity <- 1 - (euclidean_matrix / max_euc)

# Extract OSF vs each modeling group
euclidean_osf_vs_models <- data.frame(
  Group = rownames(analysis_data)[-1],
  Euclidean_Similarity = euclidean_similarity["OSF", -1],
  Euclidean_Distance = euclidean_matrix["OSF", -1]
)

cat(sprintf("    OSF vs modeling groups (Euclidean Similarity):\n"))
for (i in 1:nrow(euclidean_osf_vs_models)) {
  cat(sprintf(
    "      %s: %.3f\n",
    euclidean_osf_vs_models$Group[i],
    euclidean_osf_vs_models$Euclidean_Similarity[i]
  ))
}
cat("\n")

# =============================================================================
# 3. Combine All Metrics
# =============================================================================

cat("📊 Combining all similarity metrics...\n")

similarity_summary <- bc_osf_vs_models |>
  left_join(aitchison_osf_vs_models, by = "Group") |>
  left_join(cosine_osf_vs_models, by = "Group") |>
  left_join(euclidean_osf_vs_models, by = "Group") |>
  mutate(
    # Average similarity across all metrics
    Average_Similarity = (Bray_Curtis_Similarity +
      Aitchison_Similarity +
      Cosine_Similarity +
      Euclidean_Similarity) /
      4
  ) |>
  arrange(desc(Average_Similarity))

write_csv(
  similarity_summary,
  file.path(output_dir, "similarity_metrics_summary.csv")
)

cat("\n")
cat("=== Similarity Ranking (Higher = More Similar to OSF) ===\n")
cat(strrep("-", 80), "\n")
for (i in 1:nrow(similarity_summary)) {
  cat(sprintf(
    "%d. %-5s: Average Similarity = %.3f\n",
    i,
    similarity_summary$Group[i],
    similarity_summary$Average_Similarity[i]
  ))
  cat(sprintf(
    "   Bray-Curtis: %.3f | Aitchison: %.3f | Cosine: %.3f | Euclidean: %.3f\n",
    similarity_summary$Bray_Curtis_Similarity[i],
    similarity_summary$Aitchison_Similarity[i],
    similarity_summary$Cosine_Similarity[i],
    similarity_summary$Euclidean_Similarity[i]
  ))
}
cat(strrep("-", 80), "\n\n")

# =============================================================================
# 4. Principal Coordinates Analysis (PCoA)
# =============================================================================

cat("📉 Performing Principal Coordinates Analysis (PCoA)...\n")

# Use Bray-Curtis distance for PCoA
pcoa_result <- pcoa(bc_dist)

# Extract coordinates
pcoa_scores <- as.data.frame(pcoa_result$vectors[, 1:2])
pcoa_scores$Group <- rownames(pcoa_scores)
colnames(pcoa_scores)[1:2] <- c("PCo1", "PCo2")

# Variance explained
var_explained <- pcoa_result$values$Relative_eig[1:2] * 100

cat(sprintf("  PCo1 variance: %.2f%%\n", var_explained[1]))
cat(sprintf("  PCo2 variance: %.2f%%\n\n", var_explained[2]))

# Plot PCoA
p_pcoa <- ggplot(pcoa_scores, aes(x = PCo1, y = PCo2)) +
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
    size = 5,
    box.padding = 0.8,
    point.padding = 0.5,
    segment.color = "grey50",
    max.overlaps = 20
  ) +
  scale_shape_identity() +
  scale_size_identity() +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Principal Coordinates Analysis (PCoA)\nBased on Bray-Curtis Dissimilarity",
    subtitle = "Direct comparison of absolute cell composition",
    x = sprintf("PCo1 (%.1f%%)", var_explained[1]),
    y = sprintf("PCo2 (%.1f%%)", var_explained[2])
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggsave(
  file.path(output_dir, "1_pcoa_bray_curtis.pdf"),
  p_pcoa,
  width = 10,
  height = 8,
  device = "pdf"
)
save_pptx_gg(p_pcoa, file.path(output_dir, "1_pcoa_bray_curtis.pptx"))

cat("✓ PCoA plot saved\n\n")

# =============================================================================
# =============================================================================
# 5. NMDS (Non-metric Multidimensional Scaling)
# =============================================================================

cat("📉 Performing NMDS...\n")

set.seed(42)
nmds_result <- tryCatch(
  suppressWarnings(
    metaMDS(analysis_data, distance = "bray", k = 2, trymax = 100)
  ),
  error = function(e) {
    cat(sprintf("  ✗ NMDS Error: %s\n", e$message))
    NULL
  }
)

if (!is.null(nmds_result)) {
  stress_value <- nmds_result$stress

  cat(sprintf("  Stress: %.3f\n", stress_value))

  # Interpret stress value with special note for small sample size
  if (stress_value < 0.05) {
    stress_quality <- "Excellent (but may indicate insufficient data)"
    cat("  ⚠️  Note: Very low stress may indicate too few samples (n=7)\n")
    cat("      NMDS works best with n > 20. Consider PCoA as primary method.\n")
  } else if (stress_value < 0.1) {
    stress_quality <- "Good"
  } else if (stress_value < 0.2) {
    stress_quality <- "Acceptable"
  } else {
    stress_quality <- "Poor"
  }

  cat(sprintf("  Quality: %s\n\n", stress_quality))

  # Extract NMDS scores
  nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
  nmds_scores$Group <- rownames(nmds_scores)

  # Plot NMDS
  p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
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
      size = 5,
      box.padding = 0.8,
      point.padding = 0.5,
      segment.color = "grey50",
      max.overlaps = 20
    ) +
    scale_shape_identity() +
    scale_size_identity() +
    scale_color_manual(values = group_colors) +
    labs(
      title = "Non-metric Multidimensional Scaling (NMDS)",
      subtitle = sprintf(
        "Based on Bray-Curtis (Stress = %.3f, n=%d samples)",
        stress_value,
        nrow(analysis_data)
      ),
      x = "NMDS1",
      y = "NMDS2",
      caption = if (stress_value < 0.05) {
        "Note: Low stress may indicate insufficient sample size. Use PCoA as primary ordination."
      } else {
        NULL
      }
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
      plot.caption = element_text(
        hjust = 0.5,
        size = 10,
        color = "red",
        face = "italic"
      ),
      legend.position = "none",
      panel.grid.major = element_line(color = "grey90"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

  ggsave(
    file.path(output_dir, "2_nmds_bray_curtis.pdf"),
    p_nmds,
    width = 10,
    height = 8,
    device = "pdf"
  )
  save_pptx_gg(p_nmds, file.path(output_dir, "2_nmds_bray_curtis.pptx"))

  cat("✓ NMDS plot saved\n")
  cat(
    "  ℹ️  Recommendation: Use PCoA as the primary ordination method for n=7\n\n"
  )
} else {
  cat("✗ NMDS failed to converge, skipping this analysis\n")
  cat("  Using PCoA only\n\n")
  stress_value <- NA
  stress_quality <- "Failed"
}

# =============================================================================
# 6. Similarity Heatmap
# =============================================================================

cat("🔥 Creating similarity heatmap...\n")

# Use Bray-Curtis similarity for heatmap
pdf(file.path(output_dir, "3_similarity_heatmap.pdf"), width = 10, height = 9)

pheatmap(
  bc_similarity,
  color = colorRampPalette(c("#20a0e5", "white", "#e73838"))(100),
  breaks = seq(0, 1, length.out = 101),
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 10,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = as.dist(1 - bc_similarity),
  clustering_distance_cols = as.dist(1 - bc_similarity),
  clustering_method = "ward.D2",
  main = "Bray-Curtis Similarity Heatmap\n(Based on Absolute Cell Composition)",
  fontsize = 12,
  fontsize_row = 12,
  fontsize_col = 12,
  border_color = "white",
  cellwidth = 50,
  cellheight = 50,
  angle_col = 45
)

dev.off()

save_pptx_base(
  {
    pheatmap(
      bc_similarity,
      color = colorRampPalette(c("#20a0e5", "white", "#e73838"))(100),
      breaks = seq(0, 1, length.out = 101),
      display_numbers = TRUE,
      number_format = "%.2f",
      fontsize_number = 10,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      clustering_distance_rows = as.dist(1 - bc_similarity),
      clustering_distance_cols = as.dist(1 - bc_similarity),
      clustering_method = "ward.D2",
      main = "Bray-Curtis Similarity Heatmap\n(Based on Absolute Cell Composition)",
      fontsize = 12,
      fontsize_row = 12,
      fontsize_col = 12,
      border_color = "white",
      cellwidth = 50,
      cellheight = 50,
      angle_col = 45
    )
  },
  file.path(output_dir, "3_similarity_heatmap.pptx")
)

cat("✓ Similarity heatmap saved\n\n")

# =============================================================================
# 7. Bar Chart: Similarity Ranking
# =============================================================================

cat("📊 Creating similarity ranking bar chart...\n")

# Prepare data for plotting (long format)
similarity_long <- similarity_summary |>
  select(
    Group,
    Bray_Curtis_Similarity,
    Aitchison_Similarity,
    Cosine_Similarity,
    Euclidean_Similarity
  ) |>
  pivot_longer(cols = -Group, names_to = "Metric", values_to = "Similarity") |>
  mutate(
    Metric = str_replace(Metric, "_Similarity", ""),
    Metric = str_replace(Metric, "_", "-"),
    Group = factor(Group, levels = similarity_summary$Group)
  )

# Plot
p_ranking <- ggplot(
  similarity_long,
  aes(x = Group, y = Similarity, fill = Metric)
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    color = "black",
    linewidth = 0.3,
    width = 0.75
  ) +
  scale_fill_brewer(palette = "Set2", name = "Similarity Metric") +
  labs(
    title = "Similarity to OSF: Multi-Metric Comparison",
    subtitle = "Based on absolute cell composition (higher = more similar)",
    x = "Modeling Groups",
    y = "Similarity Score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  ylim(0, 1)

ggsave(
  file.path(output_dir, "4_similarity_ranking_bars.pdf"),
  p_ranking,
  width = 12,
  height = 8,
  device = "pdf"
)
save_pptx_gg(p_ranking, file.path(output_dir, "4_similarity_ranking_bars.pptx"))

cat("✓ Ranking bar chart saved\n\n")

# =============================================================================
# 8. PERMANOVA Test
# =============================================================================

cat("📊 Performing PERMANOVA test...\n")

# Test if modeling groups are significantly different from OSF
# Create a factor: OSF vs Others
group_factor <- ifelse(rownames(analysis_data) == "OSF", "OSF", "Models")
group_factor <- factor(group_factor)

permanova_result <- adonis2(
  analysis_data ~ group_factor,
  permutations = 999,
  method = "bray"
)

cat("\n")
cat("PERMANOVA Results:\n")
print(permanova_result)
cat("\n")

# Pairwise PERMANOVA (OSF vs each modeling group)
cat("Pairwise PERMANOVA (OSF vs each modeling group):\n")
cat(strrep("-", 70), "\n")

pairwise_results <- data.frame()

for (model_group in rownames(analysis_data)[-1]) {
  # Select OSF and this modeling group
  pair_data <- analysis_data[c("OSF", model_group), ]
  pair_factor <- factor(rownames(pair_data))

  # PERMANOVA
  pair_perm <- adonis2(
    pair_data ~ pair_factor,
    permutations = 999,
    method = "bray"
  )

  pairwise_results <- rbind(
    pairwise_results,
    data.frame(
      Group = model_group,
      R2 = pair_perm$R2[1],
      F_value = pair_perm$F[1],
      p_value = pair_perm$`Pr(>F)`[1]
    )
  )

  cat(sprintf(
    "  OSF vs %s: R² = %.3f, F = %.2f, p = %.4f %s\n",
    model_group,
    pair_perm$R2[1],
    pair_perm$F[1],
    pair_perm$`Pr(>F)`[1],
    ifelse(pair_perm$`Pr(>F)`[1] < 0.05, "*", "")
  ))
}
cat(strrep("-", 70), "\n\n")

write_csv(
  pairwise_results,
  file.path(output_dir, "permanova_pairwise_results.csv")
)

# =============================================================================
# 9. Cell Type Contribution Analysis (SIMPER)
# =============================================================================

cat("🔍 Analyzing cell type contributions (SIMPER)...\n")

# SIMPER analysis for each modeling group vs OSF
simper_results_list <- list()

for (model_group in rownames(analysis_data)[-1]) {
  tryCatch(
    {
      pair_data <- analysis_data[c("OSF", model_group), ]
      pair_groups <- factor(
        c("OSF", model_group),
        levels = c("OSF", model_group)
      )

      simper_result <- suppressWarnings(simper(
        pair_data,
        pair_groups,
        permutations = 99
      ))
      simper_summary <- summary(simper_result)[[1]]
      n_show <- min(5, nrow(simper_summary))
      top_contributors <- simper_summary[seq_len(n_show), ]

      simper_results_list[[model_group]] <- data.frame(
        Modeling_Group = model_group,
        Cell_Type = rownames(top_contributors),
        Average_Dissimilarity = top_contributors$average,
        Contribution_Percent = top_contributors$contribution,
        Cumulative_Percent = top_contributors$cumsum
      )

      cat(sprintf(
        "\n  %s vs OSF - Top %d contributing cell types:\n",
        model_group,
        n_show
      ))
      for (i in seq_len(n_show)) {
        cat(sprintf(
          "    %d. %s: %.1f%% contribution\n",
          i,
          rownames(top_contributors)[i],
          top_contributors$contribution[i]
        ))
      }
    },
    error = function(e) {
      cat(sprintf(
        "  ⚠️  SIMPER failed for %s vs OSF: %s\n",
        model_group,
        e$message
      ))
    }
  )
}

# Combine all SIMPER results
simper_combined <- bind_rows(simper_results_list)
write_csv(
  simper_combined,
  file.path(output_dir, "simper_cell_contributions.csv")
)

cat("\n✓ SIMPER analysis complete\n\n")

# =============================================================================
# 10. Generate Summary Report
# =============================================================================

cat("📝 Generating summary report...\n")

summary_report <- sprintf(
  "
================================================================================
        Absolute Cell Composition Similarity Analysis - Summary
================================================================================

【Analysis Overview】
  Comparison Method: Direct comparison of absolute cell composition
  Groups Analyzed: OSF + 6 modeling groups (Am, Pm, Av, Ai, AB, B)
  Cell Types: %d (Myocytes removed as contamination)
  
【Similarity Metrics Used】
  1. Bray-Curtis Similarity (Ecological standard for composition data)
  2. Aitchison Similarity (Compositional data specific, CLR-based)
  3. Cosine Similarity (Pattern similarity, direction-focused)
  4. Euclidean Similarity (Classic distance on percentages)

【Similarity Ranking (Most to Least Similar to OSF)】
%s

【Best Modeling Group】
  Group: %s
  Average Similarity: %.3f
  
  Individual Metrics:
    • Bray-Curtis Similarity: %.3f
    • Aitchison Similarity: %.3f
    • Cosine Similarity: %.3f
    • Euclidean Similarity: %.3f

【PCoA Results】
  Method: Principal Coordinates Analysis (Bray-Curtis based)
  PCo1 variance explained: %.2f%%
  PCo2 variance explained: %.2f%%
  Cumulative: %.2f%%

【NMDS Results】
  Method: Non-metric Multidimensional Scaling
  Stress: %.3f
  Quality: %s

【PERMANOVA Results】
  Overall test (OSF vs all modeling groups):
    R² = %.3f
    F = %.2f
    p-value = %.4f %s
  
  Interpretation: %s

【Key Findings】
  1. The most similar modeling group to OSF is %s (avg similarity = %.3f)
  2. This suggests %s best replicates OSF's cellular immune landscape
  3. The least similar group is %s (avg similarity = %.3f)
  
【Methodological Notes】
  • This analysis compares ABSOLUTE cell composition, not changes from controls
  • Avoids confounding from Normal vs PBS baseline differences
  • Multiple metrics ensure robust similarity assessment
  • Bray-Curtis is the gold standard for compositional ecology data

【Output Files】
  Visualizations:
    • 1_pcoa_bray_curtis.pdf - PCoA ordination plot
    • 2_nmds_bray_curtis.pdf - NMDS ordination plot
    • 3_similarity_heatmap.pdf - Clustered similarity heatmap
    • 4_similarity_ranking_bars.pdf - Multi-metric comparison
  
  Data Tables:
    • similarity_metrics_summary.csv - All similarity scores
    • permanova_pairwise_results.csv - Statistical tests
    • simper_cell_contributions.csv - Cell type contributions

================================================================================
✅ Absolute composition similarity analysis complete!
================================================================================
",
  # Basic info
  length(cell_types),

  # Ranking table
  paste(
    sprintf(
      "  %d. %-5s: %.3f (BC:%.2f, Ait:%.2f, Cos:%.2f, Euc:%.2f)",
      1:nrow(similarity_summary),
      similarity_summary$Group,
      similarity_summary$Average_Similarity,
      similarity_summary$Bray_Curtis_Similarity,
      similarity_summary$Aitchison_Similarity,
      similarity_summary$Cosine_Similarity,
      similarity_summary$Euclidean_Similarity
    ),
    collapse = "\n"
  ),

  # Best group
  similarity_summary$Group[1],
  similarity_summary$Average_Similarity[1],
  similarity_summary$Bray_Curtis_Similarity[1],
  similarity_summary$Aitchison_Similarity[1],
  similarity_summary$Cosine_Similarity[1],
  similarity_summary$Euclidean_Similarity[1],

  # PCoA
  var_explained[1],
  var_explained[2],
  sum(var_explained[1:2]),

  # NMDS
  nmds_result$stress,
  ifelse(
    nmds_result$stress < 0.05,
    "Excellent",
    ifelse(
      nmds_result$stress < 0.1,
      "Good",
      ifelse(nmds_result$stress < 0.2, "Acceptable", "Poor")
    )
  ),

  # PERMANOVA
  permanova_result$R2[1],
  permanova_result$F[1],
  permanova_result$`Pr(>F)`[1],
  ifelse(permanova_result$`Pr(>F)`[1] < 0.05, "*", ""),
  ifelse(
    permanova_result$`Pr(>F)`[1] < 0.05,
    "Modeling groups are significantly different from OSF in cell composition",
    "No significant difference between OSF and modeling groups overall"
  ),

  # Key findings
  similarity_summary$Group[1],
  similarity_summary$Average_Similarity[1],
  similarity_summary$Group[1],
  similarity_summary$Group[nrow(similarity_summary)],
  similarity_summary$Average_Similarity[nrow(similarity_summary)]
)

cat(summary_report)

writeLines(
  summary_report,
  file.path(output_dir, "ABSOLUTE_COMPOSITION_ANALYSIS_SUMMARY.txt")
)

cat(
  "\n================================================================================\n"
)
cat("                    🎉 All analyses complete!\n")
cat(
  "================================================================================\n"
)
