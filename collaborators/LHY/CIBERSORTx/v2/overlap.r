# =============================================================================
# CIBERSORTx 重叠与相似性分析 - 增强版（含统计分析）
# 包含：层次聚类、PCA、相关性网络、统计检验等
# 图片格式：PDF
# =============================================================================

# 加载必要的包
required_packages <- c(
  "tidyverse", # 数据处理和基础绘图
  "pheatmap", # 热图
  "ggdendro", # 树状图
  "igraph", # 网络图
  "ggrepel", # 标签防重叠
  "circlize", # Chord图
  "RColorBrewer", # 颜色
  "gridExtra", # 多图排列
  "vegan", # PERMANOVA分析
  "FactoMineR", # 增强PCA分析
  "factoextra", # PCA可视化
  "officer", # PPTX生成
  "rvg", # Office DrawingML矢量图形
  "ggalluvial" # 冲积图 / Alluvial chart
)

# 安装缺失的包
new_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]
if (length(new_packages)) {
  install.packages(new_packages)
}

# 加载包
library(tidyverse)
library(pheatmap)
library(ggdendro)
library(igraph)
library(ggrepel)
library(circlize)
library(RColorBrewer)
library(gridExtra)
library(vegan)
library(FactoMineR)
library(factoextra)
library(officer)
library(rvg)
library(ggalluvial)

# =============================================================================
# 设置工作路径
# =============================================================================

setwd(
  "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\分析\\help_others\\LHY\\CIBERSORTx\\v2"
) # 修改为您的实际路径

output_dir <- "overlap_similarity_analysis_enhanced"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 组别配色方案
group_colors <- c(
  "OSF" = "#E06681",
  "Am" = "#fdd456",
  "Pm" = "#6358cb",
  "Av" = "#4fb984",
  "Ai" = "#ed9548",
  "AB" = "#ee6161",
  "B" = "#55b7ec"
)

# 将 ggplot 对象保存为 PPTX（DrawingML 矢量格式）
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

# 将 base R 图形代码保存为 PPTX（DrawingML 矢量格式）
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
cat("           CIBERSORTx 重叠与相似性分析（增强版 - 含统计检验）\n")
cat(
  "================================================================================\n\n"
)

# =============================================================================
# 1. 数据读取和预处理
# =============================================================================

cat("📊 读取数据...\n")

job2 <- read_csv("CIBERSORTx_Job2_Results.csv", show_col_types = FALSE)
job4 <- read_csv("CIBERSORTx_Job4_Results.csv", show_col_types = FALSE)

# Remove Myocytes: muscle contamination from animal sampling
all_cell_types <- colnames(job2)[2:17]
cell_types <- all_cell_types[all_cell_types != "Myocytes"]
cat(sprintf(
  "  ⚠️  Removing Myocytes (muscle contamination): %d → %d cell types\n\n",
  length(all_cell_types),
  length(cell_types)
))

job2 <- job2 |>
  select(-any_of("Myocytes")) |>
  mutate(Group = str_extract(Mixture, "^[A-Za-z]+"))

job4 <- job4 |>
  select(-any_of("Myocytes")) |>
  mutate(Group = str_extract(Mixture, "^[A-Za-z]+"))

job2_means <- job2 |>
  group_by(Group) |>
  summarise(across(all_of(cell_types), mean, na.rm = TRUE))

job4_means <- job4 |>
  group_by(Group) |>
  summarise(across(all_of(cell_types), mean, na.rm = TRUE))

# 计算变化值
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

# 合并所有变化数据
all_changes <- bind_rows(
  data.frame(Group = "OSF", t(osf_change)),
  job4_changes |> rownames_to_column("Group")
) |>
  column_to_rownames("Group")

cat("✓ 数据准备完成\n")
cat("  分析组数:", nrow(all_changes), "\n")
cat("  细胞类型数:", ncol(all_changes), "\n\n")

# =============================================================================
# 0. Cell Composition Alluvial Chart
# =============================================================================

cat("🌊 Cell Composition: Alluvial Chart...\n")

# Cell-type color palette (15 types; any unlisted type falls back to gray)
cell_type_colors <- c(
  "Epithelial" = "#e47920",
  "Fibroblasts" = "#e58f90",
  "Macrophages" = "#5e86bd",
  "Endothelial" = "#dfd882",
  "Myofibroblasts" = "#c2acca",
  "Tprolif" = "#305e80",
  "CD8+T_cells" = "#d42c7f",
  "DCs" = "#a14e92",
  "Treg" = "#e89675",
  "Plasma" = "#bcbc26",
  "CD4+T_cells" = "#2d74b2",
  "B_cells" = "#82852c",
  "Mast" = "#68843c",
  "Monocytes" = "#d0191b",
  "NK" = "#de6463"
)

# Build absolute-composition data (Normal/OSF from job2; PBS/models from job4)
alluvial_df <- bind_rows(
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

# Group order: human (Normal → OSF) → models by efficacy → animal control (PBS)
group_order_alluvial <- c(
  "Normal",
  "OSF",
  "B",
  "AB",
  "Ai",
  "Av",
  "Am",
  "Pm",
  "PBS"
)
group_order_alluvial <- group_order_alluvial[
  group_order_alluvial %in% alluvial_df$Group
]
alluvial_df <- alluvial_df |>
  mutate(Group = factor(Group, levels = group_order_alluvial))

# Color vector: map each cell type in the data to its color (gray fallback)
ct_in_data <- unique(alluvial_df$Cell_Type)
ct_colors <- cell_type_colors[ct_in_data]
ct_colors[is.na(ct_colors)] <- "#cccccc"
names(ct_colors) <- ct_in_data

p_alluvial <- ggplot(
  alluvial_df,
  aes(
    x = Group,
    y = Percentage,
    alluvium = Cell_Type,
    stratum = Cell_Type,
    fill = Cell_Type
  )
) +
  geom_alluvium(alpha = 0.65, color = NA) +
  geom_stratum(width = 0.42, color = "white", linewidth = 0.35) +
  scale_fill_manual(values = ct_colors, name = "Cell Type") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "Cell Composition Across Groups",
    x = "Group",
    y = "Percentage (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"),
    axis.text.y = element_text(face = "bold", size = 11, color = "black"),
    axis.title = element_text(face = "bold", size = 13, color = "black"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  file.path(output_dir, "0_cell_composition_alluvial.pdf"),
  p_alluvial,
  width = 10,
  height = 8,
  device = "pdf"
)
save_pptx_gg(
  p_alluvial,
  file.path(output_dir, "0_cell_composition_alluvial.pptx")
)
cat("✓ Alluvial chart saved\n\n")

# =============================================================================
# 2. 层次聚类分析 + 统计检验
# =============================================================================

cat("🌳 分析 1: 层次聚类 + 统计检验...\n")

# 计算距离矩阵
dist_matrix <- dist(all_changes, method = "euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")

# 统计检验：Cophenetic correlation（聚类质量）
coph_corr <- cor(dist_matrix, cophenetic(hc))
cat(sprintf("  Cophenetic correlation: %.3f\n", coph_corr))
cat("  (>0.75 = 优秀, >0.60 = 良好, <0.60 = 较差)\n\n")

# 绘制树状图
pdf(
  file.path(output_dir, "1_hierarchical_clustering_dendrogram.pdf"),
  width = 12,
  height = 6
)

par(mar = c(5, 5, 4, 2))
plot(
  hc,
  main = sprintf("层次聚类树状图 (Cophenetic r = %.3f)", coph_corr),
  xlab = "组别",
  ylab = "距离（Ward方法）",
  cex.main = 1.5,
  cex.lab = 1.2,
  cex.axis = 1.1,
  hang = -1
)

rect.hclust(hc, k = 3, border = "red")
text(
  x = par("usr")[2] * 0.02,
  y = par("usr")[4] * 0.95,
  labels = sprintf("聚类质量: %.3f", coph_corr),
  adj = c(0, 1),
  cex = 1.1,
  col = "blue"
)

dev.off()

save_pptx_base(
  {
    par(mar = c(5, 5, 4, 2))
    plot(
      hc,
      main = sprintf(
        "Hierarchical Clustering Dendrogram (Cophenetic r = %.3f)",
        coph_corr
      ),
      xlab = "Group",
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
      labels = sprintf("Cluster quality: %.3f", coph_corr),
      adj = c(0, 1),
      cex = 1.1,
      col = "blue"
    )
  },
  file.path(output_dir, "1_hierarchical_clustering_dendrogram.pptx")
)

cat("✓ 树状图已保存\n\n")

# =============================================================================

# =============================================================================
# 3. 主成分分析（PCA）+ 统计检验 - 超稳健版（含数据清洗）
# =============================================================================

cat("📉 分析 2: 主成分分析（PCA）+ 统计检验...\n")

# ========== 步骤1: 数据质量检查 ==========
cat("  步骤1: 检查数据质量...\n")

# 检查NA
has_na <- apply(all_changes, 2, function(x) any(is.na(x)))
if (any(has_na)) {
  cat(sprintf(
    "  ⚠️  包含NA的列: %s\n",
    paste(names(has_na)[has_na], collapse = ", ")
  ))
}

# 检查Inf
has_inf <- apply(all_changes, 2, function(x) any(is.infinite(x)))
if (any(has_inf)) {
  cat(sprintf(
    "  ⚠️  包含Inf的列: %s\n",
    paste(names(has_inf)[has_inf], collapse = ", ")
  ))
}

# 检查方差
col_vars <- apply(all_changes, 2, var, na.rm = TRUE)
zero_var <- is.na(col_vars) | col_vars == 0 | col_vars < 1e-10

if (any(zero_var)) {
  cat(sprintf(
    "  ⚠️  方差为0或过小的列: %s\n",
    paste(names(col_vars)[zero_var], collapse = ", ")
  ))
}

# ========== 步骤2: 数据清洗 ==========
cat("  步骤2: 清洗数据...\n")

# 移除有问题的列
valid_cols <- !has_na & !has_inf & !zero_var
all_changes_clean <- all_changes[, valid_cols, drop = FALSE]

cat(sprintf("  保留的细胞类型: %d / %d\n", sum(valid_cols), length(valid_cols)))

# 如果保留的列太少，报错
if (sum(valid_cols) < 2) {
  stop("可用于PCA的细胞类型少于2个，无法进行分析")
}

# ========== 步骤3: 再次检查清洗后的数据 ==========
cat("  步骤3: 验证清洗后的数据...\n")

# 确保没有NA或Inf
if (
  any(is.na(all_changes_clean)) ||
    any(is.infinite(as.matrix(all_changes_clean)))
) {
  cat("  ⚠️  数据仍包含NA或Inf，尝试替换...\n")
  # 用列均值替换NA
  for (i in 1:ncol(all_changes_clean)) {
    col_data <- all_changes_clean[, i]
    if (any(is.na(col_data))) {
      all_changes_clean[is.na(col_data), i] <- mean(col_data, na.rm = TRUE)
    }
    if (any(is.infinite(col_data))) {
      all_changes_clean[is.infinite(col_data), i] <- 0
    }
  }
}

# 打印清洗后的数据摘要
cat("  清洗后数据摘要:\n")
print(summary(all_changes_clean))

# ========== 步骤4: PCA分析 ==========
cat("\n  步骤4: 执行PCA...\n")

# 尝试多种PCA策略
pca_result <- NULL
pca_method <- ""

# 策略1: 标准化PCA
tryCatch(
  {
    pca_result <- prcomp(all_changes_clean, center = TRUE, scale. = TRUE)
    pca_method <- "标准化PCA (center=T, scale=T)"
    cat("  ✓ 使用标准化PCA\n")
  },
  error = function(e) {
    cat("  ⚠️  标准化PCA失败，尝试其他方法...\n")
  }
)

# 策略2: 只中心化
if (is.null(pca_result)) {
  tryCatch(
    {
      pca_result <- prcomp(all_changes_clean, center = TRUE, scale. = FALSE)
      pca_method <- "中心化PCA (center=T, scale=F)"
      cat("  ✓ 使用中心化PCA（未标准化）\n")
    },
    error = function(e) {
      cat("  ⚠️  中心化PCA失败，尝试原始数据...\n")
    }
  )
}

# 策略3: 原始数据
if (is.null(pca_result)) {
  tryCatch(
    {
      pca_result <- prcomp(all_changes_clean, center = FALSE, scale. = FALSE)
      pca_method <- "原始数据PCA (center=F, scale=F)"
      cat("  ✓ 使用原始数据PCA\n")
    },
    error = function(e) {
      cat("  ✗ 所有PCA方法都失败\n")
      stop("无法执行PCA分析")
    }
  )
}

# ========== 步骤5: 提取结果 ==========
pca_scores <- as.data.frame(pca_result$x[, 1:2])
pca_scores$Group <- rownames(pca_scores)

explained_var <- summary(pca_result)$importance[2, ]

cat(sprintf("\n  PC1 解释方差: %.2f%%\n", explained_var[1] * 100))
cat(sprintf("  PC2 解释方差: %.2f%%\n", explained_var[2] * 100))
cat(sprintf("  累计解释方差: %.2f%%\n", sum(explained_var[1:2]) * 100))

# ========== 步骤6: 载荷分析 ==========
if (ncol(pca_result$rotation) >= 2) {
  loadings <- pca_result$rotation[, 1:2]

  cat("\n  PC1 贡献最大的细胞类型:\n")
  top_pc1_idx <- order(abs(loadings[, 1]), decreasing = TRUE)[
    1:min(3, nrow(loadings))
  ]
  for (i in seq_along(top_pc1_idx)) {
    idx <- top_pc1_idx[i]
    cat(sprintf(
      "    %d. %s (载荷 = %.3f)\n",
      i,
      rownames(loadings)[idx],
      loadings[idx, 1]
    ))
  }

  cat("\n  PC2 贡献最大的细胞类型:\n")
  top_pc2_idx <- order(abs(loadings[, 2]), decreasing = TRUE)[
    1:min(3, nrow(loadings))
  ]
  for (i in seq_along(top_pc2_idx)) {
    idx <- top_pc2_idx[i]
    cat(sprintf(
      "    %d. %s (载荷 = %.3f)\n",
      i,
      rownames(loadings)[idx],
      loadings[idx, 2]
    ))
  }

  # 保存载荷
  loadings_df <- as.data.frame(loadings) |>
    rownames_to_column("Cell_Type") |>
    arrange(desc(abs(PC1)))
  write_csv(loadings_df, file.path(output_dir, "pca_loadings.csv"))
}

# ========== 步骤7: PERMANOVA检验 ==========
cat("\n  PERMANOVA检验（组间差异）:\n")

perm_result <- NULL
perm_r2 <- NA
perm_p <- NA

tryCatch(
  {
    group_factor <- factor(rownames(all_changes_clean))
    perm_result <- adonis2(
      all_changes_clean ~ group_factor,
      permutations = 999,
      method = "euclidean"
    )

    perm_r2 <- perm_result$R2[1]
    perm_p <- perm_result$`Pr(>F)`[1]

    cat(sprintf("  R² = %.3f\n", perm_r2))
    cat(sprintf(
      "  p-value = %.4f %s\n",
      perm_p,
      ifelse(perm_p < 0.05, "*", "")
    ))
  },
  error = function(e) {
    cat("  ⚠️  PERMANOVA检验失败\n")
  }
)

# ========== 步骤8: 保存统计结果 ==========
pca_stats <- data.frame(
  Component = paste0("PC", 1:length(explained_var)),
  Variance_Explained = explained_var * 100,
  Cumulative_Variance = cumsum(explained_var) * 100
)
write_csv(pca_stats, file.path(output_dir, "pca_variance_explained.csv"))

# ========== 步骤9: 绘制PCA图 ==========
cat("\n  绘制PCA图...\n")

plot_title <- if (!is.na(perm_p)) {
  sprintf(
    "PCA (%s)\nPERMANOVA: R² = %.3f, p = %.3f",
    pca_method,
    perm_r2,
    perm_p
  )
} else {
  sprintf("PCA (%s)", pca_method)
}

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
  scale_color_manual(values = group_colors) +
  labs(
    title = plot_title,
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

save_pptx_gg(p_pca, file.path(output_dir, "2_pca_analysis.pptx"))

cat("\n✓ PCA分析完成并保存\n\n")

# =============================================================================
# 4. 相关性矩阵 + 显著性检验
# =============================================================================

cat("🔥 分析 3: 相关性矩阵 + 显著性检验...\n")

# 计算相关性和p值
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

# 标记显著性
sig_matrix <- pval_matrix < 0.05
sig_matrix_adj <- p.adjust(pval_matrix[upper.tri(pval_matrix)], method = "BH") <
  0.05

cat(sprintf(
  "  显著相关的配对数（p < 0.05）: %d / %d\n",
  sum(sig_matrix[upper.tri(sig_matrix)]),
  sum(upper.tri(sig_matrix))
))

# 保存相关性和p值
corr_df <- as.data.frame(corr_matrix) |> rownames_to_column("Group")
write_csv(corr_df, file.path(output_dir, "correlation_matrix.csv"))

pval_df <- as.data.frame(pval_matrix) |> rownames_to_column("Group")
write_csv(pval_df, file.path(output_dir, "correlation_pvalues.csv"))

# 绘制带显著性标记的热图
pdf(
  file.path(output_dir, "3_correlation_clustered_heatmap.pdf"),
  width = 11,
  height = 10
)

# 创建显著性标记
sig_labels <- matrix("", n_groups, n_groups)
sig_labels[pval_matrix < 0.001] <- "***"
sig_labels[pval_matrix >= 0.001 & pval_matrix < 0.01] <- "**"
sig_labels[pval_matrix >= 0.01 & pval_matrix < 0.05] <- "*"

pheatmap(
  corr_matrix,
  color = colorRampPalette(c("#20a0e5", "white", "#e73838"))(100),
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
  main = "相关性矩阵聚类热图\n(* p<0.05, ** p<0.01, *** p<0.001)",
  fontsize = 11,
  fontsize_row = 11,
  fontsize_col = 11,
  border_color = "white",
  cellwidth = 55,
  cellheight = 55
)

dev.off()

save_pptx_base(
  {
    pheatmap(
      corr_matrix,
      color = colorRampPalette(c("#20a0e5", "white", "#e73838"))(100),
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
  },
  file.path(output_dir, "3_correlation_clustered_heatmap.pptx")
)

cat("✓ 相关性聚类热图已保存\n\n")

# =============================================================================
# 5. 相关性网络图
# =============================================================================

cat("🕸️  分析 4: 相关性网络图...\n")

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

cat(sprintf("  网络边数: %d (阈值 |r| > %.2f)\n", nrow(edges), threshold))
cat(sprintf("  显著的边: %d\n", sum(edges$significant)))

# 保存边数据
write_csv(edges, file.path(output_dir, "network_edges.csv"))

# 创建网络
g <- graph_from_data_frame(
  edges,
  directed = FALSE,
  vertices = rownames(corr_matrix)
)

# 绘制网络图
pdf(file.path(output_dir, "4_correlation_network.pdf"), width = 12, height = 10)

set.seed(42)
layout <- layout_with_fr(g, niter = 1000)

node_colors <- group_colors[V(g)$name]
node_sizes <- ifelse(V(g)$name == "OSF", 30, 20)

# 边的颜色和样式
edge_colors <- ifelse(E(g)$correlation > 0, "#E06681", "#8087E2")
edge_lty <- ifelse(E(g)$significant, 1, 2) # 实线=显著，虚线=不显著

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
  main = sprintf("相关性网络图（阈值 |r| > %.1f，实线=p<0.05）", threshold)
)

legend(
  "topright",
  legend = c(
    "Positive",
    "Negative",
    "Sig. (p<0.05)",
    "n.s.",
    names(group_colors)
  ),
  col = c("#E06681", "#8087E2", "black", "black", unname(group_colors)),
  lty = c(1, 1, 1, 2, rep(NA, length(group_colors))),
  pch = c(NA, NA, NA, NA, rep(21, length(group_colors))),
  pt.bg = c(NA, NA, NA, NA, unname(group_colors)),
  pt.cex = 2,
  lwd = c(3, 3, 2, 2, rep(NA, length(group_colors))),
  bty = "n",
  cex = 1
)

dev.off()

save_pptx_base(
  {
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
        "Correlation Network (|r| > %.1f; solid = p<0.05)",
        threshold
      )
    )
    legend(
      "topright",
      legend = c(
        "Positive",
        "Negative",
        "Sig. (p<0.05)",
        "n.s.",
        names(group_colors)
      ),
      col = c("#E06681", "#8087E2", "black", "black", unname(group_colors)),
      lty = c(1, 1, 1, 2, rep(NA, length(group_colors))),
      pch = c(NA, NA, NA, NA, rep(21, length(group_colors))),
      pt.bg = c(NA, NA, NA, NA, unname(group_colors)),
      pt.cex = 2,
      lwd = c(3, 3, 2, 2, rep(NA, length(group_colors))),
      bty = "n",
      cex = 1
    )
  },
  file.path(output_dir, "4_correlation_network.pptx")
)

cat("✓ 网络图已保存\n\n")

# =============================================================================
# 6. 重叠分析 + 统计检验
# =============================================================================

cat("📊 分析 5: 显著变化细胞类型的重叠分析 + 统计...\n")

sig_threshold <- 0.01

# 识别显著变化细胞
significant_cells <- list()
for (group in rownames(all_changes)) {
  sig_cells <- colnames(all_changes)[abs(all_changes[group, ]) > sig_threshold]
  significant_cells[[group]] <- sig_cells
}

cat("\n显著变化的细胞类型（|变化| > 0.01）:\n")
cat(strrep("-", 70), "\n")
for (group in names(significant_cells)) {
  cells <- significant_cells[[group]]
  cat(sprintf(
    "%-5s: %2d 种 - %s\n",
    group,
    length(cells),
    paste(sort(cells), collapse = ", ")
  ))
}
cat(strrep("-", 70), "\n\n")

# OSF与各组的重叠 + Fisher精确检验
osf_cells <- significant_cells[["OSF"]]
overlap_data <- data.frame()

cat("Fisher精确检验（重叠是否显著）:\n")
cat(strrep("-", 70), "\n")

for (group in names(significant_cells)) {
  if (group != "OSF") {
    group_cells <- significant_cells[[group]]
    overlap <- intersect(osf_cells, group_cells)

    # Fisher精确检验
    # 2x2列联表：OSF显著 vs 非显著 × 该组显著 vs 非显著
    a <- length(overlap) # 都显著
    b <- length(setdiff(osf_cells, group_cells)) # 仅OSF显著
    c <- length(setdiff(group_cells, osf_cells)) # 仅该组显著
    d <- length(cell_types) - a - b - c # 都不显著

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
        Fisher_OR_low = fisher_test$conf.int[1],
        Fisher_OR_high = fisher_test$conf.int[2],
        Fisher_pvalue = fisher_test$p.value
      )
    )
  }
}
cat(strrep("-", 70), "\n\n")

# 保存重叠统计
write_csv(
  overlap_data,
  file.path(output_dir, "overlap_analysis_with_stats.csv")
)

# =============================================================================
# Plan A: Presence / Absence heatmap
# =============================================================================

group_order <- c("OSF", setdiff(rownames(all_changes), "OSF"))

dir_df <- do.call(
  rbind,
  lapply(rownames(all_changes), function(g) {
    data.frame(
      Group = g,
      Cell_Type = cell_types,
      Direction = sign(as.numeric(all_changes[g, cell_types])) *
        (abs(as.numeric(all_changes[g, cell_types])) > sig_threshold)
    )
  })
) |>
  mutate(
    Direction_label = factor(
      case_when(Direction == 1 ~ "Up", Direction == -1 ~ "Down", TRUE ~ "n.s."),
      levels = c("Up", "n.s.", "Down")
    )
  )

p_presence <- ggplot(
  dir_df,
  aes(
    x = factor(Group, levels = group_order),
    y = Cell_Type,
    fill = Direction_label
  )
) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_vline(xintercept = 1.5, color = "black", linewidth = 1.2) +
  scale_fill_manual(
    values = c("Up" = "#e73838", "n.s." = "#eeeeee", "Down" = "#20a0e5"),
    name = "Change direction"
  ) +
  scale_x_discrete(limits = group_order) +
  labs(
    title = "Cell Type Change Profiles Across Groups",
    subtitle = sprintf(
      "Threshold: |Δ| > %.2f   |   OSF = human disease   |   others = animal models",
      sig_threshold
    ),
    x = "Group",
    y = "Cell type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    panel.grid = element_blank()
  )

ggsave(
  file.path(output_dir, "5a_overlap_heatmap.pdf"),
  p_presence,
  width = 9,
  height = 7,
  device = "pdf"
)
save_pptx_gg(p_presence, file.path(output_dir, "5a_overlap_heatmap.pptx"))
cat("✓ Plan A: presence/absence heatmap saved\n")

# =============================================================================
# Plan D: Lollipop (overlap rate) + Fisher OR with 95% CI
# =============================================================================

overlap_d <- overlap_data |>
  arrange(Overlap_Rate_OSF) |>
  mutate(
    Group = factor(Group, levels = Group),
    sig_fill = Fisher_pvalue < 0.05,
    sig_label = case_when(
      Fisher_pvalue < 0.001 ~ "***",
      Fisher_pvalue < 0.01 ~ "**",
      Fisher_pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # cap Inf CI upper bounds for log-scale display
    Fisher_OR_high_plot = pmin(Fisher_OR_high, 500)
  )

p_lollipop <- ggplot(overlap_d, aes(y = Group, x = Overlap_Rate_OSF * 100)) +
  geom_segment(
    aes(xend = 0, yend = Group, color = sig_fill),
    linewidth = 0.8
  ) +
  geom_point(aes(size = Overlap, color = sig_fill)) +
  geom_text(
    aes(x = Overlap_Rate_OSF * 100, label = sig_label),
    nudge_x = 2,
    size = 5,
    hjust = 0,
    fontface = "bold"
  ) +
  scale_color_manual(
    values = c("TRUE" = "#9467BD", "FALSE" = "#AAAAAA"),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
    name = "Fisher exact test"
  ) +
  scale_size_continuous(name = "Shared\ncell types", range = c(3, 8)) +
  labs(title = "Overlap Rate with OSF", x = "Overlap rate (%)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

p_or <- ggplot(overlap_d, aes(y = Group, x = Fisher_OR, color = sig_fill)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = pmax(Fisher_OR_low, 0.01), xmax = Fisher_OR_high_plot),
    height = 0.2,
    linewidth = 0.8
  ) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("TRUE" = "#9467BD", "FALSE" = "#AAAAAA"),
    guide = "none"
  ) +
  scale_x_log10() +
  labs(title = "Fisher OR (95% CI)", x = "Odds ratio (log scale)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank()
  )

p_plan_d <- gridExtra::arrangeGrob(
  p_lollipop,
  p_or,
  ncol = 2,
  widths = c(2, 1.5)
)

ggsave(
  file.path(output_dir, "5b_overlap_lollipop.pdf"),
  p_plan_d,
  width = 12,
  height = 6
)
save_pptx_base(
  {
    grid::grid.draw(p_plan_d)
  },
  file.path(output_dir, "5b_overlap_lollipop.pptx")
)
cat("✓ Plan D: lollipop + OR plot saved\n\n")

# =============================================================================
# Plan C: Faceted donut chart — overlap composition per model group
# =============================================================================

overlap_pie_df <- overlap_data |>
  pivot_longer(
    cols = c(OSF_Only, Overlap, Group_Only),
    names_to = "Category",
    values_to = "Count"
  ) |>
  mutate(
    Category = factor(
      Category,
      levels = c("OSF_Only", "Overlap", "Group_Only"),
      labels = c("OSF only", "Shared", "Model only")
    ),
    sig_label = sprintf("p = %.3f", Fisher_pvalue),
    Group_label = paste0(Group, "\n(", sig_label, ")")
  )

p_pie <- ggplot(
  overlap_pie_df,
  aes(x = 2, y = Count, fill = Category)
) +
  geom_bar(
    stat = "identity",
    width = 1,
    color = "white",
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = ifelse(Count > 0, Count, "")),
    position = position_stack(vjust = 0.5),
    size = 3.8,
    color = "white",
    fontface = "bold"
  ) +
  coord_polar(theta = "y", start = 0) +
  xlim(c(0.5, 2.8)) +
  facet_wrap(~Group_label, nrow = 2) +
  scale_fill_manual(
    values = c(
      "OSF only" = "#E06681",
      "Shared" = "#9467BD",
      "Model only" = "#8087E2"
    ),
    name = "Cell type category"
  ) +
  labs(
    title = "Overlap Composition per Model Group",
    subtitle = sprintf("Threshold: |Δ| > %.2f", sig_threshold)
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom"
  )

ggsave(
  file.path(output_dir, "5c_overlap_pie.pdf"),
  p_pie,
  width = 12,
  height = 7,
  device = "pdf"
)
save_pptx_gg(p_pie, file.path(output_dir, "5c_overlap_pie.pptx"))
cat("✓ Plan C: donut chart saved\n\n")

# =============================================================================
# 7. Chord图
# =============================================================================

cat("🎵 分析 6: Chord图...\n")

chord_threshold <- 0.2
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
  grid_colors <- group_colors[groups]

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
    sprintf("Chord Diagram: Group Correlations (|r| > %.1f)", chord_threshold),
    cex.main = 1.5
  )

  dev.off()
  circos.clear()

  save_pptx_base(
    {
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
          xlim <- get.cell.meta.data("xlim")
          ylim <- get.cell.meta.data("ylim")
          sector.name <- get.cell.meta.data("sector.index")
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
          "Chord Diagram: Group Correlations (|r| > %.1f)",
          chord_threshold
        ),
        cex.main = 1.5
      )
      circos.clear()
    },
    file.path(output_dir, "6_chord_diagram.pptx")
  )

  cat("✓ Chord图已保存\n\n")
} else {
  cat("⚠️  没有满足阈值的相关性，跳过Chord图\n\n")
}


# =============================================================================
# 8. 生成综合统计报告
# =============================================================================

cat("📝 生成统计分析总结报告...\n\n")

# 安全获取PCA载荷信息
pc1_contrib_text <- if (
  exists("top_pc1") && length(top_pc1) > 0 && exists("loadings")
) {
  paste(
    sprintf(
      "    %d. %s (载荷 = %.3f)",
      1:length(top_pc1),
      top_pc1,
      loadings[top_pc1, 1]
    ),
    collapse = "\n"
  )
} else {
  "    数据不可用"
}

pc2_contrib_text <- if (
  exists("top_pc2") && length(top_pc2) > 0 && exists("loadings")
) {
  paste(
    sprintf(
      "    %d. %s (载荷 = %.3f)",
      1:length(top_pc2),
      top_pc2,
      loadings[top_pc2, 2]
    ),
    collapse = "\n"
  )
} else {
  "    数据不可用"
}

# 安全获取PERMANOVA结果
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
  "    未执行或失败"
}

# 安全获取相关性统计
total_pairs <- sum(upper.tri(corr_matrix))
sig_pairs <- if (exists("pval_matrix")) {
  sum(pval_matrix[upper.tri(pval_matrix)] < 0.05, na.rm = TRUE)
} else {
  0
}
pos_corr <- sum(corr_matrix[upper.tri(corr_matrix)] > 0, na.rm = TRUE)
neg_corr <- sum(corr_matrix[upper.tri(corr_matrix)] < 0, na.rm = TRUE)

# 安全获取重叠统计
fisher_text <- if (exists("overlap_data") && nrow(overlap_data) > 0) {
  paste(
    sprintf(
      "    %s: 重叠率=%.1f%%, OR=%.2f, p=%.4f %s",
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
  "    数据不可用"
}

# 安全获取网络统计
network_stats <- if (exists("edges") && nrow(edges) > 0) {
  sprintf(
    "  网络阈值: |r| > %.2f\n  总边数: %d\n  显著边数: %d (p<0.05)\n  平均相关强度: %.3f",
    threshold,
    nrow(edges),
    sum(edges$significant, na.rm = TRUE),
    mean(edges$weight, na.rm = TRUE)
  )
} else {
  "  数据不可用"
}

# 生成报告
summary_text <- sprintf(
  "
================================================================================
              重叠与相似性分析 - 统计分析总结报告
================================================================================

【1. 层次聚类统计】
  聚类方法: Ward's D2
  距离度量: Euclidean
  Cophenetic correlation: %.3f
  解读: %.3f %s
  
【2. 主成分分析（PCA）统计】
  PC1 解释方差: %.2f%%
  PC2 解释方差: %.2f%%
  累计解释方差: %.2f%%
  
  PC1主要贡献细胞（前3）:
%s
  
  PC2主要贡献细胞（前3）:
%s
  
  PERMANOVA检验:
%s
    
【3. 相关性分析统计】
  总配对数: %d
  显著相关对（p<0.05）: %d (%.1f%%)
  正相关: %d
  负相关: %d
  
【4. 重叠分析统计（Fisher精确检验）】
  OSF显著变化细胞数: %d
  
  各组与OSF重叠的显著性:
%s

【5. 网络分析统计】
%s
  
【统计解释】
  * p < 0.05: 显著
  ** p < 0.01: 高度显著
  *** p < 0.001: 极显著
  
  Cophenetic correlation:
    > 0.75: 聚类结构优秀
    0.60-0.75: 聚类结构良好
    < 0.60: 聚类结构较差
    
  PERMANOVA:
    检验组间细胞变化模式是否显著不同
    R²: 组间差异占总变异的比例
    
  Fisher精确检验:
    检验两组显著变化细胞的重叠是否超过随机期望
    OR > 1: 重叠高于随机期望
    OR < 1: 重叠低于随机期望

【输出文件】
  Visualizations (PDF + PPTX vector):
    • 1_hierarchical_clustering_dendrogram
    • 2_pca_analysis
    • 3_correlation_clustered_heatmap
    • 4_correlation_network
    • 5a_overlap_heatmap          (Plan A: presence/absence heatmap)
    • 5b_overlap_lollipop         (Plan D: lollipop + Fisher OR)
    • 5c_overlap_pie              (donut chart: overlap composition per group)
    • 6_chord_diagram
  
  统计数据表格:
    • correlation_matrix.csv - 相关系数矩阵
    • correlation_pvalues.csv - 相关性p值矩阵
    • pca_variance_explained.csv - PCA方差解释
    • pca_loadings.csv - PCA载荷（细胞贡献）
    • overlap_analysis_with_stats.csv - 重叠分析+Fisher检验
    • network_edges.csv - 网络边+统计信息

================================================================================
✅ 统计分析完成！
================================================================================
",
  # 层次聚类
  coph_corr,
  coph_corr,
  ifelse(
    coph_corr > 0.75,
    "（优秀）",
    ifelse(coph_corr > 0.60, "（良好）", "（较差）")
  ),

  # PCA
  explained_var[1] * 100,
  explained_var[2] * 100,
  sum(explained_var[1:2]) * 100,

  # PC贡献
  pc1_contrib_text,
  pc2_contrib_text,

  # PERMANOVA
  permanova_text,

  # 相关性
  total_pairs,
  sig_pairs,
  ifelse(total_pairs > 0, sig_pairs / total_pairs * 100, 0),
  pos_corr,
  neg_corr,

  # 重叠
  length(osf_cells),
  fisher_text,

  # 网络
  network_stats
)

cat(summary_text)

writeLines(
  summary_text,
  file.path(output_dir, "STATISTICAL_ANALYSIS_SUMMARY.txt")
)

cat("\n✓ 统计分析总结已保存\n")
cat(
  "\n================================================================================\n"
)
cat("                     🎉 所有分析完成（含统计检验）！\n")
cat(
  "================================================================================\n"
)
