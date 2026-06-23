# =============================================================================
# CIBERSORTx 细胞丰度相似性分析 - 方法2：绝对值比较
# 直接比较OSF疾病组与各建模组的细胞组成相似性
# =============================================================================

library(tidyverse)
library(tidyplots)
library(pheatmap)

# =============================================================================
# ⚙️  设置工作路径和输出目录
# =============================================================================

# 设置数据文件路径
setwd(
  "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\分析\\help_others\\LHY\\CIBERSORTx\\v2"
) # 修改为您的实际路径

# 创建输出目录
output_dir <- "method2_absolute_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# =============================================================================
# 1. 数据读取
# =============================================================================

cat("📊 读取数据...\n")

job2 <- read_csv("CIBERSORTx_Job2_Results.csv", show_col_types = FALSE)
job4 <- read_csv("CIBERSORTx_Job4_Results.csv", show_col_types = FALSE)

cell_types <- colnames(job2)[2:17]

cat("✓ 数据读取完成\n")
cat("  Job2 样本数:", nrow(job2), "\n")
cat("  Job4 样本数:", nrow(job4), "\n")
cat("  细胞类型数:", length(cell_types), "\n\n")

# =============================================================================
# 2. 数据预处理
# =============================================================================

cat("🔍 提取分组信息...\n")

job2_processed <- job2 |>
  mutate(Group = str_extract(Mixture, "^[A-Za-z]+")) |>
  select(Mixture, Group, all_of(cell_types))

job4_processed <- job4 |>
  mutate(Group = str_extract(Mixture, "^[A-Za-z]+")) |>
  select(Mixture, Group, all_of(cell_types))

job2_means <- job2_processed |>
  group_by(Group) |>
  summarise(across(all_of(cell_types), mean, na.rm = TRUE))

job4_means <- job4_processed |>
  group_by(Group) |>
  summarise(across(all_of(cell_types), mean, na.rm = TRUE))

cat("  Job2 分组:", paste(job2_means$Group, collapse = ", "), "\n")
cat("  Job4 分组:", paste(job4_means$Group, collapse = ", "), "\n\n")

# =============================================================================
# 3. 对照组一致性检查（批次效应评估）
# =============================================================================

cat("⚠️  检查对照组一致性（批次效应评估）...\n\n")

normal_abundance <- job2_means |>
  filter(Group == "Normal") |>
  select(-Group) |>
  as.numeric()
pbs_abundance <- job4_means |>
  filter(Group == "PBS") |>
  select(-Group) |>
  as.numeric()

control_comparison <- data.frame(
  Cell_Type = cell_types,
  Normal_Job2 = normal_abundance,
  PBS_Job4 = pbs_abundance,
  Difference = normal_abundance - pbs_abundance,
  Percent_Diff = round(
    (normal_abundance - pbs_abundance) / normal_abundance * 100,
    2
  )
)

# 计算对照组相似性
control_corr <- cor(normal_abundance, pbs_abundance, method = "pearson")
control_euclidean <- sqrt(sum((normal_abundance - pbs_abundance)^2))

cat("对照组对比（Normal vs PBS）:\n")
cat(strrep("-", 70), "\n")
print(control_comparison, row.names = FALSE)
cat(strrep("-", 70), "\n\n")

cat(sprintf("对照组相关性: r = %.3f\n", control_corr))
cat(sprintf("对照组欧氏距离: d = %.4f\n", control_euclidean))

if (control_corr > 0.8) {
  cat("✓ 对照组高度一致，批次效应较小，方法2结果可靠\n")
} else if (control_corr > 0.6) {
  cat("⚠ 对照组中度一致，存在一定批次效应，需谨慎解读\n")
} else {
  cat("✗ 对照组差异较大，批次效应明显，建议优先使用方法1（变化值）\n")
}

cat("\n")

# 保存对照组对比
write_csv(
  control_comparison,
  file.path(output_dir, "control_group_comparison.csv")
)

# =============================================================================
# 4. 方法2：绝对细胞丰度相似性分析
# =============================================================================

cat("🔬 方法2：计算OSF与各建模组的绝对细胞丰度相似性...\n\n")

# 提取OSF的绝对丰度
osf_absolute <- job2_means |>
  filter(Group == "OSF") |>
  select(-Group) |>
  as.numeric()
names(osf_absolute) <- cell_types

# Job4各组的绝对丰度（不包含PBS）
job4_absolute <- job4_means |>
  filter(Group != "PBS") |>
  column_to_rownames("Group")

# 计算相似性指标
correlations_abs <- sapply(rownames(job4_absolute), function(group) {
  cor(osf_absolute, as.numeric(job4_absolute[group, ]), method = "pearson")
})

euclidean_dist_abs <- sapply(rownames(job4_absolute), function(group) {
  sqrt(sum((osf_absolute - as.numeric(job4_absolute[group, ]))^2))
})

cosine_sim_abs <- sapply(rownames(job4_absolute), function(group) {
  sum(osf_absolute * as.numeric(job4_absolute[group, ])) /
    (sqrt(sum(osf_absolute^2)) *
      sqrt(sum(as.numeric(job4_absolute[group, ])^2)))
})

# 汇总结果
similarity_absolute <- data.frame(
  Modeling_Group = names(correlations_abs),
  Pearson_r = correlations_abs,
  Euclidean_Distance = euclidean_dist_abs,
  Cosine_Similarity = cosine_sim_abs
) |>
  arrange(desc(Pearson_r)) |>
  mutate(Rank = row_number())

cat("═══════════════════════════════════════════════════════════\n")
cat("📊 各建模方法与真实OSF疾病的细胞组成相似性排名\n")
cat("═══════════════════════════════════════════════════════════\n")
print(similarity_absolute, row.names = FALSE)
cat("═══════════════════════════════════════════════════════════\n\n")

# 保存结果
write_csv(
  similarity_absolute,
  file.path(output_dir, "similarity_metrics_absolute.csv")
)

# 识别最佳和最差建模组
best_model <- similarity_absolute$Modeling_Group[1]
worst_model <- similarity_absolute$Modeling_Group[nrow(similarity_absolute)]

cat(sprintf("🏆 最佳建模组: %s\n", best_model))
cat(sprintf("   相关系数: r = %.3f\n", similarity_absolute$Pearson_r[1]))
cat(sprintf(
  "   欧氏距离: d = %.4f\n\n",
  similarity_absolute$Euclidean_Distance[1]
))

cat(sprintf("❌ 最差建模组: %s\n", worst_model))
cat(sprintf(
  "   相关系数: r = %.3f\n",
  similarity_absolute$Pearson_r[nrow(similarity_absolute)]
))
cat(sprintf(
  "   欧氏距离: d = %.4f\n\n",
  similarity_absolute$Euclidean_Distance[nrow(similarity_absolute)]
))

# =============================================================================
# 5. 可视化 1: 绝对细胞丰度热图
# =============================================================================

cat("📊 生成可视化 1: 绝对细胞丰度热图...\n")

# 准备热图数据
heatmap_data <- rbind(
  data.frame(Group = "OSF_Disease", t(osf_absolute)),
  data.frame(Group = rownames(job4_absolute), job4_absolute)
) |>
  pivot_longer(
    cols = -Group,
    names_to = "Cell_Type",
    values_to = "Abundance"
  ) |>
  mutate(
    Group = factor(Group, levels = c("OSF_Disease", rownames(job4_absolute))),
    Cell_Type = factor(Cell_Type, levels = cell_types)
  )

p1 <- heatmap_data |>
  tidyplot(x = Cell_Type, y = Group, color = Abundance) |>
  add_heatmap() |>
  adjust_colors(colors_continuous_viridis) |>
  adjust_title("绝对细胞丰度热图：OSF真实疾病 vs 各建模方法") |>
  adjust_x_axis_title("细胞类型") |>
  adjust_y_axis_title("") |>
  adjust_legend_title("丰度") |>
  adjust_size(width = 150, height = 80)

print(p1)
ggsave(
  file.path(output_dir, "1_absolute_abundance_heatmap.pdf"),
  p1,
  width = 10,
  height = 5
)

cat("✓ 热图已保存\n\n")

# =============================================================================
# 6. 可视化 2: 建模质量排名
# =============================================================================

cat("📊 生成可视化 2: 建模质量排名...\n")

p2 <- similarity_absolute |>
  mutate(Modeling_Group = factor(Modeling_Group, levels = Modeling_Group)) |>
  tidyplot(x = Modeling_Group, y = Pearson_r, color = Pearson_r) |>
  add_mean_bar() |>
  adjust_colors(colors_diverging_blue2red) |>
  adjust_title("建模质量评估：各组与OSF疾病的相似性") |>
  adjust_x_axis_title("建模组") |>
  adjust_y_axis_title("Pearson 相关系数") |>
  remove_legend() |>
  adjust_size(width = 100, height = 60)

print(p2)
ggsave(
  file.path(output_dir, "2_modeling_quality_ranking.pdf"),
  p2,
  width = 8,
  height = 5
)

cat("✓ 排名图已保存\n\n")

# =============================================================================
# 7. 可视化 3: 最佳建模组 vs OSF
# =============================================================================

cat("📊 生成可视化 3: 最佳建模组 vs OSF 散点图...\n")

scatter_best <- data.frame(
  Cell_Type = cell_types,
  OSF = osf_absolute,
  Best_Model = as.numeric(job4_absolute[best_model, ])
)

p3 <- scatter_best |>
  tidyplot(x = OSF, y = Best_Model) |>
  add_data_points(size = 3, alpha = 0.7, color = "#8087E2") |>
  add_curve_fit(method = "lm", color = "#E06681", linewidth = 1.2) |>
  add_data_labels_repel(label = Cell_Type, size = 1.5) |>
  adjust_title(paste0(
    "最佳建模组 ",
    best_model,
    " vs OSF真实疾病\nr = ",
    round(correlations_abs[best_model], 3)
  )) |>
  adjust_x_axis_title("OSF疾病组细胞丰度") |>
  adjust_y_axis_title(paste0(best_model, " 建模组细胞丰度")) |>
  adjust_size(width = 100, height = 80)

print(p3)
ggsave(
  file.path(output_dir, "3_best_model_vs_osf.pdf"),
  p3,
  width = 8,
  height = 6
)

cat("✓ 散点图已保存\n\n")

# =============================================================================
# 8. 可视化 4: 最佳 vs 最差建模组对比
# =============================================================================

cat("📊 生成可视化 4: 最佳 vs 最差建模组对比...\n")

comparison_data <- data.frame(
  Cell_Type = cell_types,
  OSF = osf_absolute,
  Best = as.numeric(job4_absolute[best_model, ]),
  Worst = as.numeric(job4_absolute[worst_model, ])
) |>
  pivot_longer(
    cols = c(OSF, Best, Worst),
    names_to = "Group",
    values_to = "Abundance"
  ) |>
  mutate(
    Cell_Type = factor(Cell_Type, levels = cell_types),
    Group = factor(Group, levels = c("OSF", "Best", "Worst"))
  )

p4 <- comparison_data |>
  tidyplot(x = Cell_Type, y = Abundance, color = Group) |>
  add_mean_bar(alpha = 0.7) |>
  adjust_colors(c("#E06681", "#8087E2", "#CCCCCC")) |>
  adjust_title("最佳与最差建模组对比") |>
  adjust_x_axis_title("细胞类型") |>
  adjust_y_axis_title("细胞丰度") |>
  adjust_legend_title("") |>
  adjust_size(width = 120, height = 70)

print(p4)
ggsave(
  file.path(output_dir, "4_best_vs_worst_comparison.pdf"),
  p4,
  width = 10,
  height = 5
)

cat("✓ 对比图已保存\n\n")

# =============================================================================
# 9. 生成分析总结
# =============================================================================

cat("📝 生成分析总结...\n\n")

cat(
  "═══════════════════════════════════════════════════════════════════════════\n"
)
cat(
  "                  方法2：绝对细胞丰度相似性分析总结                          \n"
)
cat(
  "═══════════════════════════════════════════════════════════════════════════\n\n"
)

cat("【批次效应评估】\n")
cat(sprintf("  对照组相关性（Normal vs PBS）: r = %.3f\n", control_corr))
cat(sprintf("  对照组欧氏距离: d = %.4f\n", control_euclidean))
if (control_corr < 0.6) {
  cat("  ⚠️  警告: 批次效应明显，建议优先使用方法1（变化值比较）\n")
}
cat("\n")

cat("【建模质量排名】\n")
for (i in 1:nrow(similarity_absolute)) {
  cat(sprintf(
    "  %d. %-5s  r = %.3f  距离 = %.4f\n",
    i,
    similarity_absolute$Modeling_Group[i],
    similarity_absolute$Pearson_r[i],
    similarity_absolute$Euclidean_Distance[i]
  ))
}

cat("\n【关键发现】\n")
cat(sprintf(
  "  🏆 最佳建模组: %s (r = %.3f)\n",
  best_model,
  similarity_absolute$Pearson_r[1]
))
cat(sprintf(
  "  ❌ 最差建模组: %s (r = %.3f)\n",
  worst_model,
  similarity_absolute$Pearson_r[nrow(similarity_absolute)]
))

cat("\n【输出文件】\n")
cat("  所有结果已保存至:", output_dir, "\n")
cat("  • 相似性指标: similarity_metrics_absolute.csv\n")
cat("  • 对照组对比: control_group_comparison.csv\n")
cat("  • 可视化图表: 4张PDF文件\n")

cat(
  "\n═══════════════════════════════════════════════════════════════════════════\n"
)
cat("✅ 方法2分析完成！\n")
cat(
  "═══════════════════════════════════════════════════════════════════════════\n"
)
