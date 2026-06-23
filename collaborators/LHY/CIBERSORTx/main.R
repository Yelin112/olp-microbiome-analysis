# =============================================================================
# CIBERSORTx 细胞丰度变化相似性分析
# 比较Job2(OSF vs Normal)和Job4(各试验组 vs PBS)的细胞丰度变化模式
# =============================================================================

# 安装必要的包（如果尚未安装）
required_packages <- c("tidyverse", "tidyplots", "corrplot", "pheatmap")
new_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]
if (length(new_packages)) {
  install.packages(new_packages)
}

# 加载包
library(tidyverse)
library(tidyplots)
library(corrplot)
library(pheatmap)

# =============================================================================
# ⚙️  设置工作路径（请根据您的实际路径修改）
# =============================================================================

# 方法1: 设置工作目录到文件所在文件夹
setwd(
  "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\分析\\help_others\\LHY\\CIBERSORTx"
) # 修改为您的实际路径

# 方法2: 或者直接使用完整路径
# data_path <- "//help_others//LHY//CIBERSORTx//" # 修改为您的实际路径

# =============================================================================
# 1. 数据读取
# =============================================================================

cat("📊 读取数据...\n")

# 读取Job2数据 (OSF vs Normal)
job2 <- read_csv("CIBERSORTx_Job2_Results.csv", show_col_types = FALSE)

# 读取Job4数据 (各试验组 vs PBS)
job4 <- read_csv("CIBERSORTx_Job4_Results.csv", show_col_types = FALSE)

# 定义细胞类型列（排除样本名和统计列）
cell_types <- colnames(job2)[2:17]

cat("✓ 数据读取完成\n")
cat("  Job2 样本数:", nrow(job2), "\n")
cat("  Job4 样本数:", nrow(job4), "\n")
cat("  细胞类型数:", length(cell_types), "\n\n")

# =============================================================================
# 2. 数据预处理 - 提取分组信息
# =============================================================================

cat("🔍 提取分组信息...\n")

# Job2: 提取组别 (OSF vs Normal)
job2_processed <- job2 |>
  mutate(Group = str_extract(Mixture, "^[A-Za-z]+")) |>
  select(Mixture, Group, all_of(cell_types))

# Job4: 提取组别 (B, Av, AB, Am, Pm, Ai, PBS)
job4_processed <- job4 |>
  mutate(Group = str_extract(Mixture, "^[A-Za-z]+")) |>
  select(Mixture, Group, all_of(cell_types))

# 统计各组样本数
job2_summary <- table(job2_processed$Group)
job4_summary <- table(job4_processed$Group)

cat("  Job2 分组:", paste(names(job2_summary), collapse = ", "), "\n")
cat("  Job4 分组:", paste(names(job4_summary), collapse = ", "), "\n\n")

# =============================================================================
# 3. 计算各组的平均细胞丰度
# =============================================================================

cat("📈 计算各组平均细胞丰度...\n")

# Job2: 计算OSF和Normal组的均值
job2_means <- job2_processed |>
  group_by(Group) |>
  summarise(across(all_of(cell_types), mean, na.rm = TRUE))

# Job4: 计算各试验组的均值
job4_means <- job4_processed |>
  group_by(Group) |>
  summarise(across(all_of(cell_types), mean, na.rm = TRUE))

cat("✓ 均值计算完成\n\n")

# =============================================================================
# 4. 计算细胞丰度变化值 (Δ)
# =============================================================================

cat("📊 计算细胞丰度变化...\n")

# Job2: OSF - Normal
osf_change <- (job2_means |>
  filter(Group == "OSF") |>
  select(-Group) |>
  as.numeric()) -
  (job2_means |> filter(Group == "Normal") |> select(-Group) |> as.numeric())
names(osf_change) <- cell_types

# Job4: 各试验组 - PBS
pbs_mean <- job4_means |>
  filter(Group == "PBS") |>
  select(-Group) |>
  as.numeric()


job4_changes <- job4_means |>
  filter(Group != "PBS") |>
  column_to_rownames("Group") |>
  sweep(2, pbs_mean, "-")
cat("✓ 变化值计算完成\n")
cat("  Job2 OSF变化向量:", length(osf_change), "个细胞类型\n")
cat("  Job4 试验组数:", nrow(job4_changes), "\n\n")

# =============================================================================
# 5. 相关性分析 - 核心分析
# =============================================================================

cat("🔬 计算相关性矩阵...\n")

# 计算每个Job4试验组与OSF的相关性
correlations <- sapply(rownames(job4_changes), function(group) {
  cor(osf_change, as.numeric(job4_changes[group, ]), method = "pearson")
})

# 按相关性排序
correlations_sorted <- sort(correlations, decreasing = TRUE)

cat("✓ 相关性分析完成\n\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("📊 OSF疾病模型与各试验组的相似性排名 (Pearson相关系数)\n")
cat("═══════════════════════════════════════════════════════════\n")
for (i in 1:length(correlations_sorted)) {
  cat(sprintf(
    "  %d. %-5s  r = %6.3f  %s\n",
    i,
    names(correlations_sorted)[i],
    correlations_sorted[i],
    ifelse(correlations_sorted[i] > 0, "↑ 正相关 (加剧)", "↓ 负相关 (逆转)")
  ))
}
cat("═══════════════════════════════════════════════════════════\n\n")

# =============================================================================
# 6. 欧氏距离和余弦相似度
# =============================================================================

cat("📏 计算距离和相似度指标...\n")

# 欧氏距离
euclidean_dist <- sapply(rownames(job4_changes), function(group) {
  sqrt(sum((osf_change - as.numeric(job4_changes[group, ]))^2))
})

# 余弦相似度
cosine_sim <- sapply(rownames(job4_changes), function(group) {
  sum(osf_change * as.numeric(job4_changes[group, ])) /
    (sqrt(sum(osf_change^2)) * sqrt(sum(as.numeric(job4_changes[group, ])^2)))
})

# 综合指标汇总
similarity_metrics <- data.frame(
  Trial_Group = names(correlations),
  Pearson_r = correlations,
  Euclidean_Distance = euclidean_dist,
  Cosine_Similarity = cosine_sim
) |>
  arrange(desc(Pearson_r))

cat("✓ 距离指标计算完成\n\n")
print(similarity_metrics, row.names = FALSE)
cat("\n")

# =============================================================================
# 7. 识别关键细胞类型
# =============================================================================

cat("🎯 识别OSF中显著变化的细胞类型...\n")

# 计算变化幅度（绝对值）
osf_change_abs <- abs(osf_change)
osf_change_ranked <- sort(osf_change_abs, decreasing = TRUE)

cat("\nTop 5 在OSF中变化最大的细胞类型:\n")
cat("─────────────────────────────────────────\n")
for (i in 1:5) {
  cell <- names(osf_change_ranked)[i]
  change <- osf_change[cell]
  cat(sprintf("  %d. %-20s  Δ = %+.4f\n", i, cell, change))
}
cat("─────────────────────────────────────────\n\n")

# =============================================================================
# 8. 准备可视化数据
# =============================================================================

cat("🎨 准备可视化数据...\n")

# 创建长格式数据用于热图
all_changes <- rbind(
  data.frame(Group = "OSF", t(osf_change)),
  data.frame(Group = rownames(job4_changes), job4_changes)
)

# 转换为长格式
changes_long <- all_changes |>
  pivot_longer(cols = -Group, names_to = "Cell_Type", values_to = "Change") |>
  mutate(
    Group = factor(Group, levels = c("OSF", rownames(job4_changes))),
    Cell_Type = factor(Cell_Type, levels = cell_types)
  )

cat("✓ 可视化数据准备完成\n\n")

# =============================================================================
# 9. 可视化 1: 细胞丰度变化热图（tidyplots版本）
# =============================================================================

cat("📊 生成热图 1: 细胞丰度变化全景图...\n")

p1 <- changes_long |>
  tidyplot(x = Cell_Type, y = Group, color = Change) |>
  add_heatmap() |>
  adjust_colors(colors_diverging_blue2red) |>
  adjust_title("细胞丰度变化热图 (相对于对照组)") |>
  adjust_x_axis_title("细胞类型") |>
  adjust_y_axis_title("试验组") |>
  adjust_legend_title("丰度变化 (Δ)") |>
  adjust_size(width = 150, height = 80) |>
  theme_tidyplot()

print(p1)

# 保存
ggsave("heatmap_all_changes.pdf", p1, width = 10, height = 5)
cat("✓ 热图已保存: heatmap_all_changes.pdf\n\n")

# =============================================================================
# 10. 可视化 2: 相关性矩阵柱状图
# =============================================================================

cat("📊 生成热图 2: 相关性矩阵...\n")

corr_data <- data.frame(
  Trial_Group = names(correlations),
  OSF_Correlation = correlations
) |>
  arrange(desc(OSF_Correlation))

p2 <- corr_data |>
  mutate(Trial_Group = factor(Trial_Group, levels = Trial_Group)) |>
  tidyplot(x = Trial_Group, y = OSF_Correlation, color = OSF_Correlation) |>
  add_mean_bar() |>
  adjust_colors(colors_diverging_blue2red) |>
  adjust_title("各试验组与OSF疾病模型的相似性") |>
  adjust_x_axis_title("试验组") |>
  adjust_y_axis_title("Pearson 相关系数") |>
  adjust_y_axis(limits = c(-1, 1)) |>
  remove_legend() |>
  adjust_size(width = 100, height = 60) |>
  theme_tidyplot()

print(p2)

ggsave("correlation_barplot.pdf", p2, width = 8, height = 5)
cat("✓ 相关性柱状图已保存: correlation_barplot.pdf\n\n")

# =============================================================================
# 11. 可视化 3: Top细胞类型对比
# =============================================================================

cat("📊 生成可视化 3: 关键细胞类型变化对比...\n")

# 选择在OSF中变化最大的前5个细胞类型
top5_cells <- names(osf_change_ranked)[1:5]

top_changes_data <- changes_long |>
  filter(Cell_Type %in% top5_cells)

p3 <- top_changes_data |>
  tidyplot(x = Group, y = Change, color = Cell_Type) |>
  add_mean_bar(alpha = 0.7) |>
  adjust_colors(colors_discrete_friendly) |>
  adjust_title("OSF中变化最大的5种细胞类型在各组的表现") |>
  adjust_x_axis_title("试验组") |>
  adjust_y_axis_title("丰度变化 (Δ)") |>
  adjust_legend_title("细胞类型") |>
  adjust_size(width = 120, height = 70) |>
  theme_tidyplot()

print(p3)

ggsave("top5_cells_comparison.pdf", p3, width = 10, height = 5)
cat("✓ 关键细胞对比图已保存: top5_cells_comparison.pdf\n\n")

# =============================================================================
# 12. 可视化 4: OSF vs 最相似试验组对比
# =============================================================================

cat("📊 生成可视化 4: OSF vs 最相似试验组对比...\n")

# 找出相关性最高的试验组
top_group <- names(correlations_sorted)[1]

# 准备对比数据
radar_data <- data.frame(
  Cell_Type = cell_types,
  OSF = osf_change,
  Best_Match = as.numeric(job4_changes[top_group, ])
) |>
  pivot_longer(
    cols = c(OSF, Best_Match),
    names_to = "Group",
    values_to = "Change"
  )

p4 <- radar_data |>
  tidyplot(x = Cell_Type, y = Change, color = Group) |>
  add_line(group = Group) |>
  add_data_points(size = 3) |>
  adjust_colors(c("#E06681", "#8087E2")) |>
  adjust_title(paste0("OSF vs ", top_group, " (最相似组) 的细胞变化模式")) |>
  adjust_x_axis_title("细胞类型") |>
  adjust_y_axis_title("丰度变化 (Δ)") |>
  adjust_legend_title("") |>
  adjust_size(width = 150, height = 70) |>
  theme_tidyplot()

print(p4)

ggsave("osf_vs_best_match.pdf", p4, width = 12, height = 5)
cat("✓ 对比图已保存: osf_vs_best_match.pdf\n\n")

# =============================================================================
# 13. 可视化 5: 散点图矩阵
# =============================================================================

cat("📊 生成可视化 5: 所有试验组的细胞变化散点图...\n")

# 创建组合散点图数据
scatter_data <- lapply(rownames(job4_changes), function(group) {
  data.frame(
    Trial_Group = group,
    OSF_Change = osf_change,
    Trial_Change = as.numeric(job4_changes[group, ]),
    Cell_Type = cell_types,
    Correlation = correlations[group]
  )
}) |>
  bind_rows()

p5 <- scatter_data |>
  tidyplot(x = OSF_Change, y = Trial_Change) |>
  add_data_points(size = 2, alpha = 0.6) |>
  add_curve_fit(method = "lm") |>
  adjust_colors(colors_discrete_friendly) |>
  adjust_title("各试验组与OSF的细胞变化相关性散点图") |>
  adjust_x_axis_title("OSF 丰度变化") |>
  adjust_y_axis_title("试验组 丰度变化") |>
  remove_legend() |>
  adjust_size(width = 40, height = 30) |>
  # theme_tidyplot() |>
  split_plot(by = Trial_Group)

print(p5)

ggsave("scatter_correlation_matrix.pdf", p5, width = 12, height = 8)
cat("✓ 散点图矩阵已保存: scatter_correlation_matrix.pdf\n\n")

# =============================================================================
# 14. 导出结果表格
# =============================================================================

cat("💾 导出结果表格...\n")

# 导出相似性指标
write_csv(similarity_metrics, "similarity_metrics.csv")
cat("✓ 相似性指标已保存: similarity_metrics.csv\n")

# 导出所有变化值
write_csv(all_changes, "all_changes.csv")
cat("✓ 所有变化值已保存: all_changes.csv\n")

# 导出关键细胞类型
top_cells_summary <- data.frame(
  Cell_Type = names(osf_change_ranked)[1:5],
  OSF_Change = osf_change[names(osf_change_ranked)[1:5]],
  Abs_Change = osf_change_ranked[1:5]
)
write_csv(top_cells_summary, "top_changing_cells.csv")
cat("✓ 关键细胞类型已保存: top_changing_cells.csv\n\n")

# =============================================================================
# 15. 生成分析总结报告
# =============================================================================

cat("📝 生成分析总结...\n\n")

cat(
  "═══════════════════════════════════════════════════════════════════════════\n"
)
cat(
  "                      CIBERSORTx 相似性分析总结报告                         \n"
)
cat(
  "═══════════════════════════════════════════════════════════════════════════\n\n"
)

cat("【研究设计】\n")
cat("  • Job2: OSF疾病模型 vs Normal对照组\n")
cat("  • Job4: 6个药物干预组 vs PBS对照组\n")
cat("  • 分析目标: 识别哪些药物干预的细胞变化模式与OSF疾病模型最相似\n\n")

cat("【主要发现】\n")
cat("  1. 最相似的试验组 (正相关 - 加剧疾病特征):\n")
for (i in 1:min(3, sum(correlations_sorted > 0))) {
  cat(sprintf(
    "     • %s: r = %.3f\n",
    names(correlations_sorted)[i],
    correlations_sorted[i]
  ))
}

cat("\n  2. 最不相似的试验组 (负相关 - 可能逆转疾病特征):\n")
neg_corr <- correlations_sorted[correlations_sorted < 0]
if (length(neg_corr) > 0) {
  for (i in 1:min(3, length(neg_corr))) {
    idx <- length(correlations_sorted) - i + 1
    cat(sprintf(
      "     • %s: r = %.3f\n",
      names(correlations_sorted)[idx],
      correlations_sorted[idx]
    ))
  }
} else {
  cat("     • 无显著负相关组\n")
}

cat("\n  3. OSF疾病模型中变化最大的细胞类型:\n")
for (i in 1:5) {
  cell <- names(osf_change_ranked)[i]
  change <- osf_change[cell]
  direction <- ifelse(change > 0, "增加", "减少")
  cat(sprintf("     • %s: %s %.4f\n", cell, direction, abs(change)))
}

cat("\n【生物学解释建议】\n")
cat("  • 正相关组可能通过类似机制影响组织微环境\n")
cat("  • 负相关组可能具有治疗潜力，值得进一步验证\n")
cat("  • 关注在OSF中显著变化且在负相关组中被逆转的细胞类型\n\n")

cat("【输出文件】\n")
cat("  可视化图表:\n")
cat("    • heatmap_all_changes.pdf - 细胞丰度变化全景热图\n")
cat("    • correlation_barplot.pdf - 相关性柱状图\n")
cat("    • top5_cells_comparison.pdf - 关键细胞对比图\n")
cat("    • osf_vs_best_match.pdf - OSF vs 最相似组的对比图\n")
cat("    • scatter_correlation_matrix.pdf - 散点图矩阵\n")
cat("  数据表格:\n")
cat("    • similarity_metrics.csv - 相似性指标汇总\n")
cat("    • all_changes.csv - 所有组的细胞丰度变化\n")
cat("    • top_changing_cells.csv - 关键细胞类型\n\n")

cat(
  "═══════════════════════════════════════════════════════════════════════════\n"
)
cat("✅ 分析完成！所有结果已保存至当前工作目录\n")
cat(
  "═══════════════════════════════════════════════════════════════════════════\n"
)
