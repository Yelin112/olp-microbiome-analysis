# =============================================================================
# Control vs R 组对比分析
# 使用 Tidyplots 主题管理系统
# =============================================================================

library(tidyverse)
library(tidyplots)

# -----------------------------------------------------------------------------
# 1. 加载主题系统
# -----------------------------------------------------------------------------
# 将之前提供的完整主题系统代码保存为 'tidyplot_themes.R'
source("Scripts\\tidyplot_themes.R")

# 或者如果直接复制在同一个脚本中，确保主题系统的代码在这之前已经运行

# -----------------------------------------------------------------------------
# 2. 读取数据
# -----------------------------------------------------------------------------
data <- read.csv("data\\1-29流式.csv")

# -----------------------------------------------------------------------------
# 3. 数据预处理
# -----------------------------------------------------------------------------
data <- data %>%
  mutate(
    # 计算 CD8/CD4 比值
    CD8_div_CD4 = CD8_pct_CD3 / CD4_pct_CD3,
    # 确保 Group 的顺序：Control 在前，R 在后
    Group = factor(Group, levels = c("Control", "R"))
  ) %>%
  # 🔑 关键：只保留 Control 和 R 组
  filter(Group %in% c("Control", "R"))

# 检查过滤后的数据
message(paste("📊 Total samples after filtering:", nrow(data)))
message(paste(
  "📊 Groups included:",
  paste(unique(data$Group), collapse = ", ")
))

# -----------------------------------------------------------------------------
# 4. 定义需要分析的指标
# -----------------------------------------------------------------------------
metrics_map <- list(
  "CD3_pct" = "CD3+ T Cell Proportion (%)",
  "CD4_pct_CD3" = "CD4+ T Cell (% of CD3)",
  "CD8_pct_CD3" = "CD8+ T Cell (% of CD3)",
  "CD8_div_CD4" = "CD8 / CD4 Ratio",
  "IL17_pct_CD4" = "IL-17+ (% of CD4)",
  "IFNg_pct_CD4" = "IFN-γ+ (% of CD4)",
  "IL17_MFI" = "IL-17 MFI (Intensity)",
  "IFNg_MFI" = "IFN-γ MFI (Intensity)"
)

# -----------------------------------------------------------------------------
# 5. 创建主输出文件夹
# -----------------------------------------------------------------------------
main_output_dir <- "Control_vs_R_Analysis"
if (!dir.exists(main_output_dir)) {
  dir.create(main_output_dir)
  message(paste("📁 Created main output folder:", main_output_dir))
}

# -----------------------------------------------------------------------------
# 6. 获取所有组织类型
# -----------------------------------------------------------------------------
tissues <- unique(data$Tissue)
message(paste("🔬 Tissues to analyze:", paste(tissues, collapse = ", ")))

# -----------------------------------------------------------------------------
# 7. 双重循环：组织 × 指标
# -----------------------------------------------------------------------------
cat("\n")
cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║           开始 Control vs R 对比分析                      ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

for (tis in tissues) {
  # --- 步骤 A: 为每个组织创建子文件夹 ---
  folder_name <- file.path(main_output_dir, gsub(" ", "_", tis))

  if (!dir.exists(folder_name)) {
    dir.create(folder_name, recursive = TRUE)
  }

  message(paste("\n👉 Processing Tissue:", tis))
  message(paste("   Output folder:", folder_name))

  # --- 步骤 B: 筛选该组织的数据 ---
  df_tissue <- data %>%
    filter(Tissue == tis)

  # 检查该组织中 Control 和 R 组的样本量
  sample_counts <- df_tissue %>%
    count(Group) %>%
    arrange(Group)

  message("   Sample size:")
  for (i in 1:nrow(sample_counts)) {
    message(sprintf(
      "     - %s: n = %d",
      sample_counts$Group[i],
      sample_counts$n[i]
    ))
  }

  # --- 步骤 C: 循环绘制每个指标 ---
  for (metric_col in names(metrics_map)) {
    metric_label <- metrics_map[[metric_col]]

    # 检查数据是否全为空
    if (all(is.na(df_tissue[[metric_col]]))) {
      message(paste("   ⚠️  Skipping", metric_col, "- No data available"))
      next
    }

    # 过滤掉当前指标为 NA 的行
    df_plot <- df_tissue %>%
      filter(!is.na(!!sym(metric_col)))

    # 检查过滤后是否还有数据
    if (nrow(df_plot) == 0) {
      message(paste("   ⚠️  Skipping", metric_col, "- All values are NA"))
      next
    }

    # 检查两组是否都有数据
    groups_with_data <- df_plot %>%
      group_by(Group) %>%
      summarise(n = n(), .groups = 'drop') %>%
      pull(Group)

    if (length(groups_with_data) < 2) {
      message(paste(
        "   ⚠️  Skipping",
        metric_col,
        "- Insufficient groups with data"
      ))
      next
    }

    # --- 开始绘图 ---
    p <- df_plot %>%
      tidyplot(x = Group, y = !!sym(metric_col), color = Group) %>%
      add_mean_bar(alpha = 0.6) %>%
      add_sem_errorbar() %>%
      add_data_points_beeswarm(
        alpha = 0.7,
        white_border = TRUE,
        size = 2.5
      ) %>%
      # 添加统计检验（Control vs R）
      {
        tryCatch(
          add_test_pvalue(
            .,
            ref.group = "Control",
            comparisons = list(c("Control", "R")),
            padding_top = 0.12,
            label.size = 4,
            bracket.size = 0.5
          ),
          error = function(e) {
            message(paste("   ℹ️  Could not calculate p-value for", metric_col))
            message(paste("      Reason:", e$message))
            .
          }
        )
      } %>%
      adjust_x_axis_title("Treatment Group") %>%
      adjust_y_axis_title(metric_label) %>%
      adjust_title(paste(tis, ":", metric_label)) %>%
      adjust_size(width = 70, height = 85) %>%
      apply_tidyplot_theme(MY_LAB_THEME, palette = c("#E88471FF", "#39B185FF"))

    # --- 保存 PDF ---
    file_name <- file.path(folder_name, paste0(metric_col, ".pdf"))
    save_plot(p, file_name)

    message(paste("   ✓ Generated:", metric_col))
  }

  message(paste("   ✅ Completed:", tis, "\n"))
}

cat("\n")
cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║           ✅ 所有图表生成完毕！                           ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")
cat(paste("\n📂 Results saved in:", main_output_dir, "\n\n"))

# -----------------------------------------------------------------------------
# 8. 生成分析报告摘要
# -----------------------------------------------------------------------------
generate_analysis_summary <- function(data, output_dir) {
  summary_file <- file.path(output_dir, "Analysis_Summary.txt")

  sink(summary_file)

  cat("═══════════════════════════════════════════════════════════\n")
  cat("         Control vs R Group Analysis Summary\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  cat("─────────────────────────────────────────────────────────\n")
  cat("1. Sample Overview\n")
  cat("─────────────────────────────────────────────────────────\n\n")

  overall_counts <- data %>%
    count(Group) %>%
    arrange(Group)

  cat("Total Samples:\n")
  for (i in 1:nrow(overall_counts)) {
    cat(sprintf(
      "  • %s: n = %d\n",
      overall_counts$Group[i],
      overall_counts$n[i]
    ))
  }

  cat("\n─────────────────────────────────────────────────────────\n")
  cat("2. Tissue Distribution\n")
  cat("─────────────────────────────────────────────────────────\n\n")

  tissue_summary <- data %>%
    group_by(Tissue, Group) %>%
    summarise(n = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Group, values_from = n, values_fill = 0)

  print(tissue_summary)

  cat("\n─────────────────────────────────────────────────────────\n")
  cat("3. Metrics Analyzed\n")
  cat("─────────────────────────────────────────────────────────\n\n")

  for (i in 1:length(metrics_map)) {
    cat(sprintf("  %d. %s\n", i, metrics_map[[i]]))
  }

  cat("\n─────────────────────────────────────────────────────────\n")
  cat("4. Output Structure\n")
  cat("─────────────────────────────────────────────────────────\n\n")

  cat(paste("Main folder:", output_dir, "\n"))
  for (tissue in unique(data$Tissue)) {
    cat(sprintf("  ├── %s/\n", gsub(" ", "_", tissue)))
    cat(sprintf("  │   ├── CD3_pct.pdf\n"))
    cat(sprintf("  │   ├── CD4_pct_CD3.pdf\n"))
    cat(sprintf("  │   └── ... (total %d metrics)\n", length(metrics_map)))
  }

  cat("\n─────────────────────────────────────────────────────────\n")
  cat("5. Statistical Testing\n")
  cat("─────────────────────────────────────────────────────────\n\n")

  cat("  Method: Wilcoxon rank sum test (default)\n")
  cat("  Comparison: Control (reference) vs R\n")
  cat("  Significance levels:\n")
  cat("    *    p < 0.05\n")
  cat("    **   p < 0.01\n")
  cat("    ***  p < 0.001\n")
  cat("    **** p < 0.0001\n")
  cat("    ns   not significant\n")

  cat("\n═══════════════════════════════════════════════════════════\n")
  cat("                    Analysis Complete\n")
  cat("═══════════════════════════════════════════════════════════\n")

  sink()

  message(paste("\n📄 Summary report saved:", summary_file))
}

# 生成摘要报告
generate_analysis_summary(data, main_output_dir)

# -----------------------------------------------------------------------------
# 9. 可选：创建组合图（所有组织的某个指标对比）
# -----------------------------------------------------------------------------
create_multi_tissue_comparison <- function(
  data,
  metric_col,
  metric_label,
  output_dir
) {
  # 过滤有效数据
  df_plot <- data %>%
    filter(!is.na(!!sym(metric_col)))

  if (nrow(df_plot) == 0) {
    message(paste("Skipping multi-tissue plot for", metric_col, "- No data"))
    return(NULL)
  }

  # 创建分面图
  p <- df_plot %>%
    tidyplot(x = Group, y = !!sym(metric_col), color = Group) %>%
    add_mean_bar(alpha = 0.6) %>%
    add_sem_errorbar() %>%
    add_data_points_beeswarm(alpha = 0.7, white_border = TRUE, size = 2) %>%
    adjust_x_axis_title("Treatment Group") %>%
    adjust_y_axis_title(metric_label) %>%
    adjust_title(paste("Multi-Tissue Comparison:", metric_label)) %>%
    adjust_size(width = 180, height = 120) %>%
    theme_nature(font_size = 8, palette = c("#E64B35", "#4DBBD5"))

  # 添加分面
  p <- p + facet_wrap(~Tissue, scales = "free_y", ncol = 3)

  # 保存
  file_name <- file.path(output_dir, paste0("MultiTissue_", metric_col, ".pdf"))
  save_plot(p, file_name)

  message(paste("✓ Created multi-tissue comparison for", metric_col))
}

# 可选：为关键指标创建多组织对比图
cat("\n📊 Creating multi-tissue comparison plots...\n\n")

comparison_dir <- file.path(main_output_dir, "Multi_Tissue_Comparisons")
if (!dir.exists(comparison_dir)) {
  dir.create(comparison_dir, recursive = TRUE)
}

# 选择几个关键指标进行多组织对比
key_metrics <- c("CD3_pct", "CD8_div_CD4", "IL17_pct_CD4", "IFNg_pct_CD4")

for (metric in key_metrics) {
  if (metric %in% names(metrics_map)) {
    create_multi_tissue_comparison(
      data,
      metric,
      metrics_map[[metric]],
      comparison_dir
    )
  }
}

message("\n✅ All analyses completed successfully!\n")
