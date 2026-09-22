# =============================================================================
# 流式 Mn 数据：不同组别模型小鼠细胞类型占比比较
# 数据：projects/flow_cytometry/2026-07/data/流式Mn.csv
# 使用 r_functions/lib/compare_plot 组间比较函数 + r-pub-toolkit 的 utils 主题/配色系统
# =============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(rstatix)

# 定位仓库根目录（本地/服务器通用）：优先 R_TOOLKIT_ROOT 环境变量，否则用 git
PROJ <- Sys.getenv("R_TOOLKIT_ROOT", unset = "")
if (!nzchar(PROJ)) PROJ <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE), error = function(e) "")
if (length(PROJ) != 1 || !nzchar(PROJ)) stop("找不到仓库根目录：请设置 R_TOOLKIT_ROOT 环境变量，或在仓库目录内运行")
source(file.path(PROJ, "init.R"))
source(file.path(PROJ, "utils/helpers.R"))
source(file.path(PROJ, "utils/palette_system.R"))
source(file.path(PROJ, "utils/theme_system.R"))
source(file.path(PROJ, "utils/panel_fix.R"))
source(file.path(PROJ, "utils/export_pptx.R"))
source(file.path(PROJ, "r_functions/lib/compare_plot/compare_plot_optimized.R"))

DATA_DIR   <- file.path(PROJ, "projects/flow_cytometry/2026-07/data")
OUTPUT_DIR <- file.path(PROJ, "projects/flow_cytometry/2026-07/output")
FIG_DIR    <- file.path(OUTPUT_DIR, "figures")
STAT_DIR   <- file.path(OUTPUT_DIR, "stats")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. 读取数据（文件为 GBK 编码）
# -----------------------------------------------------------------------------
raw <- read.csv(
  file.path(DATA_DIR, "流式Mn.csv"),
  fileEncoding = "GBK",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

group_levels <- c("Control", "OXA", "OXA_Mn", "OXA_Mn_Rm", "Rm", "Rm_HI")

metric_info <- tibble::tribble(
  ~raw_col,                 ~new_col,         ~label,
  "活细胞（%）",            "live_pct",       "活细胞占比 (%)",
  "CD45+（%）",              "cd45_pct",       "CD45+ (%)",
  "CD3+（%）",                "cd3_pct",        "CD3+ T细胞 (%)",
  "CD4+（%）",                "cd4_pct",        "CD4+ T细胞 (%)",
  "CD8+（%）",                "cd8_pct",        "CD8+ T细胞 (%)",
  "CD4+IL-17A+（%）",       "cd4_il17a_pct",  "CD4+IL-17A+ (Th17) (%)",
  "CD4+IFN-γ+（%）",        "cd4_ifng_pct",   "CD4+IFN-γ+ (Th1) (%)"
)

stopifnot(all(c("组别", "编号", metric_info$raw_col) %in% names(raw)))

data <- raw %>%
  rename(group = 组别, sample_id = 编号) %>%
  rename(!!!setNames(metric_info$raw_col, metric_info$new_col)) %>%
  mutate(group = factor(group, levels = group_levels))

message("样本量：")
print(table(data$group))

# -----------------------------------------------------------------------------
# 2. 长表化
# -----------------------------------------------------------------------------
data_long <- data %>%
  pivot_longer(
    cols = all_of(metric_info$new_col),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = metric_info$new_col),
    label  = metric_info$label[match(metric, metric_info$new_col)]
  )

# -----------------------------------------------------------------------------
# 3. 关键成对比较 + 逐指标绘图
# -----------------------------------------------------------------------------
key_pairs <- list(
  c("Control", "OXA"),
  c("Control", "OXA_Mn"),
  c("Control", "OXA_Mn_Rm"),
  c("Control", "Rm"),
  c("Control", "Rm_HI")
)

# 细胞占比不可能超过 100%：先用标准间距 (step_increase/y_expand) 出图，
# 若 5 层比较括号把 Y 轴自动撑到 100 以上（高基线指标如 CD8+ 容易出现），
# 再收紧间距重画一次，保证轴上限不超过 100。
build_stat_plot <- function(df_m, m_label, comparisons, cap = 100) {
  make_plot <- function(step, y_exp) {
    compare_plot(
      df_m,
      value.var = "value", group.by = "group",
      strategy = "auto",
      add_stat = "wilcox.test", comparisons = comparisons,
      stat_label = "p.signif", hide_ns = FALSE,
      step_increase = step, y_expand = y_exp,
      palette = "NPG", theme_use = theme_pub_base,
      title = m_label, xlab = "组别", ylab = m_label
    )
  }

  p <- make_plot(0.14, 0.65)
  top <- ggplot_build(p)$layout$panel_params[[1]]$y.range[2]
  if (top > cap) {
    # 收紧间距后，用硬上限顶到 cap，为最上层括号留出实际余量（而不是再按比例扩展）
    p <- make_plot(0.09, 0.02) +
      scale_y_continuous(limits = c(0, cap), expand = expansion(mult = c(0.05, 0)))
  }
  p
}

plot_list <- list()
kruskal_results <- list()
wilcox_results <- list()

for (i in seq_len(nrow(metric_info))) {
  m_col   <- metric_info$new_col[i]
  m_label <- metric_info$label[i]

  df_m <- data_long %>% filter(metric == m_col)

  # --- Kruskal-Wallis 整体检验 ---
  kw <- df_m %>% rstatix::kruskal_test(value ~ group)
  kruskal_results[[m_col]] <- tibble(
    metric = m_col, label = m_label,
    statistic = kw$statistic, df = kw$df,
    p_value = kw$p
  )

  # --- 关键成对 Wilcoxon 检验 ---
  wt <- purrr::map_dfr(key_pairs, function(pr) {
    d1 <- df_m$value[df_m$group == pr[1]]
    d2 <- df_m$value[df_m$group == pr[2]]
    res <- suppressWarnings(wilcox.test(d1, d2))
    tibble(
      metric = m_col, label = m_label,
      group1 = pr[1], group2 = pr[2],
      n1 = length(d1), n2 = length(d2),
      statistic = unname(res$statistic),
      p_value = res$p.value
    )
  })
  wilcox_results[[m_col]] <- wt

  # --- 绘图 ---
  kw_subtitle <- sprintf("Kruskal-Wallis, p = %.4f", kw$p)

  p <- build_stat_plot(df_m, m_label, key_pairs, cap = 100) +
    labs(subtitle = kw_subtitle) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none",
      plot.subtitle = element_text(size = 9, color = "grey30")
    )

  plot_list[[m_col]] <- p

  g <- fix_panel_size(p, width = 9, height = 9.5)
  dw <- grid::convertWidth(sum(g$widths), "cm", valueOnly = TRUE)
  dh <- grid::convertHeight(sum(g$heights), "cm", valueOnly = TRUE)

  file_stub <- sprintf("%02d_%s", i, m_col)
  ggsave(file.path(FIG_DIR, paste0(file_stub, ".pdf")), g, width = dw, height = dh, units = "cm")
  ggsave(file.path(FIG_DIR, paste0(file_stub, ".png")), g, width = dw, height = dh, units = "cm", dpi = 300)
  save_plot_pptx(
    p, file.path(FIG_DIR, paste0(file_stub, ".pptx")),
    panel_width = 9, panel_height = 9.5, title = m_label
  )

  message("✓ 已生成图: ", file_stub)
}

# -----------------------------------------------------------------------------
# 4. 总览拼图
# -----------------------------------------------------------------------------
overview <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(title = "各组细胞类型占比总览")

ggsave(
  file.path(FIG_DIR, "00_overview_combined.pdf"),
  overview, width = 36, height = 18, units = "cm"
)
ggsave(
  file.path(FIG_DIR, "00_overview_combined.png"),
  overview, width = 36, height = 18, units = "cm", dpi = 300
)
message("✓ 已生成总览拼图")

# 所有指标打包为一份多页可编辑 PPTX（每图一页，可在 PowerPoint 中继续调整）
named_plots <- setNames(plot_list, metric_info$label)
save_plots_pptx(
  named_plots, file.path(FIG_DIR, "all_metrics.pptx"),
  panel_width = 9, panel_height = 9.5
)
message("✓ 已生成合并 PPTX: all_metrics.pptx")

# -----------------------------------------------------------------------------
# 5. 统计结果表（含 BH 多重校正）
# -----------------------------------------------------------------------------
kruskal_tbl <- bind_rows(kruskal_results) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

wilcox_tbl <- bind_rows(wilcox_results) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_value < 0.0001 ~ "****",
      p_value < 0.001  ~ "***",
      p_value < 0.01   ~ "**",
      p_value < 0.05   ~ "*",
      TRUE             ~ "ns"
    )
  ) %>%
  arrange(metric, group1, group2)

write.csv(kruskal_tbl, file.path(STAT_DIR, "kruskal_overall.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(wilcox_tbl, file.path(STAT_DIR, "pairwise_wilcox.csv"), row.names = FALSE, fileEncoding = "UTF-8")

message("✓ 已生成统计结果表")

# -----------------------------------------------------------------------------
# 6. 生成 Markdown 分析报告
# -----------------------------------------------------------------------------
report_path <- file.path(OUTPUT_DIR, "report.md")

sig_star <- function(p) {
  case_when(
    p < 0.0001 ~ "\\*\\*\\*\\*",
    p < 0.001  ~ "\\*\\*\\*",
    p < 0.01   ~ "\\*\\*",
    p < 0.05   ~ "\\*",
    TRUE       ~ "ns"
  )
}

sample_n_tbl <- data %>% count(group)

lines <- c(
  "# 流式 Mn 数据分析报告：不同组别模型小鼠细胞类型占比比较",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. 实验背景",
  "",
  "本批流式数据为几种建模小鼠口腔组织（上唇/舌）的免疫细胞亚群占比检测结果，",
  "组别包括：",
  "",
  "- **Control**：正常对照",
  "- **OXA**：OXA 诱导造模",
  "- **OXA_Mn**：OXA 造模 + Mn 处理",
  "- **OXA_Mn_Rm**：OXA 造模 + Mn + Rm 联合处理",
  "- **Rm**：Rm（活菌）处理",
  "- **Rm_HI**：Rm 热灭活（Heat-Inactivated）对照",
  "",
  "检测指标（7 项，均为占比 %）：活细胞占比、CD45+、CD3+、CD4+、CD8+、",
  "CD4+IL-17A+（Th17）、CD4+IFN-γ+（Th1）。",
  "",
  "各组样本量：",
  "",
  "| 组别 | n |",
  "|---|---|",
  paste0("| ", sample_n_tbl$group, " | ", sample_n_tbl$n, " |"),
  "",
  "## 2. 分析方法",
  "",
  "- 绘图使用项目库函数 `compare_plot()`（`r_functions/lib/compare_plot/compare_plot_optimized.R`），",
  "  按每组 n=6 自动选用 `bar_points` 策略（均值柱 + SEM 误差线 + 全部散点）。",
  "- 配色使用 `palette_system.R` 的 `NPG` 期刊配色，主题使用 `theme_pub_base()`。",
  "- 每张图同时导出 PDF（矢量，投稿/印刷）、PNG（位图，预览/嵌入文档）、",
  "  PPTX（`export_pptx.R` 用 officer + rvg::dml() 生成，图形与文字在 PowerPoint 中可编辑）三种格式。",
  "- 统计检验：",
  "  - 整体差异：Kruskal-Wallis 检验（7 指标独立进行，BH 法多重校正）。",
  "  - 关键组间比较：Wilcoxon 秩和检验，各处理组分别与 Control 比较",
  "    （`OXA`、`OXA_Mn`、`OXA_Mn_Rm`、`Rm`、`Rm_HI` vs `Control`），",
  "    评估各建模/处理条件相对正常对照的整体偏离程度。",
  "  - 图中显著性星号为**未经多重校正**的探索性展示（常规惯例）；",
  "    `output/stats/pairwise_wilcox.csv` 中同时给出 BH 校正后 p 值供严谨判断。",
  "",
  "## 3. 整体差异（Kruskal-Wallis）",
  "",
  "| 指标 | 统计量 | p 值 | BH 校正 p 值 |",
  "|---|---|---|---|",
  sprintf(
    "| %s | %.3f | %.4f | %.4f |",
    kruskal_tbl$label, kruskal_tbl$statistic, kruskal_tbl$p_value, kruskal_tbl$p_adj_BH
  ),
  "",
  "## 4. 各指标结果与解读",
  ""
)

for (i in seq_len(nrow(metric_info))) {
  m_col   <- metric_info$new_col[i]
  m_label <- metric_info$label[i]
  file_stub <- sprintf("%02d_%s", i, m_col)

  df_m <- data_long %>% filter(metric == m_col)
  means_tbl <- df_m %>%
    group_by(group) %>%
    summarise(mean = mean(value), sem = sd(value) / sqrt(n()), .groups = "drop")

  kw_row <- kruskal_tbl %>% filter(metric == m_col)
  wt_rows <- wilcox_tbl %>% filter(metric == m_col)

  lines <- c(lines,
    paste0("### 4.", i, " ", m_label),
    "",
    paste0(
      "整体 Kruskal-Wallis: p = ", sprintf("%.4f", kw_row$p_value),
      "（BH 校正 p = ", sprintf("%.4f", kw_row$p_adj_BH), "）"
    ),
    "",
    "各组均值 ± SEM：",
    "",
    "| 组别 | 均值 (%) | SEM |",
    "|---|---|---|",
    sprintf("| %s | %.2f | %.2f |", means_tbl$group, means_tbl$mean, means_tbl$sem),
    "",
    "关键组间比较：",
    "",
    "| 比较 | p 值 | BH 校正 p 值 | 显著性 |",
    "|---|---|---|---|",
    sprintf(
      "| %s vs %s | %.4f | %.4f | %s |",
      wt_rows$group1, wt_rows$group2, wt_rows$p_value, wt_rows$p_adj_BH, wt_rows$significance
    ),
    "",
    paste0("图: `figures/", file_stub, ".pdf`（同目录下含 `.png` 位图版、`.pptx` 可编辑矢量版）"),
    ""
  )
}

lines <- c(lines,
  "## 5. 文件清单",
  "",
  "- `figures/00_overview_combined.pdf`（+ `.png`） — 7 指标总览拼图",
  "- `figures/all_metrics.pptx` — 7 指标可编辑矢量图合集（每图一页，PowerPoint 打开后图形/文字均可编辑）",
  paste0(
    "- `figures/", sprintf("%02d_%s", seq_len(nrow(metric_info)), metric_info$new_col),
    ".pdf`（+ `.png` + `.pptx`） — ", metric_info$label, " 单指标图"
  ),
  "- `stats/kruskal_overall.csv` — 整体 Kruskal-Wallis 检验结果",
  "- `stats/pairwise_wilcox.csv` — 关键组间 Wilcoxon 检验结果",
  ""
)

writeLines(lines, report_path, useBytes = TRUE)
message("✓ 已生成分析报告: ", report_path)

message("\n全部完成。输出目录: ", OUTPUT_DIR)
