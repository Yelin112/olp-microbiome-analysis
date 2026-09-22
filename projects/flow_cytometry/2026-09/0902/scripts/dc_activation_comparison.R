# =============================================================================
# DC 活化/成熟状态分析（IECM 刺激实验）
# 数据：projects/flow_cytometry/2026-09/0902/raw_data/dc_activation_tidy.csv
# 门控: Live/CD11c+ (DC-like) -> CD80/CD86/MHC-II
# 使用 r_functions/lib/compare_plot 组间比较函数 + r-pub-toolkit 的 utils 主题/配色系统
#
# 实验设计（两套，分开处理）：
#   1) 分化时间进程：D0 / D3 / D7（各条件仅 1-2 个重复，纯描述性，不做正式统计）
#   2) 刺激处理反应：KSFM（对照）/ LPS / TNF / HI / M50_6h / M50_24h（vs KSFM）
#
# batch 说明：original 与 NewSample 数值高度接近，经与用户确认，
# 合并作为同一条件下的生物学重复处理（不再按 batch 拆分子图/统计）。
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

# 本仓库（olp-microbiome-analysis）根目录：脚本本身在这个仓库里，用 git 定位
# （不能用 PROJ——PROJ 现在指向独立的 r-pub-toolkit 仓库，没有 projects/ 目录）
REPO_ROOT <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE), error = function(e) "")
if (length(REPO_ROOT) != 1 || !nzchar(REPO_ROOT)) stop("找不到本仓库根目录：请在仓库目录内运行脚本")

DATA_DIR   <- file.path(REPO_ROOT, "projects/flow_cytometry/2026-09/0902/raw_data")
OUTPUT_DIR <- file.path(REPO_ROOT, "projects/flow_cytometry/2026-09/0902/output")
FIG_DIR    <- file.path(OUTPUT_DIR, "figures")
STAT_DIR   <- file.path(OUTPUT_DIR, "stats")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. 读取 tidy 数据
# -----------------------------------------------------------------------------
raw <- read.csv(
  file.path(DATA_DIR, "dc_activation_tidy.csv"),
  stringsAsFactors = FALSE
)

# batch(original/NewSample) + replicate 合并成唯一重复编号
# (用户确认：两个 batch 合并作为生物学重复，不拆分展示)
raw <- raw %>% mutate(rep_id = paste(batch, replicate, sep = "_"))

# group_label：M50 把 timepoint 拼进条件名 (M50_6h / M50_24h)，
# 避免和分化时间点组的 timepoint 语义混淆
raw <- raw %>%
  mutate(
    group_label = if_else(condition == "M50", paste0("M50_", timepoint), condition)
  )

time_levels <- c("D0", "D3", "D7")
stim_levels <- c("KSFM", "LPS", "TNF", "HI", "M50_6h", "M50_24h")

marker_label_map <- c(
  CD80 = "CD80 gMFI (a.u.)",
  CD86 = "CD86 gMFI (a.u.)",
  MHC2 = "MHC-II gMFI (a.u.)"
)

# gMFI 长表（单标志物几何均值）
gmfi <- raw %>%
  filter(metric == "geometric_mean") %>%
  mutate(
    marker = factor(marker, levels = names(marker_label_map)),
    marker_label = factor(marker_label_map[as.character(marker)], levels = marker_label_map)
  )

# 双阳性百分比长表（成熟度指标：CD80+CD86+ / CD86+MHC-II+）
dpos <- raw %>%
  filter(readout %in% c("CD80pos_CD86pos_pct", "CD86pos_MHC2pos_pct")) %>%
  mutate(
    pair_label = factor(
      case_when(
        readout == "CD80pos_CD86pos_pct" ~ "CD80+CD86+ (%)",
        readout == "CD86pos_MHC2pos_pct" ~ "CD86+MHC-II+ (%)"
      ),
      levels = c("CD80+CD86+ (%)", "CD86+MHC-II+ (%)")
    )
  )

# 拆分为 时间进程 / 刺激反应 两个子集
gmfi_time <- gmfi %>%
  filter(condition %in% time_levels) %>%
  mutate(group_label = factor(group_label, levels = time_levels))

gmfi_stim <- gmfi %>%
  filter(group_label %in% stim_levels) %>%
  mutate(group_label = factor(group_label, levels = stim_levels))

dpos_time <- dpos %>%
  filter(condition %in% time_levels) %>%
  mutate(group_label = factor(group_label, levels = time_levels))

dpos_stim <- dpos %>%
  filter(group_label %in% stim_levels) %>%
  mutate(group_label = factor(group_label, levels = stim_levels))

message("时间进程各组样本量 (gMFI, 每 marker):")
print(table(gmfi_time$group_label) / nlevels(gmfi_time$marker))
message("刺激反应各组样本量 (gMFI, 每 marker):")
print(table(gmfi_stim$group_label) / nlevels(gmfi_stim$marker))

# -----------------------------------------------------------------------------
# 2. 图 A：分化时间进程 D0->D3->D7（纯描述性，n=2，不做正式统计）
# -----------------------------------------------------------------------------
p_time_gmfi <- compare_plot(
  gmfi_time,
  value.var = "value", group.by = "group_label",
  split.by = "marker_label",
  strategy = "auto", add_trend = TRUE,
  palette = "NPG", theme_use = theme_pub_base,
  title = "DC maturation markers across differentiation (D0-D3-D7)",
  xlab = "Timepoint", ylab = "gMFI (a.u.)"
) + theme(legend.position = "none")

p_time_dpos <- compare_plot(
  dpos_time,
  value.var = "value", group.by = "group_label",
  split.by = "pair_label",
  strategy = "auto", add_trend = TRUE,
  palette = "NPG", theme_use = theme_pub_base,
  title = "Double-positive maturation markers across differentiation",
  xlab = "Timepoint", ylab = "% of parent"
) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0.05, 0.05))) +
  theme(legend.position = "none")

save_fig <- function(p, stub, width = 18, height = 9, title = NULL) {
  g <- fix_panel_size(p, width = width, height = height)
  dw <- grid::convertWidth(sum(g$widths), "cm", valueOnly = TRUE)
  dh <- grid::convertHeight(sum(g$heights), "cm", valueOnly = TRUE)
  ggsave(file.path(FIG_DIR, paste0(stub, ".pdf")), g, width = dw, height = dh, units = "cm")
  ggsave(file.path(FIG_DIR, paste0(stub, ".png")), g, width = dw, height = dh, units = "cm", dpi = 300)
  save_plot_pptx(p, file.path(FIG_DIR, paste0(stub, ".pptx")), panel_width = width, panel_height = height, title = title)
  message("已生成图: ", stub)
}

save_fig(p_time_gmfi, "06_time_course_gmfi", width = 20, height = 8, title = "Time course - gMFI")
save_fig(p_time_dpos, "07_time_course_double_positive", width = 15, height = 8, title = "Time course - double positive %")

# -----------------------------------------------------------------------------
# 3. 图 B：刺激处理反应，各组 vs KSFM 对照
# -----------------------------------------------------------------------------
key_pairs_stim <- list(
  c("KSFM", "LPS"),
  c("KSFM", "TNF"),
  c("KSFM", "HI"),
  c("KSFM", "M50_6h"),
  c("KSFM", "M50_24h")
)

# 百分比指标 Y 轴上限不能超过 100%；非百分比（gMFI）不设上限
build_stat_plot <- function(df_m, m_label, comparisons, cap = Inf) {
  make_plot <- function(step, y_exp) {
    compare_plot(
      df_m,
      value.var = "value", group.by = "group_label",
      strategy = "auto",
      add_stat = "t.test", comparisons = comparisons,
      stat_label = "p.signif", hide_ns = FALSE,
      step_increase = step, y_expand = y_exp,
      palette = "NPG", theme_use = theme_pub_base,
      title = m_label, xlab = "Group", ylab = m_label
    )
  }

  p <- make_plot(0.14, 0.5)
  top <- ggplot_build(p)$layout$panel_params[[1]]$y.range[2]
  if (is.finite(cap) && top > cap) {
    p <- make_plot(0.09, 0.02) +
      scale_y_continuous(limits = c(0, cap), expand = expansion(mult = c(0.05, 0)))
  }
  p
}

stim_metric_info <- tibble::tribble(
  ~key,              ~data_source, ~filter_col,    ~filter_val,             ~label,                    ~cap,  ~show_brackets,
  "cd80_gmfi",       "gmfi",       "marker",       "CD80",                  "CD80 gMFI (a.u.)",       Inf,   TRUE,
  "cd86_gmfi",       "gmfi",       "marker",       "CD86",                  "CD86 gMFI (a.u.)",       Inf,   TRUE,
  "mhc2_gmfi",       "gmfi",       "marker",       "MHC2",                  "MHC-II gMFI (a.u.)",     Inf,   TRUE,
  "cd80cd86_dpos",   "dpos",       "readout",      "CD80pos_CD86pos_pct",   "CD80+CD86+ (%)",         100,   FALSE,
  "cd86mhc2_dpos",   "dpos",       "readout",      "CD86pos_MHC2pos_pct",   "CD86+MHC-II+ (%)",       100,   FALSE
)
# show_brackets = FALSE：双阳性 % 指标里多数条件已接近 90-95%，硬顶 100% 的前提下
# 5 条 vs-KSFM 比较括号完全没有堆叠空间（会被裁切/整条丢失，实测确认），
# 所以这两个指标图上只保留散点+均值+KW 总体检验 subtitle，
# 完整的两两 Wilcoxon p 值仍然进 stats CSV 和 report.md 表格。

plot_list <- list()
anova_results <- list()
ttest_results <- list()

for (i in seq_len(nrow(stim_metric_info))) {
  key     <- stim_metric_info$key[i]
  m_label <- stim_metric_info$label[i]
  cap     <- stim_metric_info$cap[i]

  src <- if (stim_metric_info$data_source[i] == "gmfi") gmfi_stim else dpos_stim
  df_m <- src %>% filter(.data[[stim_metric_info$filter_col[i]]] == stim_metric_info$filter_val[i])

  # --- 单因素 ANOVA 整体检验 ---
  av <- rstatix::anova_test(data = df_m, formula = value ~ group_label)
  anova_results[[key]] <- tibble(
    metric = key, label = m_label,
    statistic = av$F, df1 = av$DFn, df2 = av$DFd, p_value = av$p
  )

  # --- 关键成对 t 检验（各组 vs KSFM，Welch's t-test）---
  wt <- purrr::map_dfr(key_pairs_stim, function(pr) {
    d1 <- df_m$value[df_m$group_label == pr[1]]
    d2 <- df_m$value[df_m$group_label == pr[2]]
    res <- suppressWarnings(t.test(d1, d2))
    tibble(
      metric = key, label = m_label,
      group1 = pr[1], group2 = pr[2],
      n1 = length(d1), n2 = length(d2),
      statistic = unname(res$statistic), p_value = res$p.value
    )
  })
  ttest_results[[key]] <- wt

  # --- 绘图 ---
  kw_subtitle <- sprintf("One-way ANOVA, p = %.4f", av$p)

  if (stim_metric_info$show_brackets[i]) {
    p <- build_stat_plot(df_m, m_label, key_pairs_stim, cap = cap)
  } else {
    # 不叠加两两比较括号（100% 硬顶下没有堆叠空间），仅散点+均值+KW subtitle
    p <- compare_plot(
      df_m,
      value.var = "value", group.by = "group_label",
      strategy = "auto",
      palette = "NPG", theme_use = theme_pub_base,
      title = m_label, xlab = "Group", ylab = m_label
    ) +
      scale_y_continuous(limits = c(0, cap), expand = expansion(mult = c(0.05, 0.05)))
  }

  p <- p +
    labs(subtitle = kw_subtitle) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none",
      plot.subtitle = element_text(size = 9, color = "grey30")
    )

  plot_list[[key]] <- p
  save_fig(p, sprintf("%02d_%s_stim", i, key), width = 9, height = 9.5, title = m_label)
}

# 总览拼图（5 个指标）
overview <- wrap_plots(plot_list, ncol = 3) +
  plot_annotation(title = "DC stimulation response vs KSFM control - overview")

ggsave(file.path(FIG_DIR, "00_overview_stim.pdf"), overview, width = 32, height = 20, units = "cm")
ggsave(file.path(FIG_DIR, "00_overview_stim.png"), overview, width = 32, height = 20, units = "cm", dpi = 300)
message("已生成总览拼图: 00_overview_stim")

named_plots <- setNames(plot_list, stim_metric_info$label)
save_plots_pptx(named_plots, file.path(FIG_DIR, "all_stim_metrics.pptx"), panel_width = 9, panel_height = 9.5)
message("已生成合并 PPTX: all_stim_metrics.pptx")

# -----------------------------------------------------------------------------
# 4. 统计结果表（原始 p 值，不做多重检验校正）
# -----------------------------------------------------------------------------
anova_tbl <- bind_rows(anova_results) %>%
  arrange(p_value)

ttest_tbl <- bind_rows(ttest_results) %>%
  mutate(
    significance = case_when(
      p_value < 0.0001 ~ "****",
      p_value < 0.001  ~ "***",
      p_value < 0.01   ~ "**",
      p_value < 0.05   ~ "*",
      TRUE             ~ "ns"
    )
  ) %>%
  arrange(metric, group1, group2)

write.csv(anova_tbl, file.path(STAT_DIR, "anova_overall_stim.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(ttest_tbl, file.path(STAT_DIR, "pairwise_ttest_vs_ksfm.csv"), row.names = FALSE, fileEncoding = "UTF-8")
message("已生成统计结果表")

# -----------------------------------------------------------------------------
# 5. 生成 Markdown 分析报告
# -----------------------------------------------------------------------------
report_path <- file.path(OUTPUT_DIR, "report.md")

sample_n_time <- gmfi_time %>% distinct(group_label, rep_id) %>% count(group_label)
sample_n_stim <- gmfi_stim %>% distinct(group_label, rep_id) %>% count(group_label)

lines <- c(
  "# DC 活化/成熟状态分析报告",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. 实验背景",
  "",
  "评估 IECM 对 DC（树突状细胞）成熟/活化的刺激能力。",
  "门控路径 `Live -> CD11c+ (DC-like)`，检测 CD80、CD86、MHC-II 表达（gMFI），",
  "以及 CD80/CD86 与 CD86/MHC-II 的双阳性象限占比（成熟度指标）。",
  "",
  "数据包含两套独立设计：",
  "",
  "- **分化时间进程**：D0 / D3 / D7，观察未刺激分化过程中标志物的基线变化。",
  "- **刺激处理反应**：KSFM（基础培养基对照）/ LPS（TLR4 阳性对照）/ TNF / HI / M50_6h / M50_24h，",
  "  各处理组与 KSFM 对照比较，评估处理是否促进 DC 成熟。",
  "",
  "**batch 处理说明**：原始数据含 `original` 与 `NewSample` 两个 batch，数值高度接近，",
  "经确认后**合并作为同一条件下的生物学重复**处理，不再按 batch 拆分展示或统计。",
  "",
  "各组样本量（合并 batch 后）：",
  "",
  "时间进程组：",
  "",
  "| 组别 | n |",
  "|---|---|",
  paste0("| ", sample_n_time$group_label, " | ", sample_n_time$n, " |"),
  "",
  "刺激反应组：",
  "",
  "| 组别 | n |",
  "|---|---|",
  paste0("| ", sample_n_stim$group_label, " | ", sample_n_stim$n, " |"),
  "",
  "## 2. 分析方法",
  "",
  "- 绘图使用项目库函数 `compare_plot()`，各组 n 较小（时间进程 n=2，刺激反应 n=3），",
  "  自动选用 `pure_scatter` 策略（散点 + 均值线）。",
  "- 配色使用 `NPG` 期刊配色，主题使用 `theme_pub_base()`。",
  "- 每张图同时导出 PDF / PNG(300dpi) / PPTX 三种格式。",
  "- **时间进程组（D0/D3/D7）为纯描述性展示**，不做正式统计检验（n=2 不足以支撑）。",
  "- **刺激反应组**统计检验（参数检验，均为原始 p 值，不做多重检验校正）：",
  "  - 整体差异：单因素 ANOVA（5 项指标独立进行）。",
  "  - 关键组间比较：独立样本 t 检验（Welch's t-test，不假设方差齐性），",
  "    各处理组分别与 KSFM 对照比较。",
  "  - **提示**：各组 n=3，t 检验自由度很低（Welch 近似 df 通常 2-4），",
  "    效应量估计和 p 值本身波动较大，请结合均值/散点分布与 ANOVA 整体检验判断，",
  "    不宜仅凭单个 p 值下结论。",
  "  - **双阳性 % 两项指标（CD80+CD86+、CD86+MHC-II+）图上不叠加两两比较括号**：",
  "    多数条件的双阳性比例已接近 90-95%，在 100% 硬顶约束下 5 条 vs-KSFM 括号",
  "    没有堆叠空间（实测会被裁切/丢失），因此这两张图仅保留散点+均值+ANOVA 总体检验，",
  "    完整两两 t 检验结果见统计表。",
  "",
  "## 3. 整体差异（单因素 ANOVA，刺激反应组）",
  "",
  "| 指标 | F 统计量 | p 值 |",
  "|---|---|---|",
  sprintf(
    "| %s | %.3f | %.4f |",
    anova_tbl$label, anova_tbl$statistic, anova_tbl$p_value
  ),
  "",
  "## 4. 文件清单",
  "",
  "- `figures/06_time_course_gmfi.pdf/.png/.pptx` — 分化时间进程 gMFI（CD80/CD86/MHC-II 分面）",
  "- `figures/07_time_course_double_positive.pdf/.png/.pptx` — 分化时间进程双阳性 %",
  "- `figures/00_overview_stim.pdf/.png` — 刺激反应 5 项指标总览拼图",
  "- `figures/all_stim_metrics.pptx` — 刺激反应 5 项指标可编辑合集",
  paste0(
    "- `figures/", sprintf("%02d_%s_stim", seq_len(nrow(stim_metric_info)), stim_metric_info$key),
    ".pdf/.png/.pptx` — ", stim_metric_info$label, " 单指标图（vs KSFM）"
  ),
  "- `stats/anova_overall_stim.csv` — 刺激反应组整体单因素 ANOVA 检验结果",
  "- `stats/pairwise_ttest_vs_ksfm.csv` — 刺激反应组各处理 vs KSFM 的 t 检验结果",
  ""
)

writeLines(lines, report_path, useBytes = TRUE)
message("已生成分析报告: ", report_path)

message("\n全部完成。输出目录: ", OUTPUT_DIR)
