# =============================================================================
# CD8+ T 细胞活化分析（CD25/CD69 panel）
# 数据：projects/flow_cytometry/2026-09/0908/raw_data/CD25_CD69.xlsx
# 门控: Live/CD8+ -> CD25 (晚期活化) / CD69 (早期活化)
# 使用 r_functions/lib/compare_plot 组间比较函数 + r-pub-toolkit 的 utils 主题/配色系统
#
# panel 目的（见 projects/flow_cytometry/2026-09/0908/readme.md）：
# 评估 OT-I CD8+ T 细胞早期活化。DC-T 共培养后，排除死细胞、圈定 CD8+ T 细胞，
# 检测 CD69（早期活化）与 CD25（IL-2 受体上调，进一步活化）。
# 若某处理组 CD69+%、CD25+% 或 CD69+CD25+% 升高，说明该组 DC（对应 0902 那批
# DC 成熟度实验里的处理条件）更能诱导 T 细胞活化——这批数据是 0902 DC 刺激
# 实验的下游读出（DC 功能验证）。
#
# 实验设计：单一刺激反应比较，CONTROL / KSFM / LPS / TNF / HI / M50_6h / M50_24h，
# 各条件 n=3（无 batch 问题）。经与用户确认，vs-Control 统计比较以 CONTROL 为参照组。
# 统计检验不做多重检验校正（经用户确认），仅报告原始 p 值。
# =============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(rstatix)
library(readxl)

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

DATA_DIR   <- file.path(PROJ, "projects/flow_cytometry/2026-09/0908/raw_data")
OUTPUT_DIR <- file.path(PROJ, "projects/flow_cytometry/2026-09/0908/output")
FIG_DIR    <- file.path(OUTPUT_DIR, "figures")
STAT_DIR   <- file.path(OUTPUT_DIR, "stats")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. 读取 FlowJo 导出的宽表并整理成 tidy 长表
# -----------------------------------------------------------------------------
raw_wide <- read_excel(file.path(DATA_DIR, "CD25_CD69.xlsx"), sheet = 1)
names(raw_wide)[1] <- "sample_id"

# 样本名解析：CONTROL.fcs / CONTROL-2.fcs / CONTROL-3.fcs -> condition=CONTROL, replicate=1/2/3
# M50-24.fcs / M50-24-2.fcs -> condition=M50, timepoint=24h；M50-6.fcs 同理 timepoint=6h
parse_sample <- function(id) {
  stub <- sub("\\.fcs$", "", id)
  if (grepl("^M50-24", stub)) {
    condition <- "M50"
    timepoint <- "24h"
    rep_str <- sub("^M50-24-?", "", stub)
    replicate <- if (rep_str == "") 1L else as.integer(rep_str)
  } else if (grepl("^M50-6", stub)) {
    condition <- "M50"
    timepoint <- "6h"
    rep_str <- sub("^M50-6-?", "", stub)
    replicate <- if (rep_str == "") 1L else as.integer(rep_str)
  } else {
    parts <- strsplit(stub, "-")[[1]]
    condition <- parts[1]
    timepoint <- NA_character_
    replicate <- if (length(parts) == 1) 1L else as.integer(parts[2])
  }
  tibble(condition = condition, timepoint = timepoint, replicate = replicate)
}

meta <- purrr::map_dfr(raw_wide$sample_id, parse_sample)
raw_wide <- bind_cols(raw_wide, meta)

group_label_levels <- c("CONTROL", "KSFM", "LPS", "TNF", "HI", "M50_6h", "M50_24h")

raw_wide <- raw_wide %>%
  mutate(
    group_label = if_else(condition == "M50", paste0("M50_", timepoint), condition),
    group_label = factor(group_label, levels = group_label_levels)
  )

col_map <- tibble::tribble(
  ~raw_col,                                                    ~key,             ~data_type,   ~label,
  "T/Live/CD8+ | Geometric Mean (VL1-A :: CD25-A)",            "cd25_gmfi",      "gmfi",       "CD25 gMFI (a.u.)",
  "T/Live/CD8+ | Geometric Mean (RL1-A :: CD69-A)",            "cd69_gmfi",      "gmfi",       "CD69 gMFI (a.u.)",
  "T/Live/CD8+/Q1: CD25-A- , CD69-A+ | Freq. of Parent (%)",   "cd69_single_pos","pct",        "CD69+CD25- early activation (%)",
  "T/Live/CD8+/Q2: CD25-A+ , CD69-A+ | Freq. of Parent (%)",   "cd25cd69_dpos",  "pct",        "CD25+CD69+ full activation (%)"
)

stopifnot(all(col_map$raw_col %in% names(raw_wide)))

data_long <- raw_wide %>%
  select(sample_id, group_label, condition, timepoint, replicate, all_of(col_map$raw_col)) %>%
  pivot_longer(cols = all_of(col_map$raw_col), names_to = "raw_col", values_to = "value") %>%
  left_join(col_map, by = "raw_col")

message("各组样本量：")
print(data_long %>% distinct(group_label, sample_id) %>% count(group_label))

# -----------------------------------------------------------------------------
# 2. 各组 vs CONTROL 的关键比较
# -----------------------------------------------------------------------------
key_pairs <- list(
  c("CONTROL", "KSFM"),
  c("CONTROL", "LPS"),
  c("CONTROL", "TNF"),
  c("CONTROL", "HI"),
  c("CONTROL", "M50_6h"),
  c("CONTROL", "M50_24h")
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

save_fig <- function(p, stub, width = 9, height = 9.5, title = NULL) {
  g <- fix_panel_size(p, width = width, height = height)
  dw <- grid::convertWidth(sum(g$widths), "cm", valueOnly = TRUE)
  dh <- grid::convertHeight(sum(g$heights), "cm", valueOnly = TRUE)
  ggsave(file.path(FIG_DIR, paste0(stub, ".pdf")), g, width = dw, height = dh, units = "cm")
  ggsave(file.path(FIG_DIR, paste0(stub, ".png")), g, width = dw, height = dh, units = "cm", dpi = 300)
  save_plot_pptx(p, file.path(FIG_DIR, paste0(stub, ".pptx")), panel_width = width, panel_height = height, title = title)
  message("已生成图: ", stub)
}

plot_list <- list()
anova_results <- list()
ttest_results <- list()

for (i in seq_len(nrow(col_map))) {
  key     <- col_map$key[i]
  m_label <- col_map$label[i]
  cap     <- if (col_map$data_type[i] == "pct") 100 else Inf

  df_m <- data_long %>% filter(key == col_map$key[i])

  # --- 单因素 ANOVA 整体检验 ---
  av <- rstatix::anova_test(data = df_m, formula = value ~ group_label)
  anova_results[[key]] <- tibble(
    metric = key, label = m_label,
    statistic = av$F, df1 = av$DFn, df2 = av$DFd, p_value = av$p
  )

  # --- 关键成对 t 检验（各组 vs CONTROL，Welch's t-test）---
  wt <- purrr::map_dfr(key_pairs, function(pr) {
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

  p <- build_stat_plot(df_m, m_label, key_pairs, cap = cap) +
    labs(subtitle = kw_subtitle) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none",
      plot.subtitle = element_text(size = 9, color = "grey30")
    )

  plot_list[[key]] <- p
  save_fig(p, sprintf("%02d_%s", i, key), width = 9, height = 9.5, title = m_label)
}

# 总览拼图（4 个指标）
overview <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(title = "CD8+ T cell activation (CD25/CD69) vs CONTROL - overview")

ggsave(file.path(FIG_DIR, "00_overview.pdf"), overview, width = 24, height = 20, units = "cm")
ggsave(file.path(FIG_DIR, "00_overview.png"), overview, width = 24, height = 20, units = "cm", dpi = 300)
message("已生成总览拼图: 00_overview")

named_plots <- setNames(plot_list, col_map$label)
save_plots_pptx(named_plots, file.path(FIG_DIR, "all_metrics.pptx"), panel_width = 9, panel_height = 9.5)
message("已生成合并 PPTX: all_metrics.pptx")

# -----------------------------------------------------------------------------
# 3. 统计结果表（原始 p 值，不做多重检验校正）
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

write.csv(anova_tbl, file.path(STAT_DIR, "anova_overall.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(ttest_tbl, file.path(STAT_DIR, "pairwise_ttest_vs_control.csv"), row.names = FALSE, fileEncoding = "UTF-8")
message("已生成统计结果表")

# -----------------------------------------------------------------------------
# 4. 生成 Markdown 分析报告
# -----------------------------------------------------------------------------
report_path <- file.path(OUTPUT_DIR, "report.md")

sample_n_tbl <- data_long %>% distinct(group_label, sample_id) %>% count(group_label)

lines <- c(
  "# CD8+ T 细胞活化分析报告（CD25/CD69）",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. 实验背景",
  "",
  "本 panel（CD25/CD69）用于评估 **OT-I CD8+ T 细胞早期活化**。",
  "DC-T 共培养后，先排除死细胞，再圈定 CD8+ T 细胞，检测 CD69（早期活化标志物）",
  "和 CD25（IL-2 受体上调，反映进一步活化状态）表达（gMFI），以及 CD25/CD69 双染象限占比：",
  "",
  "- Q1（CD69+CD25-）：早期活化，未进入晚期活化阶段",
  "- Q2（CD25+CD69+）：完全活化（双阳性）",
  "",
  "若处理组的 CD69+%、CD25+% 或 CD69+CD25+% 升高，说明该组 DC 更能诱导 T 细胞活化——",
  "这批数据是 `projects/flow_cytometry/2026-09/0902/`（DC 成熟度）刺激实验的下游功能读出，",
  "两批数据共用同一套处理条件（KSFM/LPS/TNF/HI/M50_6h/M50_24h）。",
  "",
  "条件：**CONTROL**（未处理对照）/ KSFM（基础培养基对照）/ LPS（TLR4 阳性对照）/ TNF / HI / M50_6h / M50_24h，",
  "各组 n=3（无 batch 问题）。经确认，**vs-Control 统计比较以 CONTROL 为参照组**",
  "（KSFM 作为普通处理组之一与其他组一起比较，不作为参照）。",
  "",
  "各组样本量：",
  "",
  "| 组别 | n |",
  "|---|---|",
  paste0("| ", sample_n_tbl$group_label, " | ", sample_n_tbl$n, " |"),
  "",
  "## 2. 分析方法",
  "",
  "- 绘图使用项目库函数 `compare_plot()`，各组 n=3，自动选用 `pure_scatter` 策略（散点 + 均值线）。",
  "- 配色使用 `NPG` 期刊配色，主题使用 `theme_pub_base()`。",
  "- 每张图同时导出 PDF / PNG(300dpi) / PPTX 三种格式。",
  "- 统计检验（参数检验，均为原始 p 值，不做多重检验校正）：",
  "  - 整体差异：单因素 ANOVA（4 项指标独立进行）。",
  "  - 关键组间比较：独立样本 t 检验（Welch's t-test，不假设方差齐性），各组分别与 CONTROL 比较。",
  "  - **提示**：各组 n=3，t 检验自由度很低（Welch 近似 df 通常 2-4），",
  "    效应量估计和 p 值本身波动较大，请结合均值/散点分布与 ANOVA 整体检验判断，",
  "    不宜仅凭单个 p 值下结论。",
  "- 双阳性 %（Q2, CD25+CD69+）数值本身较低（约 0.7-3.5%），未接近 100% 硬顶，",
  "  比较括号无裁切问题，正常叠加显示。",
  "",
  "## 3. 整体差异（单因素 ANOVA）",
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
  "- `figures/00_overview.pdf/.png` — 4 项指标总览拼图",
  "- `figures/all_metrics.pptx` — 4 项指标可编辑合集",
  paste0(
    "- `figures/", sprintf("%02d_%s", seq_len(nrow(col_map)), col_map$key),
    ".pdf/.png/.pptx` — ", col_map$label, " 单指标图（vs CONTROL）"
  ),
  "- `stats/anova_overall.csv` — 整体单因素 ANOVA 检验结果",
  "- `stats/pairwise_ttest_vs_control.csv` — 各组 vs CONTROL 的 t 检验结果",
  ""
)

writeLines(lines, report_path, useBytes = TRUE)
message("已生成分析报告: ", report_path)

message("\n全部完成。输出目录: ", OUTPUT_DIR)
