# =============================================================================
# IHC CD8 染色批量分析：不同组别 CD8 染色强度组间比较
# 数据：projects/ihc/2026-07/data/FieldResults.csv（Fiji 宏批量输出，每图 9 视野）
# 只使用 Valid_Field == 1（组织占比 >= ~20%）的视野，每个视野作为一个观测
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
source(file.path(
  PROJ,
  "r_functions/lib/compare_plot/compare_plot_optimized.R"
))

# 本仓库（olp-microbiome-analysis）根目录：脚本本身在这个仓库里，用 git 定位
# （不能用 PROJ——PROJ 现在指向独立的 r-pub-toolkit 仓库，没有 projects/ 目录）
REPO_ROOT <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE), error = function(e) "")
if (length(REPO_ROOT) != 1 || !nzchar(REPO_ROOT)) stop("找不到本仓库根目录：请在仓库目录内运行脚本")

DATA_DIR <- file.path(REPO_ROOT, "projects/ihc/2026-07/data")
OUTPUT_DIR <- file.path(REPO_ROOT, "projects/ihc/2026-07/output")
FIG_DIR <- file.path(OUTPUT_DIR, "figures")
STAT_DIR <- file.path(OUTPUT_DIR, "stats")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. 读取数据 + 筛选有效视野
# -----------------------------------------------------------------------------
raw <- read.csv(
  file.path(DATA_DIR, "FieldResults.csv"),
  stringsAsFactors = FALSE
)
raw <- raw %>% filter(!is.na(Image_Name) & Image_Name != "")

group_levels <- c("C", "OXA", "OXA_Mn", "OXA_Rm")

data <- raw %>%
  filter(Valid_Field == 1) %>%
  mutate(
    group = sub("_xc_20X_.*", "", Image_Name),
    group = factor(group, levels = group_levels)
  )

message("有效视野样本量（每视野 = 一个观测）：")
print(table(data$group))
message(sprintf(
  "\n全部视野 %d 个，有效视野 %d 个（%.1f%%）",
  nrow(raw),
  nrow(data),
  100 * nrow(data) / nrow(raw)
))

# -----------------------------------------------------------------------------
# 2. 指标定义（Positive_Area_pct 为百分比，理论上限 100%；其余两项无固定上限）
# -----------------------------------------------------------------------------
metric_info <- tibble::tribble(
  ~raw_col            , ~new_col        , ~label              , ~cap ,
  "Positive_Area_pct" , "positive_pct"  , "DAB+ Area (%)"     ,  100 ,
  "Intensity"         , "dab_intensity" , "DAB Intensity"     , Inf  ,
  "Corrected_IOD"     , "corrected_iod" , "Corrected DAB IOD" , Inf
)

stopifnot(all(metric_info$raw_col %in% names(data)))

data <- data %>% rename(!!!setNames(metric_info$raw_col, metric_info$new_col))

data_long <- data %>%
  pivot_longer(
    cols = all_of(metric_info$new_col),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = metric_info$new_col),
    label = metric_info$label[match(metric, metric_info$new_col)]
  )

# -----------------------------------------------------------------------------
# 3. 关键成对比较（各组 vs Control）+ 自适应轴上限的绘图函数
# -----------------------------------------------------------------------------
key_pairs <- list(
  c("C", "OXA"),
  c("C", "OXA_Mn"),
  c("C", "OXA_Rm")
)

# 百分比类指标不可能超过 100%：先用标准间距出图，若自动扩展的 Y 轴上限超过 cap
# （对非百分比指标 cap = Inf，恒不触发），收紧间距重画一次并卡死硬上限。
build_stat_plot <- function(df_m, m_label, comparisons, cap = Inf) {
  make_plot <- function(step, y_exp) {
    compare_plot(
      df_m,
      value.var = "value",
      group.by = "group",
      strategy = "auto",
      add_stat = "wilcox.test",
      comparisons = comparisons,
      stat_label = "p.signif",
      hide_ns = FALSE,
      step_increase = step,
      y_expand = y_exp,
      palette = "NPG",
      theme_use = theme_pub_base,
      title = m_label,
      xlab = "Group",
      ylab = m_label
    )
  }

  p <- make_plot(0.14, 0.5)
  top <- ggplot_build(p)$layout$panel_params[[1]]$y.range[2]
  if (is.finite(cap) && top > cap) {
    p <- make_plot(0.09, 0.02) +
      scale_y_continuous(
        limits = c(0, cap),
        expand = expansion(mult = c(0.05, 0))
      )
  }
  p
}

plot_list <- list()
kruskal_results <- list()
wilcox_results <- list()

for (i in seq_len(nrow(metric_info))) {
  m_col <- metric_info$new_col[i]
  m_label <- metric_info$label[i]
  m_cap <- metric_info$cap[i]

  df_m <- data_long %>% filter(metric == m_col)

  # --- Kruskal-Wallis 整体检验 ---
  kw <- df_m %>% rstatix::kruskal_test(value ~ group)
  kruskal_results[[m_col]] <- tibble(
    metric = m_col,
    label = m_label,
    statistic = kw$statistic,
    df = kw$df,
    p_value = kw$p
  )

  # --- 关键成对 Wilcoxon 检验（各组 vs C）---
  wt <- purrr::map_dfr(key_pairs, function(pr) {
    d1 <- df_m$value[df_m$group == pr[1]]
    d2 <- df_m$value[df_m$group == pr[2]]
    res <- suppressWarnings(wilcox.test(d1, d2))
    tibble(
      metric = m_col,
      label = m_label,
      group1 = pr[1],
      group2 = pr[2],
      n1 = length(d1),
      n2 = length(d2),
      statistic = unname(res$statistic),
      p_value = res$p.value
    )
  })
  wilcox_results[[m_col]] <- wt

  # --- 绘图 ---
  kw_subtitle <- sprintf("Kruskal-Wallis, p = %.4f", kw$p)

  p <- build_stat_plot(df_m, m_label, key_pairs, cap = m_cap) +
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
  ggsave(
    file.path(FIG_DIR, paste0(file_stub, ".pdf")),
    g,
    width = dw,
    height = dh,
    units = "cm"
  )
  ggsave(
    file.path(FIG_DIR, paste0(file_stub, ".png")),
    g,
    width = dw,
    height = dh,
    units = "cm",
    dpi = 300
  )
  save_plot_pptx(
    p,
    file.path(FIG_DIR, paste0(file_stub, ".pptx")),
    panel_width = 9,
    panel_height = 9.5,
    title = m_label
  )

  message("✓ 已生成图: ", file_stub)
}

# -----------------------------------------------------------------------------
# 4. 总览拼图 + 合并 PPTX
# -----------------------------------------------------------------------------
overview <- wrap_plots(plot_list, ncol = 3) +
  plot_annotation(
    title = "CD8 IHC Staining Metrics Overview (Valid Fields Only)"
  )

ggsave(
  file.path(FIG_DIR, "00_overview_combined.pdf"),
  overview,
  width = 27,
  height = 18,
  units = "cm"
)
ggsave(
  file.path(FIG_DIR, "00_overview_combined.png"),
  overview,
  width = 27,
  height = 18,
  units = "cm",
  dpi = 300
)
message("✓ 已生成总览拼图")

named_plots <- setNames(plot_list, metric_info$label)
save_plots_pptx(
  named_plots,
  file.path(FIG_DIR, "all_metrics.pptx"),
  panel_width = 9,
  panel_height = 9.5
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
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(metric, group1, group2)

write.csv(
  kruskal_tbl,
  file.path(STAT_DIR, "kruskal_overall.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  wilcox_tbl,
  file.path(STAT_DIR, "pairwise_wilcox.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("✓ 已生成统计结果表")

# -----------------------------------------------------------------------------
# 6. 生成 Markdown 分析报告
# -----------------------------------------------------------------------------
report_path <- file.path(OUTPUT_DIR, "report.md")

sample_n_tbl <- data %>% count(group)
field_qc_tbl <- raw %>%
  mutate(
    group = factor(sub("_xc_20X_.*", "", Image_Name), levels = group_levels)
  ) %>%
  count(group, Valid_Field) %>%
  pivot_wider(
    names_from = Valid_Field,
    values_from = n,
    values_fill = 0,
    names_prefix = "valid_"
  )

lines <- c(
  "# IHC CD8 染色分析报告：不同组别染色指标组间比较",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. 数据背景",
  "",
  "第三批小鼠模型（OXA_Mn 实验）CD8 免疫组化染色，Fiji 宏批量分析。",
  "每张切片图像分为 9 个视野（F01–F09）逐视野统计 DAB 染色，",
  "本次分析以**视野**为观测单位（非按图像取均值），仅保留 `Valid_Field == 1`",
  "（组织占比 ≥ ~20% 的视野，过滤掉组织覆盖不足、量化不可靠的视野）。",
  "",
  "组别：",
  "",
  "- **C**：Control 对照",
  "- **OXA**：OXA 诱导造模",
  "- **OXA_Mn**：OXA 造模 + Mn 处理",
  "- **OXA_Rm**：OXA 造模 + Rm 处理",
  "",
  "各组有效视野样本量：",
  "",
  "| 组别 | 有效视野数 |",
  "|---|---|",
  paste0("| ", sample_n_tbl$group, " | ", sample_n_tbl$n, " |"),
  "",
  "各组视野有效性明细（`valid_0`=无效视野数，`valid_1`=有效视野数）：",
  "",
  "| 组别 | valid_0 | valid_1 |",
  "|---|---|---|",
  sprintf(
    "| %s | %d | %d |",
    field_qc_tbl$group,
    field_qc_tbl$valid_0,
    field_qc_tbl$valid_1
  ),
  "",
  "## 2. 分析方法",
  "",
  "- 绘图使用项目库函数 `compare_plot()`（`r_functions/lib/compare_plot/compare_plot_optimized.R`），",
  "  按各组有效视野数自动选择可视化策略（本数据集 max n 落在 10–20 区间，自动选用",
  "  `boxplot_points`：箱线图 + 全部散点）。",
  "- 配色使用 `palette_system.R` 的 `NPG` 期刊配色，主题使用 `theme_pub_base()`。",
  "- 每张图同时导出 PDF（矢量）、PNG（位图）、PPTX（`export_pptx.R`，可在 PowerPoint 中编辑）三种格式。",
  "- `Positive_Area_pct`（百分比指标）Y 轴上限硬性不超过 100%；`DAB_Intensity`、",
  "  `Corrected_DAB_IOD` 无固定上限，按数据自动缩放。",
  "- 统计检验：",
  "  - 整体差异：Kruskal-Wallis 检验（3 指标独立进行，BH 法多重校正）。",
  "  - 关键组间比较：Wilcoxon 秩和检验，各处理组分别与 Control（`C`）比较",
  "    （`OXA`、`OXA_Mn`、`OXA_Rm` vs `C`）。",
  "  - 图中显著性星号为**未经多重校正**的探索性展示；",
  "    `output/stats/pairwise_wilcox.csv` 中同时给出 BH 校正后 p 值。",
  "",
  "## 3. 整体差异（Kruskal-Wallis）",
  "",
  "| 指标 | 统计量 | p 值 | BH 校正 p 值 |",
  "|---|---|---|---|",
  sprintf(
    "| %s | %.3f | %.4f | %.4f |",
    kruskal_tbl$label,
    kruskal_tbl$statistic,
    kruskal_tbl$p_value,
    kruskal_tbl$p_adj_BH
  ),
  "",
  "## 4. 各指标结果",
  ""
)

for (i in seq_len(nrow(metric_info))) {
  m_col <- metric_info$new_col[i]
  m_label <- metric_info$label[i]
  file_stub <- sprintf("%02d_%s", i, m_col)

  df_m <- data_long %>% filter(metric == m_col)
  means_tbl <- df_m %>%
    group_by(group) %>%
    summarise(mean = mean(value), sem = sd(value) / sqrt(n()), .groups = "drop")

  kw_row <- kruskal_tbl %>% filter(metric == m_col)
  wt_rows <- wilcox_tbl %>% filter(metric == m_col)

  lines <- c(
    lines,
    paste0("### 4.", i, " ", m_label),
    "",
    paste0(
      "整体 Kruskal-Wallis: p = ",
      sprintf("%.4f", kw_row$p_value),
      "（BH 校正 p = ",
      sprintf("%.4f", kw_row$p_adj_BH),
      "）"
    ),
    "",
    "各组均值 ± SEM：",
    "",
    "| 组别 | 均值 | SEM |",
    "|---|---|---|",
    sprintf(
      "| %s | %.2f | %.2f |",
      means_tbl$group,
      means_tbl$mean,
      means_tbl$sem
    ),
    "",
    "关键组间比较：",
    "",
    "| 比较 | p 值 | BH 校正 p 值 | 显著性 |",
    "|---|---|---|---|",
    sprintf(
      "| %s vs %s | %.4f | %.4f | %s |",
      wt_rows$group1,
      wt_rows$group2,
      wt_rows$p_value,
      wt_rows$p_adj_BH,
      wt_rows$significance
    ),
    "",
    paste0("图: `figures/", file_stub, ".pdf`（同目录下含 `.png`、`.pptx`）"),
    ""
  )
}

lines <- c(
  lines,
  "## 5. 文件清单",
  "",
  "- `figures/00_overview_combined.pdf`（+ `.png`） — 3 指标总览拼图",
  "- `figures/all_metrics.pptx` — 3 指标可编辑矢量图合集（每图一页）",
  paste0(
    "- `figures/",
    sprintf("%02d_%s", seq_len(nrow(metric_info)), metric_info$new_col),
    ".pdf`（+ `.png` + `.pptx`） — ",
    metric_info$label,
    " 单指标图"
  ),
  "- `stats/kruskal_overall.csv` — 整体 Kruskal-Wallis 检验结果",
  "- `stats/pairwise_wilcox.csv` — 关键组间 Wilcoxon 检验结果",
  ""
)

writeLines(lines, report_path, useBytes = TRUE)
message("✓ 已生成分析报告: ", report_path)

message("\n全部完成。输出目录: ", OUTPUT_DIR)
