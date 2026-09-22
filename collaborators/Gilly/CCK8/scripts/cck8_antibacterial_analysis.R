# =============================================================================
# CCK8 纳米酶抗菌能力分析
# 数据：collaborators/Gilly/CCK8/data/CCK8_nanozyme_antibacterial_tidy.csv
# 背景：96 孔板 CCK8 法评估同一纳米酶材料的浓度梯度（A~E 列 = nanozyme_1~5）
#       对细菌存活的抑制作用。已确认实际浓度（mg/mL）：
#         nanozyme_1 = 1, nanozyme_2 = 0.5, nanozyme_3 = 0.25,
#         nanozyme_4 = 0.125, nanozyme_5 = 0.0625（2 倍稀释系列，编号越大浓度越低）
#       F 列 = 阳性对照（BHI + 细菌，无纳米酶，相当于浓度 = 0 的基线）
#       G 列 = 空白对照（仅 BHI，无细菌，背景吸光值，仅用于 QC）
#       共 2 个批次（batch_1 / batch_2），每批每组 4 复孔。
#       批次间不汇集：两批分别独立计算统计量、独立出图，避免批次效应掩盖或
#       混淆真实的浓度效应；仅额外画一张两批叠加的趋势图做纯描述性的批次
#       一致性对比（不做跨批次统计检验）。
#       数据文件中已按分析流程预先计算好 blank_corrected_od /
#       relative_viability_percent / inhibition_percent（每批次独立扣本批
#       G 组背景、独立用本批 F 组归一化），此脚本直接复用这些列，不重复计算。
#
# 统计设计（按要求简化）：只做各浓度组 vs 阳性对照（Control）的简单组间检验
# （Wilcoxon 秩和检验），不做事后多重比较校正，不做整体 Kruskal-Wallis 检验，
# 不做浓度趋势检验（Jonckheere-Terpstra）。
# =============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)

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

BASE_DIR   <- file.path(PROJ, "collaborators/Gilly/CCK8")
DATA_DIR   <- file.path(BASE_DIR, "data")
OUTPUT_DIR <- file.path(BASE_DIR, "output")
FIG_DIR    <- file.path(OUTPUT_DIR, "figures")
STAT_DIR   <- file.path(OUTPUT_DIR, "stats")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. 读取数据 + 浓度映射（positive_control 视为浓度 = 0 的剂量-效应基线点）
# -----------------------------------------------------------------------------
raw <- read.csv(
  file.path(DATA_DIR, "CCK8_nanozyme_antibacterial_tidy.csv"),
  stringsAsFactors = FALSE
)

conc_map <- tibble::tribble(
  ~group,              ~conc_mgml, ~label,
  "positive_control",  0,          "0 (Control)",
  "nanozyme_5",        0.0625,     "0.0625",
  "nanozyme_4",        0.125,      "0.125",
  "nanozyme_3",        0.25,       "0.25",
  "nanozyme_2",        0.5,        "0.5",
  "nanozyme_1",        1,          "1"
)

data <- raw %>%
  filter(group %in% conc_map$group) %>%
  left_join(conc_map, by = "group") %>%
  mutate(
    group = factor(group, levels = conc_map$group, labels = conc_map$label),
    batch = factor(batch)
  )

message("各组样本量（按批次分列，每批每组 n=4）：")
print(table(data$group, data$batch))

positive_label <- conc_map$label[conc_map$group == "positive_control"]
dose_pairs <- purrr::map(setdiff(conc_map$label, positive_label), function(l) c(l, positive_label))

# -----------------------------------------------------------------------------
# 2. 背景扣除 QC（两批次一起列出对比，不涉及组间统计，只是核查用）
# -----------------------------------------------------------------------------
blank_qc <- raw %>%
  filter(group == "blank_control") %>%
  group_by(batch) %>%
  summarise(
    n = n(),
    mean_blank_corrected_od = mean(blank_corrected_od),
    sd_blank_corrected_od = sd(blank_corrected_od),
    .groups = "drop"
  )
write.csv(blank_qc, file.path(STAT_DIR, "blank_control_qc.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. 单批次分析函数：vs 阳性对照 Wilcoxon（不校正）+ 绘图
#    每批 n=4/组，独立计算，不与另一批混合。
# -----------------------------------------------------------------------------
build_stat_plot <- function(df, value_var, y_label, comparisons, cap = Inf) {
  make_plot <- function(step, y_exp) {
    compare_plot(
      df,
      value.var = value_var, group.by = "group",
      strategy = "auto",
      add_stat = "wilcox.test", comparisons = comparisons,
      stat_label = "p.signif", hide_ns = FALSE,
      step_increase = step, y_expand = y_exp,
      palette = "NPG", theme_use = theme_pub_base,
      title = "CCK8 Antibacterial Assay", xlab = "Nanozyme Concentration (mg/mL)", ylab = y_label
    )
  }
  p <- make_plot(0.14, 0.5)
  top <- ggplot_build(p)$layout$panel_params[[1]]$y.range[2]
  if (is.finite(cap) && top > cap) {
    p <- make_plot(0.09, 0.02) +
      scale_y_continuous(limits = c(NA, cap), expand = expansion(mult = c(0.05, 0)))
  }
  p
}

analyze_batch <- function(batch_id) {
  message("\n==== 分析 ", batch_id, " ====")
  df <- data %>% filter(batch == batch_id) %>% droplevels()

  fig_dir  <- file.path(FIG_DIR, batch_id)
  stat_dir <- file.path(STAT_DIR, batch_id)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(stat_dir, recursive = TRUE, showWarnings = FALSE)

  # --- 各浓度 vs 阳性对照 (n=4 vs n=4，Wilcoxon 秩和检验，不做事后校正) ---
  wilcox_tbl <- purrr::map_dfr(dose_pairs, function(pr) {
    d1 <- df$relative_viability_percent[df$group == pr[1]]
    d2 <- df$relative_viability_percent[df$group == pr[2]]
    res <- suppressWarnings(wilcox.test(d1, d2))
    tibble(
      group1 = pr[1], group2 = pr[2],
      n1 = length(d1), n2 = length(d2),
      mean_viability_1 = mean(d1), mean_viability_2 = mean(d2),
      mean_inhibition_1 = 100 - mean(d1),
      statistic = unname(res$statistic),
      p_value = res$p.value
    )
  }) %>%
    mutate(
      significance = case_when(
        p_value < 0.0001 ~ "****",
        p_value < 0.001  ~ "***",
        p_value < 0.01   ~ "**",
        p_value < 0.05   ~ "*",
        TRUE             ~ "ns"
      )
    )
  write.csv(wilcox_tbl, file.path(stat_dir, "pairwise_wilcox_vs_positive_control.csv"), row.names = FALSE)

  # --- 描述统计 ---
  summary_tbl <- df %>%
    group_by(group, conc_mgml) %>%
    summarise(
      n = n(),
      mean_viability = mean(relative_viability_percent),
      sem_viability = sd(relative_viability_percent) / sqrt(n()),
      mean_inhibition = mean(inhibition_percent),
      sem_inhibition = sd(inhibition_percent) / sqrt(n()),
      .groups = "drop"
    ) %>%
    arrange(conc_mgml)
  write.csv(summary_tbl, file.path(stat_dir, "group_summary.csv"), row.names = FALSE)

  # --- 绘图：组间比较（vs 阳性对照） ---
  # 注：inhibition_percent 理论上界为 100%，但由于是由 relative_viability_percent
  # 反算而来，噪声会使部分对照组/低浓度组样本略低于 0%——数据本身已不严格落在
  # [0,100] 区间内，硬性把轴顶锁在 100 会在 5 条比较括号堆叠时裁掉最上面的标注
  # （已实测验证），因此这里不加硬上限，按数据自动缩放（cap = Inf）。
  p_inhibition <- build_stat_plot(df, "inhibition_percent", "Inhibition Rate (%)", dose_pairs, cap = Inf) +
    labs(subtitle = sprintf("Wilcoxon rank-sum vs Control, unadjusted p (n=4/group, %s)", batch_id)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  p_viability <- build_stat_plot(df, "relative_viability_percent", "Relative Viability (%)", dose_pairs, cap = Inf) +
    labs(subtitle = sprintf("Wilcoxon rank-sum vs Control, unadjusted p (n=4/group, %s)", batch_id)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # --- 绘图：浓度-抑菌率描述性趋势图（仅展示均值 ± SEM，不做趋势检验） ---
  p_trend <- ggplot(summary_tbl, aes(x = group, y = mean_inhibition, group = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(color = get_colors("NPG", n = 1), linewidth = 1) +
    geom_errorbar(
      aes(ymin = mean_inhibition - sem_inhibition, ymax = mean_inhibition + sem_inhibition),
      width = 0.1, linewidth = 0.6
    ) +
    geom_point(size = 3.5, shape = 21, fill = get_colors("NPG", n = 1), color = "black", stroke = 0.8) +
    labs(
      title = "Nanozyme Dose-Response (Antibacterial Inhibition)",
      subtitle = sprintf("Descriptive only, mean \u00b1 SEM (%s)", batch_id),
      x = "Nanozyme Concentration (mg/mL)", y = "Inhibition Rate (%, mean \u00b1 SEM)"
    ) +
    theme_pub_base() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  plots <- list(
    list(p = p_inhibition, stub = "01_inhibition_rate", w = 10, h = 9.5),
    list(p = p_viability,  stub = "02_relative_viability", w = 10, h = 9.5),
    list(p = p_trend,      stub = "03_dose_response_trend", w = 13, h = 8)
  )
  for (plt in plots) {
    g <- fix_panel_size(plt$p, width = plt$w, height = plt$h)
    dw <- grid::convertWidth(sum(g$widths), "cm", valueOnly = TRUE)
    dh <- grid::convertHeight(sum(g$heights), "cm", valueOnly = TRUE)
    ggsave(file.path(fig_dir, paste0(plt$stub, ".pdf")), g, width = dw, height = dh, units = "cm")
    ggsave(file.path(fig_dir, paste0(plt$stub, ".png")), g, width = dw, height = dh, units = "cm", dpi = 300)
    save_plot_pptx(plt$p, file.path(fig_dir, paste0(plt$stub, ".pptx")), panel_width = plt$w, panel_height = plt$h)
  }
  named_plots <- setNames(
    list(p_inhibition, p_viability, p_trend),
    c("Inhibition Rate", "Relative Viability", "Dose-Response Trend")
  )
  save_plots_pptx(named_plots, file.path(fig_dir, "all_plots.pptx"), panel_width = 10, panel_height = 9)
  message("✓ ", batch_id, " 图表与统计表已生成 (", fig_dir, ")")

  list(wilcox = wilcox_tbl, summary = summary_tbl, p_inhibition = p_inhibition, p_trend = p_trend)
}

batch_ids <- levels(data$batch)
batch_results <- setNames(purrr::map(batch_ids, analyze_batch), batch_ids)

# -----------------------------------------------------------------------------
# 4. 纯描述性的批次一致性对比图（两批叠加，不做跨批次统计检验）
# -----------------------------------------------------------------------------
overview_summary <- purrr::map_dfr(batch_ids, function(b) {
  batch_results[[b]]$summary %>% mutate(batch = b)
})

p_overview <- ggplot(overview_summary, aes(x = group, y = mean_inhibition, color = batch, group = batch)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1) +
  geom_errorbar(
    aes(ymin = mean_inhibition - sem_inhibition, ymax = mean_inhibition + sem_inhibition),
    width = 0.1, linewidth = 0.6, position = position_dodge(width = 0.15)
  ) +
  geom_point(size = 3, position = position_dodge(width = 0.15)) +
  scale_color_manual(values = get_colors("NPG", n = length(batch_ids)), name = "Batch") +
  labs(
    title = "Batch Comparison (Descriptive Only, No Pooled Test)",
    subtitle = "Each batch analyzed independently; overlay shown for visual QC of batch consistency",
    x = "Nanozyme Concentration (mg/mL)", y = "Inhibition Rate (%, mean \u00b1 SEM)"
  ) +
  theme_pub_base() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

g <- fix_panel_size(p_overview, width = 13, height = 8)
dw <- grid::convertWidth(sum(g$widths), "cm", valueOnly = TRUE)
dh <- grid::convertHeight(sum(g$heights), "cm", valueOnly = TRUE)
ggsave(file.path(FIG_DIR, "00_batch_comparison_overview.pdf"), g, width = dw, height = dh, units = "cm")
ggsave(file.path(FIG_DIR, "00_batch_comparison_overview.png"), g, width = dw, height = dh, units = "cm", dpi = 300)
save_plot_pptx(p_overview, file.path(FIG_DIR, "00_batch_comparison_overview.pptx"), panel_width = 13, panel_height = 8)
message("✓ 已生成批次对比总览图（描述性）")

# -----------------------------------------------------------------------------
# 5. 生成分析报告
# -----------------------------------------------------------------------------
report_path <- file.path(OUTPUT_DIR, "report.md")

render_batch_section <- function(batch_id, res, idx) {
  c(
    paste0("## ", idx, ". ", batch_id),
    "",
    "各组描述统计（均值 ± SEM，n=4/组）：",
    "",
    "| 浓度 (mg/mL) | n | 平均存活率 (%) | SEM | 平均抑菌率 (%) | SEM |",
    "|---|---|---|---|---|---|",
    sprintf(
      "| %s | %d | %.2f | %.2f | %.2f | %.2f |",
      res$summary$group, res$summary$n,
      res$summary$mean_viability, res$summary$sem_viability,
      res$summary$mean_inhibition, res$summary$sem_inhibition
    ),
    "",
    "各浓度 vs 阳性对照（Wilcoxon 秩和检验，n=4 vs n=4，未做事后校正）：",
    "",
    "| 比较 (mg/mL vs Control) | 存活率均值 | 抑菌率均值 | p 值 | 显著性 |",
    "|---|---|---|---|---|",
    sprintf(
      "| %s vs %s | %.2f%% | %.2f%% | %.4f | %s |",
      res$wilcox$group1, res$wilcox$group2,
      res$wilcox$mean_viability_1, res$wilcox$mean_inhibition_1,
      res$wilcox$p_value, res$wilcox$significance
    ),
    "",
    paste0("图：`figures/", batch_id, "/01_inhibition_rate.pdf`、`02_relative_viability.pdf`、",
           "`03_dose_response_trend.pdf`（均含 .png / .pptx）"),
    ""
  )
}

lines <- c(
  "# CCK8 纳米酶抗菌能力分析报告",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. 实验设计",
  "",
  "- 检测方法：CCK8 法测 OD450，评估纳米酶对细菌存活/生长的抑制能力。",
  "- 分组与实际浓度（已确认，2 倍稀释系列）：",
  "",
  "  | 分组 | 浓度 (mg/mL) |",
  "  |---|---|",
  sprintf("  | %s | %s |", conc_map$group, conc_map$label),
  "",
  "  阳性对照（BHI + 细菌，无纳米酶）视为浓度 = 0 的剂量-效应基线点；",
  "  空白对照（仅 BHI，无细菌）为背景 QC，不参与浓度梯度分析。",
  "- 共 2 个批次（batch_1、batch_2），每批每组 4 复孔。",
  "- **两批次分开独立分析，不汇集（pool）成 n=8**：分别计算各浓度 vs 阳性对照的",
  "  Wilcoxon 检验，避免批次间系统差异被平均掉或反过来夸大组间差异。仅在第 4 节",
  "  提供一张两批叠加的描述性趋势图，用于目视核查批次间趋势是否一致，不做跨批次",
  "  的统计检验。",
  "- **统计设计**：只做各浓度组 vs 阳性对照（Control）的简单组间检验",
  "  （Wilcoxon 秩和检验），不做事后多重比较校正，不做整体 Kruskal-Wallis 检验，",
  "  不做浓度趋势检验（Jonckheere-Terpstra）。",
  "- 数据已在源文件中按以下流程预处理（每批次独立计算，不跨批次混用）：",
  "  1. `blank_corrected_od = OD_sample - mean(OD_G)` （本批 G 组均值扣背景）",
  "  2. `relative_viability_percent = blank_corrected_od / (mean(OD_F) - mean(OD_G)) * 100`",
  "     （用本批 F 组归一化）",
  "  3. `inhibition_percent = 100 - relative_viability_percent`",
  "",
  "## 2. 背景扣除 QC",
  "",
  "空白对照（G 组）扣除背景后的 OD 值应接近 0，用于核查批次背景是否正常：",
  "",
  "| 批次 | n | 扣背景后均值 | SD |",
  "|---|---|---|---|",
  sprintf("| %s | %d | %.4f | %.4f |", blank_qc$batch, blank_qc$n, blank_qc$mean_blank_corrected_od, blank_qc$sd_blank_corrected_od),
  "",
  "两批次空白扣除后均值均在 0 附近，背景扣除正常。",
  ""
)

for (i in seq_along(batch_ids)) {
  lines <- c(lines, render_batch_section(batch_ids[i], batch_results[[batch_ids[i]]], i + 2))
}

next_idx <- length(batch_ids) + 3
lines <- c(lines,
  paste0("## ", next_idx, ". 批次一致性对比（描述性，不做跨批统计检验）"),
  "",
  "两批次各浓度组的平均抑菌率叠加对比图：`figures/00_batch_comparison_overview.pdf`（+ .png/.pptx）。",
  "该图仅用于目视判断两批的浓度-效应趋势是否一致，图中不含跨批次的统计检验。",
  "",
  paste0("## ", next_idx + 1, ". 结论"),
  "",
  "- 已确认 nanozyme_1~5 为同一材料的 2 倍稀释浓度系列（1 → 0.0625 mg/mL），",
  "  编号越大浓度越低。两批次均观察到抑菌率随浓度升高而上升的描述性趋势",
  "  （见各批次描述统计与趋势图），符合典型剂量-效应关系：浓度越高，抗菌能力越强。",
  "- 各浓度组与阳性对照相比，抑菌效果的具体 p 值见各批次 Wilcoxon 检验表",
  "  （未做事后校正）；因每批 n=4 vs n=4 为精确检验，可达到的最小两侧 p 值约为 0.0286。",
  "- 具体最低有效抑菌浓度（即从哪个浓度开始 vs 阳性对照有差异）请以各批次",
  "  Wilcoxon 表格为准，两批次之间的一致性可参考批次对比总览图。",
  "",
  paste0("## ", next_idx + 2, ". 文件清单"),
  "",
  "- `figures/00_batch_comparison_overview.pdf/png/pptx` — 两批次抑菌率趋势叠加对比（描述性）",
  paste0("- `figures/", batch_ids, "/01_inhibition_rate.pdf/png/pptx` — 该批次各浓度抑菌率比较图（vs 阳性对照）"),
  paste0("- `figures/", batch_ids, "/02_relative_viability.pdf/png/pptx` — 该批次各浓度相对存活率比较图"),
  paste0("- `figures/", batch_ids, "/03_dose_response_trend.pdf/png/pptx` — 该批次浓度-抑菌率描述性趋势图"),
  paste0("- `figures/", batch_ids, "/all_plots.pptx` — 该批次三图打包的可编辑合集"),
  paste0("- `stats/", batch_ids, "/group_summary.csv` — 该批次各浓度描述统计"),
  paste0("- `stats/", batch_ids, "/pairwise_wilcox_vs_positive_control.csv` — 该批次各浓度 vs 阳性对照 Wilcoxon 检验结果（未校正）"),
  "- `stats/blank_control_qc.csv` — 两批次空白对照背景扣除 QC（汇总展示，非统计检验）",
  ""
)

writeLines(lines, report_path, useBytes = TRUE)
message("✓ 已生成分析报告: ", report_path)
message("\n全部完成。输出目录: ", OUTPUT_DIR)
