# =============================================================================
# CCK8 纳米酶抗菌能力分析：Pg / Fn 两种细菌数据
# 数据：collaborators/Gilly/CCK8/data/CCK8_Pg_Fn_without_A_tidy.csv
# 背景：单板实验，每组 3 复孔，不含批次重复（与 cck8_antibacterial_analysis.R
#       分析的另一份 2 批次数据是独立实验，不合并）。nanozyme_1（A 组，1 mg/mL）
#       在源数据中已被上游剔除（note 列注明"因疑似浓度/设置问题移除"）。
#       浓度映射（与另一份数据一致的 2 倍稀释系列，只是缺 1 mg/mL 这一级）：
#         nanozyme_2 = 0.5, nanozyme_3 = 0.25, nanozyme_4 = 0.125,
#         nanozyme_5 = 0.0625 mg/mL；positive_control 视为浓度 = 0。
#       Pg 和 Fn 是两种不同细菌，浓度-效应关系不能混在一起算，必须分开——
#       本脚本按 bacteria 列分别处理。
#
# Pg 数据质量问题（已与用户确认）：
#   nanozyme_3 复孔2、nanozyme_4 复孔2/3、nanozyme_5 复孔3 共 4 个孔的 OD450
#   是该批次阳性对照均值的 2.3~2.9 倍，纳米酶处理组的活菌量理论上不应超过
#   "无处理"阳性对照，判断为技术性错误（串孔污染/读数错误等）。用户决定：
#   本次不对 Pg 出正式结果，只保留问题记录，Pg 需要重新做实验。
#   本脚本因此只对 Fn 做完整分析，Pg 仅生成一份异常孔记录表。
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

dir.create(file.path(FIG_DIR, "fn"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(STAT_DIR, "fn"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(STAT_DIR, "pg"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. 读取数据
# -----------------------------------------------------------------------------
raw <- read.csv(
  file.path(DATA_DIR, "CCK8_Pg_Fn_without_A_tidy.csv"),
  stringsAsFactors = FALSE
)

conc_map <- tibble::tribble(
  ~group,              ~conc_mgml, ~label,
  "positive_control",  0,          "0 (Control)",
  "nanozyme_5",        0.0625,     "0.0625",
  "nanozyme_4",        0.125,      "0.125",
  "nanozyme_3",        0.25,       "0.25",
  "nanozyme_2",        0.5,        "0.5"
)

# -----------------------------------------------------------------------------
# 2. Pg：仅记录异常孔，不出正式统计结果（用户已确认这批 Pg 数据需重做）
# -----------------------------------------------------------------------------
pg_raw <- raw %>% filter(bacteria == "Pg", group %in% conc_map$group)

pg_flag <- pg_raw %>%
  mutate(
    od_to_positive_ratio = od450 / positive_mean,
    flagged_implausible = od_to_positive_ratio > 1.5
  ) %>%
  select(bacteria, group, replicate, well, od450, positive_mean, od_to_positive_ratio, flagged_implausible) %>%
  arrange(desc(od_to_positive_ratio))

write.csv(pg_flag, file.path(STAT_DIR, "pg", "pg_data_quality_flag.csv"), row.names = FALSE)
message("Pg 数据质量记录已生成（本次不出正式分析结果，等待重做实验）：")
print(pg_flag %>% filter(flagged_implausible))

# -----------------------------------------------------------------------------
# 3. Fn：完整分析（背景 QC + 描述统计 + vs 阳性对照 Wilcoxon + 绘图）
#    统计设计（按要求简化）：只做各浓度组 vs 阳性对照的简单组间检验，不做事后
#    多重比较校正，不做整体 Kruskal-Wallis 检验，不做浓度趋势检验。
# -----------------------------------------------------------------------------
fn_data <- raw %>%
  filter(bacteria == "Fn", group %in% conc_map$group) %>%
  left_join(conc_map, by = "group") %>%
  mutate(group = factor(group, levels = conc_map$group, labels = conc_map$label))

message("\nFn 各组样本量（n=3/组，单板实验，无批次重复）：")
print(table(fn_data$group))

fn_blank_qc <- raw %>%
  filter(bacteria == "Fn", group == "blank_control") %>%
  summarise(
    n = n(),
    mean_blank_corrected_od = mean(blank_corrected_od),
    sd_blank_corrected_od = sd(blank_corrected_od)
  )
write.csv(fn_blank_qc, file.path(STAT_DIR, "fn", "blank_control_qc.csv"), row.names = FALSE)

positive_label <- conc_map$label[conc_map$group == "positive_control"]
dose_pairs <- purrr::map(setdiff(conc_map$label, positive_label), function(l) c(l, positive_label))

wilcox_tbl <- purrr::map_dfr(dose_pairs, function(pr) {
  d1 <- fn_data$relative_viability_percent[fn_data$group == pr[1]]
  d2 <- fn_data$relative_viability_percent[fn_data$group == pr[2]]
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
write.csv(wilcox_tbl, file.path(STAT_DIR, "fn", "pairwise_wilcox_vs_positive_control.csv"), row.names = FALSE)

summary_tbl <- fn_data %>%
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
write.csv(summary_tbl, file.path(STAT_DIR, "fn", "group_summary.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. 绘图
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
      title = "CCK8 Antibacterial Assay (Fn)", xlab = "Nanozyme Concentration (mg/mL)", ylab = y_label
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

# 与 nanozyme_1~5 数据集同理：inhibition_percent 由 viability 反算，噪声可能
# 使部分点略低于 0%，硬顶在 100 会裁掉堆叠比较括号里最上面的标注，故 cap = Inf。
p_inhibition <- build_stat_plot(fn_data, "inhibition_percent", "Inhibition Rate (%)", dose_pairs, cap = Inf) +
  labs(subtitle = "Wilcoxon rank-sum vs Control, unadjusted p (n=3/group, Fn, single plate)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_viability <- build_stat_plot(fn_data, "relative_viability_percent", "Relative Viability (%)", dose_pairs, cap = Inf) +
  labs(subtitle = "Wilcoxon rank-sum vs Control, unadjusted p (n=3/group, Fn, single plate)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_trend <- ggplot(summary_tbl, aes(x = group, y = mean_inhibition, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(color = get_colors("NPG", n = 1), linewidth = 1) +
  geom_errorbar(
    aes(ymin = mean_inhibition - sem_inhibition, ymax = mean_inhibition + sem_inhibition),
    width = 0.1, linewidth = 0.6
  ) +
  geom_point(size = 3.5, shape = 21, fill = get_colors("NPG", n = 1), color = "black", stroke = 0.8) +
  labs(
    title = "Nanozyme Dose-Response (Fn Inhibition)",
    subtitle = "Descriptive only, mean \u00b1 SEM (n=3/group)",
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
  ggsave(file.path(FIG_DIR, "fn", paste0(plt$stub, ".pdf")), g, width = dw, height = dh, units = "cm")
  ggsave(file.path(FIG_DIR, "fn", paste0(plt$stub, ".png")), g, width = dw, height = dh, units = "cm", dpi = 300)
  save_plot_pptx(plt$p, file.path(FIG_DIR, "fn", paste0(plt$stub, ".pptx")), panel_width = plt$w, panel_height = plt$h)
}
named_plots <- setNames(
  list(p_inhibition, p_viability, p_trend),
  c("Inhibition Rate", "Relative Viability", "Dose-Response Trend")
)
save_plots_pptx(named_plots, file.path(FIG_DIR, "fn", "all_plots.pptx"), panel_width = 10, panel_height = 9)
message("✓ Fn 图表与统计表已生成 (", file.path(FIG_DIR, "fn"), ")")

# -----------------------------------------------------------------------------
# 5. 生成分析报告（独立文件，与 nanozyme_1~5 那份 2 批次数据的报告分开）
# -----------------------------------------------------------------------------
report_path <- file.path(OUTPUT_DIR, "report_fn_pg.md")

pg_flagged <- pg_flag %>% filter(flagged_implausible)

lines <- c(
  "# CCK8 纳米酶抗菌能力分析报告：Pg / Fn 细菌特异性数据",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. 实验设计",
  "",
  "- 数据来源：`data/CCK8_Pg_Fn_without_A_tidy.csv`，与另一份 `report.md` 分析的",
  "  2 批次 nanozyme_1~5 数据是**独立实验**，本报告不与其合并。",
  "- 单板实验，每组 3 复孔，无批次重复。nanozyme_1（A 组，1 mg/mL）已在源数据中",
  "  被上游剔除（\"因疑似浓度/设置问题移除\"），本次浓度梯度只覆盖 0.5/0.25/0.125/0.0625 mg/mL",
  "  + 阳性对照（视为浓度 = 0）。",
  "- **Pg（牙龈卟啉单胞菌）和 Fn（具核梭杆菌）是两种不同细菌，浓度-效应关系分开计算，",
  "  不汇集**。",
  "",
  "## 2. Pg 数据质量问题（本次不出正式结果）",
  "",
  "**⚠️ Pg 组数据存在明显的技术性异常，已与用户确认：本次不对 Pg 做正式统计分析，",
  "该批 Pg 实验建议重做。**",
  "",
  sprintf(
    "以下 %d 个孔的 OD450 是当批阳性对照均值（%.4f）的 1.5 倍以上，纳米酶处理组的",
    nrow(pg_flagged), pg_flag$positive_mean[1]
  ),
  "活菌量理论上不应超过\"无处理\"阳性对照，判断为串孔污染、读数错误等技术性问题：",
  "",
  "| 分组 | 复孔 | 孔位 | OD450 | 阳性对照均值 | OD/阳性对照 倍数 |",
  "|---|---|---|---|---|---|",
  sprintf(
    "| %s | %d | %s | %.4f | %.4f | %.2f× |",
    pg_flagged$group, pg_flagged$replicate, pg_flagged$well,
    pg_flagged$od450, pg_flagged$positive_mean, pg_flagged$od_to_positive_ratio
  ),
  "",
  "完整的 Pg 原始数据 + 标记结果见 `stats/pg/pg_data_quality_flag.csv`",
  "（`flagged_implausible` 列标注是否超过 1.5 倍阈值），供核对原始读数或重做实验时参考。",
  "",
  "## 3. Fn 背景扣除 QC",
  "",
  "空白对照扣除背景后的 OD 值应接近 0：",
  "",
  "| n | 扣背景后均值 | SD |",
  "|---|---|---|",
  sprintf("| %d | %.4f | %.4f |", fn_blank_qc$n, fn_blank_qc$mean_blank_corrected_od, fn_blank_qc$sd_blank_corrected_od),
  "",
  "背景扣除正常。",
  "",
  "## 4. Fn 各组描述统计",
  "",
  "| 浓度 (mg/mL) | n | 平均存活率 (%) | SEM | 平均抑菌率 (%) | SEM |",
  "|---|---|---|---|---|---|",
  sprintf(
    "| %s | %d | %.2f | %.2f | %.2f | %.2f |",
    summary_tbl$group, summary_tbl$n,
    summary_tbl$mean_viability, summary_tbl$sem_viability,
    summary_tbl$mean_inhibition, summary_tbl$sem_inhibition
  ),
  "",
  "## 5. Fn 各浓度 vs 阳性对照（Wilcoxon 秩和检验，n=3 vs n=3，未做事后校正）",
  "",
  "**注意**：n=3 vs n=3 的 Wilcoxon 精确检验，两侧最小可达 p 值约为 0.1",
  "（完全分离时 2×C(6,3)⁻¹ = 0.1），达不到传统 0.05 显著性阈值属检验力不足，",
  "不代表处理无效，解读时需结合效应量（均值差异）一并判断。",
  "",
  "| 比较 (mg/mL vs Control) | 存活率均值 | 抑菌率均值 | p 值 | 显著性 |",
  "|---|---|---|---|---|",
  sprintf(
    "| %s vs %s | %.2f%% | %.2f%% | %.4f | %s |",
    wilcox_tbl$group1, wilcox_tbl$group2,
    wilcox_tbl$mean_viability_1, wilcox_tbl$mean_inhibition_1,
    wilcox_tbl$p_value, wilcox_tbl$significance
  ),
  "",
  "## 6. 结论",
  "",
  "- **Pg**：数据存在明显技术性异常（4/12 处理孔 OD450 超阳性对照 1.5 倍以上），",
  "  本次不出正式抗菌效果结论，建议重做实验后再分析。",
  "- **Fn**：与另一份 nanozyme_1~5 数据集里观察到的干净单调剂量-效应关系不同，",
  "  本次 Fn 的各浓度组抑菌率均值依次为 0.0625→约 80%、0.125→约 41%、0.25→约 90%、",
  "  0.5→约 50%，呈现\"高-低-高-低\"的振荡模式（非单调）：0.25 mg/mL（row C）三个复孔",
  "  的 OD450 均接近空白对照水平、明显低于相邻的 0.5 mg/mL（row B）和 0.125 mg/mL",
  "  （row D），即抑菌效果反而比更高浓度的 0.5 mg/mL 组更强。这不属于第 2 节那种",
  "  \"单孔读数超过阳性对照\"的物理不可能情形，因此本报告仍保留全部数据出正式结果，",
  "  但建议核对 row C（nanozyme_3, 0.25 mg/mL）与 row B（nanozyme_2, 0.5 mg/mL）",
  "  的原始配液/加样记录，确认是否存在浓度配制或加样顺序上的问题，",
  "  否则这批 Fn 数据不能得出\"浓度越高抑菌越强\"的结论。",
  "- 本次 Fn 实验只有单板 n=3、且缺少最高浓度 1 mg/mL 数据点，各浓度 vs 阳性对照",
  "  两两比较由于样本量小、检验力有限（Wilcoxon 精确检验两侧最小 p≈0.1），",
  "  第 5 节里全部比较均为 ns，解读时需结合效应量（均值抑菌率）而非单纯看 p 值。",
  "",
  "## 7. 文件清单",
  "",
  "- `figures/fn/01_inhibition_rate.pdf/png/pptx` — Fn 各浓度抑菌率比较图（vs 阳性对照）",
  "- `figures/fn/02_relative_viability.pdf/png/pptx` — Fn 各浓度相对存活率比较图",
  "- `figures/fn/03_dose_response_trend.pdf/png/pptx` — Fn 浓度-抑菌率描述性趋势图",
  "- `figures/fn/all_plots.pptx` — Fn 三图打包的可编辑合集",
  "- `stats/fn/group_summary.csv` — Fn 各浓度描述统计",
  "- `stats/fn/pairwise_wilcox_vs_positive_control.csv` — Fn 各浓度 vs 阳性对照 Wilcoxon 检验结果（未校正）",
  "- `stats/fn/blank_control_qc.csv` — Fn 空白对照背景扣除 QC",
  "- `stats/pg/pg_data_quality_flag.csv` — Pg 全部原始数据 + 异常孔标记（未出正式结果）",
  ""
)

writeLines(lines, report_path, useBytes = TRUE)
message("✓ 已生成分析报告: ", report_path)
message("\n全部完成。输出目录: ", OUTPUT_DIR)
