# =============================================================================
# 生长曲线（OD600）纳米酶抗菌能力分析 —— 单批次
#
# 用法：
#   1) 改下面的 BATCH_ID，然后整个 source 运行；或者
#   2) 命令行：Rscript 02_analyze_batch.R <batch_id>（参数会覆盖 BATCH_ID）
# 前提：该批次已经跑过 01_prepare_data.R（即 data/<BATCH_ID>/ 下已有三张表）。
#
# 只处理单个批次的出图 + 统计，不生成跨批次报告。所有批次都跑完 02 后，运行
# 03_compare_batches.R 生成批次对比图 + 一致性核查 + 汇总 report.md
# （03 会自动扫描 data/ 下已处理好的批次，不需要在脚本里手动列批次名单）。
#
# 统计设计（与项目既有惯例一致，见 r-pub-toolkit 的 PLOTTING_CONVENTIONS.md）：
#   - 各浓度组 vs 阳性对照（两独立样本 t 检验，Welch，不假设方差齐性），逐孔 n=3 vs n=3
#   - 报告原始 p 值 + BH 校正 p 值；图上星号为未校正的探索性展示
#   - 注意：n=3/组的小样本 t 检验对正态性假设很敏感，检验效能也低，结果仍应以
#     趋势参考为主，不宜仅凭 p<0.05 下达"显著/不显著"的结论。
# =============================================================================

# 定位仓库根目录（本地/服务器通用）：优先 R_TOOLKIT_ROOT 环境变量，否则用 git
PROJ <- Sys.getenv("R_TOOLKIT_ROOT", unset = "")
if (!nzchar(PROJ)) PROJ <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE), error = function(e) "")
if (length(PROJ) != 1 || !nzchar(PROJ)) stop("找不到仓库根目录：请设置 R_TOOLKIT_ROOT 环境变量，或在仓库目录内运行")
source(file.path(PROJ, "init.R"))
source(file.path(PROJ, "collaborators/Gilly/Growth_curve/SM/scripts/_common.R"))

BATCH_ID <- "batch3" # <-- 手动改成要处理的批次，例如 "batch2" / "batch3"
BATCH_ID <- resolve_batch_id(BATCH_ID)

data_dir <- file.path(DATA_DIR, BATCH_ID)
fig_dir <- file.path(FIG_DIR, BATCH_ID)
stat_dir <- file.path(STAT_DIR, BATCH_ID)

stopifnot(
  "data/<BATCH_ID>/ 不存在，请先跑 01_prepare_data.R（同一个 BATCH_ID）" = dir.exists(
    data_dir
  )
)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stat_dir, recursive = TRUE, showWarnings = FALSE)

summary_data <- read.csv(
  file.path(data_dir, "growth_curve_summary.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(group = factor(group, levels = conc_map$label))
replicate_metrics <- read.csv(
  file.path(data_dir, "replicate_endpoint_auc.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(group = factor(group, levels = conc_map$label))

message("\n==== 分析 ", BATCH_ID, " ====")

# --- 1. 生长曲线主图：raw 与 blank-corrected ---
p_raw <- plot_growth_curve(
  summary_data,
  "od600_mean",
  "od600_sd",
  "OD600 (raw)",
  paste0("Growth Curve (Raw OD600) - ", BATCH_ID)
)
p_blank <- plot_growth_curve(
  summary_data,
  "od600_blank_corrected_mean",
  "od600_blank_corrected_sd",
  "OD600 (blank-corrected)",
  paste0("Growth Curve (Blank-Corrected OD600) - ", BATCH_ID)
)
save_fig(p_raw, fig_dir, "00_growth_curve_raw", 13, 9)
save_fig(p_blank, fig_dir, "01_growth_curve_blank_corrected", 13, 9)

# --- 2/3. 组间比较：终点OD600、AUC，逐孔 n=3 vs 阳性对照 ---
p_endpoint <- build_stat_plot(
  replicate_metrics,
  "final_od600_blank_corrected",
  "Final OD600 (blank-corrected, 13.5h)",
  paste0("Endpoint Growth (13.5h) - ", BATCH_ID)
)
p_auc <- build_stat_plot(
  replicate_metrics,
  "auc_od600_blank_corrected",
  "AUC (blank-corrected OD600 \u00d7 h)",
  paste0("Growth AUC (0-13.5h) - ", BATCH_ID)
)
save_fig(p_endpoint, fig_dir, "02_endpoint_od600_comparison", 11, 9.5)
save_fig(p_auc, fig_dir, "03_auc_comparison", 11, 9.5)

# --- 4. 统计表：终点OD600、AUC 的 t 检验（Welch，原始p + BH校正p） ---
ttest_tbl <- bind_rows(
  run_ttest_table(replicate_metrics, "final_od600_blank_corrected"),
  run_ttest_table(replicate_metrics, "auc_od600_blank_corrected")
)
write.csv(
  ttest_tbl,
  file.path(stat_dir, "pairwise_ttest_vs_control.csv"),
  row.names = FALSE
)

# --- 5. 描述统计：各组终点/AUC 均值 ± SEM ---
group_summary <- replicate_metrics %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean_final_od600 = mean(final_od600_blank_corrected),
    sem_final_od600 = sd(final_od600_blank_corrected) / sqrt(n()),
    mean_auc = mean(auc_od600_blank_corrected),
    sem_auc = sd(auc_od600_blank_corrected) / sqrt(n()),
    .groups = "drop"
  ) %>%
  left_join(conc_map, by = c("group" = "label"))
write.csv(
  group_summary,
  file.path(stat_dir, "group_summary.csv"),
  row.names = FALSE
)

# --- 6. 剂量-效应描述性趋势图：终点生长抑制率（相对阳性对照 Control 组均值） ---
control_mean_final <- group_summary$mean_final_od600[
  group_summary$group == control_label
]
inhibition_summary <- replicate_metrics %>%
  mutate(
    inhibition_percent = 100 *
      (1 - final_od600_blank_corrected / control_mean_final)
  ) %>%
  group_by(group, conc_ugml) %>%
  summarise(
    n = n(),
    mean_inhibition = mean(inhibition_percent),
    sem_inhibition = sd(inhibition_percent) / sqrt(n()),
    .groups = "drop"
  ) %>%
  arrange(conc_ugml)
write.csv(
  inhibition_summary,
  file.path(stat_dir, "inhibition_percent_summary.csv"),
  row.names = FALSE
)

p_trend <- ggplot(
  inhibition_summary,
  aes(x = group, y = mean_inhibition, group = 1)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(color = get_colors("NPG", n = 1), linewidth = 1) +
  geom_errorbar(
    aes(
      ymin = mean_inhibition - sem_inhibition,
      ymax = mean_inhibition + sem_inhibition
    ),
    width = 0.1,
    linewidth = 0.6
  ) +
  geom_point(
    size = 3.5,
    shape = 21,
    fill = get_colors("NPG", n = 1),
    color = "black",
    stroke = 0.8
  ) +
  labs(
    title = paste0("Nanozyme Dose-Response (Growth Inhibition) - ", BATCH_ID),
    subtitle = "Descriptive only, mean \u00b1 SEM, endpoint (13.5h) vs Control",
    x = "Nanozyme Concentration (ug/mL)",
    y = "Growth Inhibition (%, mean \u00b1 SEM)"
  ) +
  theme_pub_base() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_trend, fig_dir, "04_dose_response_trend", 12, 8)

named_plots <- setNames(
  list(p_raw, p_blank, p_endpoint, p_auc, p_trend),
  c(
    "Growth Curve (Raw)",
    "Growth Curve (Blank-Corrected)",
    "Endpoint OD600",
    "AUC",
    "Dose-Response Trend"
  )
)
save_plots_pptx(
  named_plots,
  file.path(fig_dir, "all_plots.pptx"),
  panel_width = 12,
  panel_height = 8.5
)

message("✓ [", BATCH_ID, "] 图表与统计表已生成")
message("  - figures/", BATCH_ID, "/  (5 张图 x pdf/png/pptx + all_plots.pptx)")
message(
  "  - stats/",
  BATCH_ID,
  "/  (pairwise_ttest_vs_control.csv, group_summary.csv, inhibition_percent_summary.csv)"
)
message(
  "全部批次都跑完后，运行 03_compare_batches.R 生成批次对比图与汇总报告。"
)
