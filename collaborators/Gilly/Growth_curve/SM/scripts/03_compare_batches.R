# =============================================================================
# 生长曲线（OD600）纳米酶抗菌能力分析 —— 批次对比 + 汇总报告
#
# 用法：不需要改脚本、不需要传参数，直接 source 或 Rscript 运行即可。
# 自动扫描 data/ 目录下"已经跑过 02_analyze_batch.R"的批次文件夹（即
# output/stats/<batch>/ 下三张统计表齐全的批次），生成：
#   - 批次抑制率趋势叠加对比图（纯描述性，不做跨批次统计检验，见
#     [[project-instructions]] 里两批次不汇集的惯例，与 CCK8 项目一致）
#   - 浓度-抑制率 Spearman 相关的批次一致性核查
#   - 汇总 report.md（每个批次一节 + 批次对比小节 + 文件清单）
#
# 新增批次时不需要改这个脚本：只要该批次跑完了 01/02，这里重新跑一次就会自动
# 把它纳入对比和报告。
# =============================================================================

# 定位仓库根目录（本地/服务器通用）：优先 R_TOOLKIT_ROOT 环境变量，否则用 git
PROJ <- Sys.getenv("R_TOOLKIT_ROOT", unset = "")
if (!nzchar(PROJ)) PROJ <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE), error = function(e) "")
if (length(PROJ) != 1 || !nzchar(PROJ)) stop("找不到仓库根目录：请设置 R_TOOLKIT_ROOT 环境变量，或在仓库目录内运行")
source(file.path(PROJ, "init.R"))

# 本仓库（olp-microbiome-analysis）根目录：脚本本身在这个仓库里，用 git 定位
# （不能用 PROJ——PROJ 现在指向独立的 r-pub-toolkit 仓库，没有 collaborators/ 目录）
REPO_ROOT <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE), error = function(e) "")
if (length(REPO_ROOT) != 1 || !nzchar(REPO_ROOT)) stop("找不到本仓库根目录：请在仓库目录内运行脚本")
source(file.path(REPO_ROOT, "collaborators/Gilly/Growth_curve/SM/scripts/_common.R"))

# -----------------------------------------------------------------------------
# 自动发现已处理好的批次：data/<batch>/ 存在，且 output/stats/<batch>/ 下
# 三张统计表齐全（说明 02_analyze_batch.R 已经跑过）
# -----------------------------------------------------------------------------
candidate_batches <- list.dirs(DATA_DIR, full.names = FALSE, recursive = FALSE)

required_stat_files <- c(
  "group_summary.csv",
  "pairwise_ttest_vs_control.csv",
  "inhibition_percent_summary.csv"
)
is_ready <- purrr::map_lgl(candidate_batches, function(b) {
  all(file.exists(file.path(STAT_DIR, b, required_stat_files)))
})
batch_ids <- candidate_batches[is_ready]

if (length(candidate_batches) > length(batch_ids)) {
  skipped <- setdiff(candidate_batches, batch_ids)
  message(
    "跳过尚未跑 02_analyze_batch.R 的批次: ",
    paste(skipped, collapse = ", ")
  )
}
stopifnot(
  "没有找到任何已完成 02_analyze_batch.R 的批次，请先跑 01 + 02" = length(
    batch_ids
  ) >
    0
)
message("纳入对比的批次: ", paste(batch_ids, collapse = ", "))

read_batch_stats <- function(batch_id) {
  stat_dir <- file.path(STAT_DIR, batch_id)
  list(
    group_summary = read.csv(
      file.path(stat_dir, "group_summary.csv"),
      stringsAsFactors = FALSE
    ) %>%
      mutate(group = factor(group, levels = conc_map$label)),
    ttest_tbl = read.csv(
      file.path(stat_dir, "pairwise_ttest_vs_control.csv"),
      stringsAsFactors = FALSE
    ),
    inhibition_summary = read.csv(
      file.path(stat_dir, "inhibition_percent_summary.csv"),
      stringsAsFactors = FALSE
    ) %>%
      mutate(group = factor(group, levels = conc_map$label))
  )
}
batch_results <- setNames(purrr::map(batch_ids, read_batch_stats), batch_ids)

# -----------------------------------------------------------------------------
# 批次一致性对比：各批次抑制率趋势叠加（纯描述性，不做跨批次统计检验）
# -----------------------------------------------------------------------------
overview_summary <- purrr::map_dfr(batch_ids, function(b) {
  batch_results[[b]]$inhibition_summary %>% mutate(batch = b)
})

p_overview <- ggplot(
  overview_summary,
  aes(x = group, y = mean_inhibition, color = batch, group = batch)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1) +
  geom_errorbar(
    aes(
      ymin = mean_inhibition - sem_inhibition,
      ymax = mean_inhibition + sem_inhibition
    ),
    width = 0.1,
    linewidth = 0.6,
    position = position_dodge(width = 0.15)
  ) +
  geom_point(size = 3, position = position_dodge(width = 0.15)) +
  scale_color_manual(
    values = get_colors("NPG", n = length(batch_ids)),
    name = "Batch"
  ) +
  labs(
    title = "Batch Comparison (Descriptive Only, No Pooled Test)",
    subtitle = "Each batch analyzed independently; overlay shown for visual QC of batch consistency",
    x = "Nanozyme Concentration (ug/mL)",
    y = "Growth Inhibition (%, mean \u00b1 SEM)"
  ) +
  theme_pub_base() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_overview, FIG_DIR, "00_batch_comparison_overview", 13, 8)
message("✓ 已生成批次对比总览图（描述性）")

# --- 批次一致性定量核查：浓度与抑制率的 Spearman 相关（仅 B~F 剂量组，不含 Control） ---
batch_consistency <- overview_summary %>%
  filter(group != control_label) %>%
  group_by(batch) %>%
  summarise(
    spearman_rho = suppressWarnings(cor(
      conc_ugml,
      mean_inhibition,
      method = "spearman"
    )),
    inhibition_at_min_conc = mean_inhibition[which.min(conc_ugml)],
    inhibition_at_max_conc = mean_inhibition[which.max(conc_ugml)],
    .groups = "drop"
  )
write.csv(
  batch_consistency,
  file.path(STAT_DIR, "batch_consistency_check.csv"),
  row.names = FALSE
)

rho_range <- range(batch_consistency$spearman_rho)
consistency_note <- if (
  diff(sign(rho_range)) != 0 || (rho_range[1] < 0 && rho_range[2] >= 0)
) {
  paste(
    "**核查发现批次间剂量-效应方向不一致**（有的批次 rho 为正、有的为负，见下表）。",
    "这是数据本身呈现的现象，不是绘图/统计问题。建议核实各批次孔位与浓度的对应",
    "关系是否一致（会不会摆盘/标签顺序有出入）；若确认摆盘无误，需要额外实验",
    "（如 DLS 粒径、更细的浓度梯度复测）才能解释，不能仅凭本次结果下结论。"
  )
} else {
  "各批次浓度-抑制率的相关方向基本一致（rho 同号），未发现明显的批次间矛盾趋势。"
}

# -----------------------------------------------------------------------------
# 汇总报告（按批次分节 + 批次对比小节）
# -----------------------------------------------------------------------------
report_path <- file.path(OUTPUT_DIR, "report.md")

render_batch_section <- function(batch_id, res, idx) {
  gs <- res$group_summary
  tt <- res$ttest_tbl
  c(
    paste0("## ", idx, ". ", batch_id),
    "",
    "各组描述统计（终点13.5h blank-corrected OD600 与 AUC，mean \u00b1 SEM，n=3）：",
    "",
    "| 组别 (ug/mL) | n | 终点OD600均值 | SEM | AUC均值 | SEM |",
    "|---|---|---|---|---|---|",
    sprintf(
      "| %s | %d | %.4f | %.4f | %.3f | %.3f |",
      gs$group,
      gs$n,
      gs$mean_final_od600,
      gs$sem_final_od600,
      gs$mean_auc,
      gs$sem_auc
    ),
    "",
    "各浓度 vs 阳性对照 Welch t 检验（原始p + BH校正p）：",
    "",
    "| 指标 | 比较 (ug/mL vs Control) | mean1 | mean2 | p 值 | BH校正p | 显著性 |",
    "|---|---|---|---|---|---|---|",
    sprintf(
      "| %s | %s vs %s | %.4f | %.4f | %.4f | %.4f | %s |",
      tt$metric,
      tt$group1,
      tt$group2,
      tt$mean1,
      tt$mean2,
      tt$p_value,
      tt$p_adj_BH,
      tt$significance
    ),
    "",
    paste0(
      "图：`figures/",
      batch_id,
      "/00_growth_curve_raw.pdf`、`01_growth_curve_blank_corrected.pdf`、",
      "`02_endpoint_od600_comparison.pdf`、`03_auc_comparison.pdf`、`04_dose_response_trend.pdf`",
      "（均含 .png / .pptx）"
    ),
    ""
  )
}

lines <- c(
  "# 生长曲线（OD600）纳米酶抗菌能力分析报告",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. 实验设计",
  "",
  "- 检测方法：96 孔板 OD600 生长曲线，每 30 min 读数一次，共 28 个时间点（0~13.5 h）。",
  "- 分组与浓度（用户口头确认，2 倍稀释系列，所有批次共用同一套浓度体系）：",
  "",
  "  | 组别 | 浓度 (ug/mL) |",
  "  |---|---|",
  sprintf("  | %s | %s |", conc_map$group, conc_map$label),
  "",
  "  A 组为阳性对照（细菌 + 无纳米酶处理）；G01-G03 为空白孔（仅培养基，无细菌），",
  "  用于背景扣除，不参与浓度梯度分析。每组 3 复孔（01-03）。",
  "- **重要提示**：浓度映射基于用户口头确认，raw 原始数据文件本身不含浓度标签；",
  "  若与仪器/实验记录有出入，请核对 `scripts/_common.R` 里的 conc_map 并重跑脚本。",
  "- 背景扣除：`od600_blank_corrected = od600 - mean(G组同时间点OD600)`。",
  paste0(
    "- **本次纳入 ",
    length(batch_ids),
    " 个批次**：",
    paste(batch_ids, collapse = "、"),
    "。每个批次的原始数据来源/预处理方式见各自 `raw/<batch>/` 下的说明文件（如有）。"
  ),
  "- **各批次分开独立分析，不汇集（pool）**：分别计算各浓度 vs 阳性对照的 t 检验，",
  "  避免批次间系统差异被平均掉或反过来夸大组间差异（与 collaborators/Gilly/CCK8",
  "  的批次处理方式一致）。仅在最后提供一张各批次叠加的描述性趋势图，用于目视",
  "  核查批次间趋势是否一致，不做跨批次的统计检验。",
  "",
  "## 2. 统计设计与局限",
  "",
  "- 组间比较：各浓度组 vs 阳性对照（两独立样本 t 检验，Welch，不假设方差齐性），",
  "  逐孔 n=3 vs n=3。",
  "- 统计表同时给出原始 p 值和 BH 校正 p 值；图上显著性星号为**未校正**的探索性展示。",
  "- **关键局限**：每组仅 n=3，t 检验对正态性假设敏感、检验效能也低，结果仍应以",
  "  剂量-效应趋势参考为主，不宜仅凭单次 p<0.05 下达\"显著/不显著\"的结论。",
  "- OD600 同时受细菌生长、细胞碎片、材料浊度、纳米颗粒聚集/沉降等因素影响；",
  "  对于纳米酶抗菌实验，OD600 应作为生长抑制的证据，而非直接杀菌证据。若需要",
  "  杀菌结论，应结合 CFU 计数、活死染色、SEM/TEM 膜损伤或 ROS 实验佐证",
  "  （见 `collaborators/Gilly/CFU_killing_assay`）。",
  ""
)

for (i in seq_along(batch_ids)) {
  lines <- c(
    lines,
    render_batch_section(batch_ids[i], batch_results[[batch_ids[i]]], i + 2)
  )
}

next_idx <- length(batch_ids) + 3
lines <- c(
  lines,
  paste0("## ", next_idx, ". 批次一致性对比（描述性，不做跨批统计检验）"),
  "",
  "各批次的平均生长抑制率叠加对比图：`figures/00_batch_comparison_overview.pdf`",
  "（+ .png/.pptx）。该图仅用于目视判断各批次的剂量-效应趋势是否一致，图中不含",
  "跨批次的统计检验。",
  "",
  "浓度-抑制率 Spearman 相关核查（`stats/batch_consistency_check.csv`，仅计算",
  "B~F 剂量组，不含 Control）：",
  "",
  "| 批次 | 浓度-抑制率 Spearman rho | 最低浓度抑制率 | 最高浓度抑制率 |",
  "|---|---|---|---|",
  sprintf(
    "| %s | %.2f | %.1f%% | %.1f%% |",
    batch_consistency$batch,
    batch_consistency$spearman_rho,
    batch_consistency$inhibition_at_min_conc,
    batch_consistency$inhibition_at_max_conc
  ),
  "",
  consistency_note,
  "",
  paste0("## ", next_idx + 1, ". 文件清单"),
  "",
  "- `figures/00_batch_comparison_overview.pdf/png/pptx` — 各批次抑制率趋势叠加对比（描述性）",
  "- `stats/batch_consistency_check.csv` — 各批次浓度-抑制率 Spearman 相关核查",
  paste0(
    "- `figures/",
    batch_ids,
    "/00_growth_curve_raw.pdf/png/pptx` — 该批次原始OD600生长曲线"
  ),
  paste0(
    "- `figures/",
    batch_ids,
    "/01_growth_curve_blank_corrected.pdf/png/pptx` — 该批次扣背景OD600生长曲线（主图）"
  ),
  paste0(
    "- `figures/",
    batch_ids,
    "/02_endpoint_od600_comparison.pdf/png/pptx` — 该批次终点(13.5h)OD600组间比较（vs Control）"
  ),
  paste0(
    "- `figures/",
    batch_ids,
    "/03_auc_comparison.pdf/png/pptx` — 该批次AUC组间比较（vs Control）"
  ),
  paste0(
    "- `figures/",
    batch_ids,
    "/04_dose_response_trend.pdf/png/pptx` — 该批次剂量-效应描述性趋势图"
  ),
  paste0(
    "- `figures/",
    batch_ids,
    "/all_plots.pptx` — 该批次五图打包的可编辑合集"
  ),
  paste0(
    "- `stats/",
    batch_ids,
    "/pairwise_ttest_vs_control.csv` — 该批次终点OD600与AUC的 Welch t 检验（原始p+BH校正p）"
  ),
  paste0(
    "- `stats/",
    batch_ids,
    "/group_summary.csv` — 该批次各组终点OD600/AUC描述统计"
  ),
  paste0(
    "- `stats/",
    batch_ids,
    "/inhibition_percent_summary.csv` — 该批次各组生长抑制率描述统计"
  ),
  "- `data/concentration_map.csv` — 组别-浓度映射表（如浓度需更正，改 `scripts/_common.R` 后重跑即可）",
  paste0(
    "- `data/",
    batch_ids,
    "/replicate_endpoint_auc.csv` — 该批次逐孔终点OD600与AUC（用于统计检验的原始明细）"
  ),
  ""
)

writeLines(lines, report_path, useBytes = TRUE)
message("✓ 已生成分析报告: ", report_path)
message("\n全部完成。输出目录: ", OUTPUT_DIR)
