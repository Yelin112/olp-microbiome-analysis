# =============================================================================
# TNF-alpha ELISA：不同处理组 DC 上清液 TNF-α 浓度组间比较
# 数据：projects/elisa/20260725/data/TNF_alpha_ELISA_two_batches_tidy.csv
#       （原始 batch_2 曾把 OD450 信号读数误标成第二批浓度数据，
#       已由用户修正为真实的第二批浓度读数，2 batch × 3 技术重复 = n=6/组）
# 说明：原始分组 G1~G9 已按用户最新说明重新定义 (G1、G4 与本次分析无关，已剔除)：
#         G2 = Uninfected (4MOSC1 无感染培养基)
#         G3 = LPS (LPS 阳性处理组)
#         G5 = Rm_HI (Rm 热灭活)
#         G6 = Rm (Rm 在细胞培养基中的上清)
#         G7/G8/G9 = iECM 感染 6h/12h/24h 上清
#       （与同批 IL-12 数据的分组定义一致）每组已有 n=6（2 batch × 3 技术
#       重复），无需再抽样，直接用全部 6 孔做组间比较，检验方式与 IL-12
#       最新分析一致：Welch ANOVA + Welch 两样本 t 检验。图上只标注
#       "各组 vs Uninfected"（6 条比较），统计表输出全部两两比较（21 条），
#       按用户要求两两比较不做多重检验校正。
# 使用 r_functions/lib/compare_plot 组间比较函数 + r-pub-toolkit 的 utils 主题/配色系统
# =============================================================================

library(tidyverse)
library(ggpubr)
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

DATA_DIR   <- file.path(REPO_ROOT, "projects/elisa/20260725/data")
OUTPUT_DIR <- file.path(REPO_ROOT, "projects/elisa/20260725/outputs")
FIG_DIR    <- file.path(OUTPUT_DIR, "figures")
STAT_DIR   <- file.path(OUTPUT_DIR, "stats")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 0. 分组重新定义（G1、G4 与本次分析无关，剔除；与 IL-12 脚本保持一致）
# -----------------------------------------------------------------------------
group_map <- c(
  G2 = "Uninfected",
  G3 = "LPS",
  G5 = "Rm_HI",
  G6 = "Rm",
  G7 = "iECM_6h",
  G8 = "iECM_12h",
  G9 = "iECM_24h"
)
group_order   <- unname(group_map)
control_group <- "Uninfected"

# -----------------------------------------------------------------------------
# 1. 读取数据
# -----------------------------------------------------------------------------
raw <- read.csv(
  file.path(DATA_DIR, "TNF_alpha_ELISA_two_batches_tidy.csv"),
  stringsAsFactors = FALSE
)

stopifnot(all(c("batch", "sample", "replicate", "concentration") %in% names(raw)))

well_data <- raw %>%
  filter(sample %in% names(group_map)) %>%
  mutate(group = factor(unname(group_map[sample]), levels = group_order))

message("各组样本量 (全部孔, n 应为 6 = 2 batch x 3 技术重复)：")
print(table(well_data$group))

# batch 内技术重复取均值 -> 每组每 batch 一个值，仅用于下方诊断图
batch_means <- well_data %>%
  group_by(batch, group) %>%
  summarise(mean_conc = mean(concentration, na.rm = TRUE), .groups = "drop")

# -----------------------------------------------------------------------------
# 2. Batch 一致性诊断图（不参与主统计，仅供肉眼核查两批趋势是否一致）
# -----------------------------------------------------------------------------
p_batch <- ggplot(batch_means, aes(x = group, y = mean_conc, color = batch, group = batch)) +
  geom_line(linewidth = 0.6, alpha = 0.7) +
  geom_point(size = 2.5) +
  scale_color_manual(values = get_colors("NPG", n = 2, type = "discrete")) +
  labs(
    title = "TNF-alpha concentration across batches",
    subtitle = "Diagnostic check: consistency of group trend across 2 batches",
    x = "Group", y = "TNF-alpha (pg/mL)", color = "Batch"
  ) +
  theme_pub_base() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(FIG_DIR, "10_tnf_batch_consistency.png"), p_batch, width = 18, height = 10, units = "cm", dpi = 300)
message("✓ 已生成图: 10_tnf_batch_consistency")

# -----------------------------------------------------------------------------
# 3. 组间比较主图（图上只标注 vs Uninfected）+ 统计检验（表格含全部两两比较）
# -----------------------------------------------------------------------------
key_pairs <- lapply(setdiff(group_order, control_group), function(g) c(control_group, g))
all_pairs <- combn(group_order, 2, simplify = FALSE)

# --- 整体检验 (Welch ANOVA) ---
aov_res <- oneway.test(concentration ~ group, data = well_data, var.equal = FALSE)
anova_tbl <- tibble(
  metric = "TNF-alpha", statistic = unname(aov_res$statistic),
  df = unname(aov_res$parameter[1]), p_value = aov_res$p.value
)

# --- 全部两两 Welch t 检验（21 对，按用户要求不做多重检验校正）---
ttest_tbl <- purrr::map_dfr(all_pairs, function(pr) {
  d1 <- well_data$concentration[well_data$group == pr[1]]
  d2 <- well_data$concentration[well_data$group == pr[2]]
  res <- suppressWarnings(t.test(d1, d2))
  tibble(
    group1 = pr[1], group2 = pr[2],
    n1 = length(d1), n2 = length(d2),
    statistic = unname(res$statistic),
    p_value = res$p.value,
    significance = case_when(
      p_value < 0.0001 ~ "****",
      p_value < 0.001  ~ "***",
      p_value < 0.01   ~ "**",
      p_value < 0.05   ~ "*",
      TRUE             ~ "ns"
    )
  )
})

# 图上标注用的子集：vs Uninfected
ttest_vs_control <- ttest_tbl %>% filter(group1 == control_group | group2 == control_group)

# --- 绘图（n=6/组；只标注 vs Uninfected 的比较，减少留白）---
aov_subtitle <- sprintf("Welch ANOVA, p = %.4f", anova_tbl$p_value)

p_main <- compare_plot(
  well_data,
  value.var = "concentration", group.by = "group",
  strategy = "auto",
  add_stat = "t.test", comparisons = key_pairs,
  stat_label = "p.signif", hide_ns = FALSE,
  step_increase = 0.10, y_expand = 0.2,
  palette = "NPG", theme_use = theme_pub_base,
  title = "TNF-alpha concentration", xlab = "Group", ylab = "TNF-alpha (pg/mL)"
) +
  labs(subtitle = aov_subtitle) +
  theme(
    legend.position = "none", plot.subtitle = element_text(size = 9, color = "grey30"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

g <- fix_panel_size(p_main, width = 14, height = 9.5)
dw <- grid::convertWidth(sum(g$widths), "cm", valueOnly = TRUE)
dh <- grid::convertHeight(sum(g$heights), "cm", valueOnly = TRUE)

ggsave(file.path(FIG_DIR, "11_tnf_group_comparison.pdf"), g, width = dw, height = dh, units = "cm")
ggsave(file.path(FIG_DIR, "11_tnf_group_comparison.png"), g, width = dw, height = dh, units = "cm", dpi = 300)
save_plot_pptx(
  p_main, file.path(FIG_DIR, "11_tnf_group_comparison.pptx"),
  panel_width = 14, panel_height = 9.5, title = "TNF-alpha concentration"
)
message("✓ 已生成图: 11_tnf_group_comparison")

# -----------------------------------------------------------------------------
# 4. 统计结果表
# -----------------------------------------------------------------------------
write.csv(anova_tbl, file.path(STAT_DIR, "tnf_welch_anova_overall.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(ttest_tbl, file.path(STAT_DIR, "tnf_pairwise_ttest_all.csv"), row.names = FALSE, fileEncoding = "UTF-8")
message("✓ 已生成统计结果表")

# -----------------------------------------------------------------------------
# 5. 生成 Markdown 分析报告（独立文件，与 IL-12 的 report.md 分开）
# -----------------------------------------------------------------------------
report_path <- file.path(OUTPUT_DIR, "report_tnf_alpha.md")

means_tbl <- well_data %>%
  group_by(group) %>%
  summarise(mean = mean(concentration), sem = sd(concentration) / sqrt(n()), n = n(), .groups = "drop")

lines <- c(
  "# TNF-α ELISA 分析报告：不同处理组 DC 上清液 TNF-α 浓度比较",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. 实验背景",
  "",
  "本批数据为同一批 DC 处理实验上清液中 TNF-α 浓度检测结果。原始分组为",
  "`G1`~`G9`，按最新说明重新定义并剔除了 `G1`、`G4`（与本次分析无关），",
  "剩余 7 组含义与 IL-12 数据一致：",
  "",
  "- **Uninfected**（原 G2）：4MOSC1 无感染培养基，作为基准对照",
  "- **LPS**（原 G3）：LPS 阳性处理组",
  "- **Rm_HI**（原 G5）：Rm 热灭活",
  "- **Rm**（原 G6）：Rm 在细胞培养基中的上清",
  "- **iECM_6h / iECM_12h / iECM_24h**（原 G7/G8/G9）：iECM 感染 6h/12h/24h 上清",
  "",
  "- 数据来自 2 个独立 batch，每 batch 内每组 3 个技术重复孔，",
  "  合计每组 n=6。",
  "- **数据质量说明**：原始 `TNF_alpha_ELISA_two_batches_tidy.csv` 中的",
  "  `batch_2` 曾被发现是 OD450 信号读数被误标成了第二批浓度数据",
  "  （数量级、来源截图标签都对不上 `batch_1` 的真实浓度），已由用户",
  "  修正为真实的第二批浓度读数后再用于本次分析。",
  "- 每组已有 n=6（2 batch × 3 技术重复），与 IL-12 最终分析采用的",
  "  样本量一致，因此直接使用全部 6 孔，未做随机抽样。",
  "- 技术重复并非完全独立的生物学重复，检验功效仍是以",
  "  『假定各孔独立』为前提换来的；下方诊断图显示 2 个 batch 组间趋势",
  "  是否一致，供解读时参考。",
  "",
  "各组样本量：",
  "",
  "| 组别 | n |",
  "|---|---|",
  paste0("| ", means_tbl$group, " | ", means_tbl$n, " |"),
  "",
  "## 2. 分析方法",
  "",
  "- 绘图使用项目库函数 `compare_plot()`（`r_functions/lib/compare_plot/compare_plot_optimized.R`），",
  "  每组 n=6 自动选用 `bar_points` 策略（均值柱 + SEM 误差线 + 全部散点）。",
  "- 配色使用 `palette_system.R` 的 `NPG` 期刊配色，主题使用 `theme_pub_base()`。",
  "- 图同时导出 PDF（矢量）、PNG（预览/嵌入文档）、PPTX（`export_pptx.R` 生成，可编辑）三种格式。",
  "- 统计检验（与 IL-12 最新分析方法一致）：",
  "  - 整体差异：Welch ANOVA（`oneway.test(var.equal = FALSE)`）。",
  "  - **图上只标注『各组 vs Uninfected』**（6 条比较），避免括号过多导致图面拥挤。",
  "  - **统计表（`stats/tnf_pairwise_ttest_all.csv`）包含全部两两比较**（7 组共 21 对，",
  "    Welch 两样本 t 检验），按用户要求**不做多重检验校正**（无 BH 校正列，",
  "    p 值均为原始 p 值）。",
  "",
  "## 3. 整体差异（Welch ANOVA）",
  "",
  "| 统计量 | p 值 |",
  "|---|---|",
  sprintf("| %.3f | %.4f |", anova_tbl$statistic, anova_tbl$p_value),
  "",
  "## 4. 各组均值 ± SEM",
  "",
  "| 组别 | 均值 (pg/mL) | SEM |",
  "|---|---|---|",
  sprintf("| %s | %.3f | %.3f |", means_tbl$group, means_tbl$mean, means_tbl$sem),
  "",
  "## 5. 各组 vs Uninfected（图上标注的比较）",
  "",
  "| 比较 | p 值 | 显著性 |",
  "|---|---|---|",
  sprintf(
    "| %s vs %s | %.4f | %s |",
    ttest_vs_control$group1, ttest_vs_control$group2, ttest_vs_control$p_value, ttest_vs_control$significance
  ),
  "",
  "完整的全部两两比较（21 对，未做多重校正）见 `stats/tnf_pairwise_ttest_all.csv`。",
  "",
  "## 6. 文件清单",
  "",
  "- `figures/10_tnf_batch_consistency.png` — 2 个 batch 组间趋势诊断图",
  "- `figures/11_tnf_group_comparison.pdf`（+ `.png` + `.pptx`） — 主比较图（图上标注 vs Uninfected）",
  "- `stats/tnf_welch_anova_overall.csv` — 整体 Welch ANOVA 检验结果",
  "- `stats/tnf_pairwise_ttest_all.csv` — 全部两两 Welch t 检验结果（21 对，未做多重检验校正）",
  ""
)

writeLines(lines, report_path, useBytes = TRUE)
message("✓ 已生成 TNF-α 分析报告: ", report_path)

message("\n全部完成。输出目录: ", OUTPUT_DIR)
