# =============================================================================
# CFU 杀菌实验：可视化 + 组间比较
# 数据：collaborators/Gilly/CFU_killing_assay/data/cfu_clean_data.csv
#       （已剔除 nanozyme_1mg 与 exp_1，见 01_clean_data.R 顶部说明）
# 设计：control vs nanozyme_0.5mg，两组，3 次独立实验（exp_2~exp_4）同批次
#       配对测定（同一次实验、同一起始菌液），故按配对方式检验。
#       每组 n=3，compare_plot 的 strategy="auto" 会自动选用 pure_scatter
#       （全部散点 + 均值 ± SE），不依赖其内置的非配对检验，改为图外单独计算
#       配对 t 检验（t.test(paired=TRUE)）并以 subtitle 形式标注（遵循
#       PLOTTING_CONVENTIONS.md "整体检验放 subtitle、图内星号仅做探索性展示"
#       的既有约定）。
# 指标：killing_percent, survival_percent, log10_cfu_mid, cfu_mid。
#       前两者是物理上限为 100% 的比例指标，套用 % 轴硬上限逻辑（先建图查看
#       实际 y 轴范围，超过 100 才收紧重画）；后两者为浓度类指标，不设上限。
# =============================================================================

library(tidyverse)

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

BASE_DIR   <- file.path(PROJ, "collaborators/Gilly/CFU_killing_assay")
DATA_DIR   <- file.path(BASE_DIR, "data")
OUTPUT_DIR <- file.path(BASE_DIR, "output")
FIG_DIR    <- file.path(OUTPUT_DIR, "figures")
STAT_DIR   <- file.path(OUTPUT_DIR, "stats")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. 读取干净数据 + 英文展示标签（图内文字统一英文，脚本注释用中文）
# -----------------------------------------------------------------------------
clean <- read.csv(file.path(DATA_DIR, "cfu_clean_data.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    experiment = factor(experiment),
    group = factor(group, levels = c("control", "nanozyme_0.5mg")),
    group_label = factor(
      recode(as.character(group), control = "Control", `nanozyme_0.5mg` = "Nanozyme 0.5 mg/mL"),
      levels = c("Control", "Nanozyme 0.5 mg/mL")
    )
  )

message("各组样本量（n=3/组，3 次实验配对）：")
print(table(clean$group))

# -----------------------------------------------------------------------------
# 2. 配对 t 检验（按 experiment 配对，4 个指标一起做 BH 校正）
# -----------------------------------------------------------------------------
paired_test <- function(metric) {
  wide <- clean %>%
    select(experiment, group, all_of(metric)) %>%
    pivot_wider(names_from = group, values_from = all_of(metric))
  ctrl  <- wide$control
  treat <- wide$`nanozyme_0.5mg`
  res <- t.test(ctrl, treat, paired = TRUE)
  tibble(
    metric = metric,
    n_pairs = length(ctrl),
    mean_control = mean(ctrl),
    mean_treatment = mean(treat),
    statistic = unname(res$statistic),
    p_raw = res$p.value
  )
}

metrics <- c("killing_percent", "survival_percent", "log10_cfu_mid", "cfu_mid")
ttest_tbl <- purrr::map_dfr(metrics, paired_test) %>%
  mutate(p_BH = p.adjust(p_raw, method = "BH"))

write.csv(ttest_tbl, file.path(STAT_DIR, "paired_ttest_control_vs_0.5mg.csv"), row.names = FALSE)
message("✓ 已保存配对 t 检验结果")
print(ttest_tbl)

subtitle_for <- function(metric) {
  row <- ttest_tbl %>% filter(metric == !!metric)
  sprintf("Paired t-test, n=%d pairs, p = %.3f (BH-adj p = %.3f)", row$n_pairs, row$p_raw, row$p_BH)
}

# -----------------------------------------------------------------------------
# 3. 绘图函数：pure_scatter（n<5 自动选用）+ subtitle 标注配对检验结果
#    百分比指标套用硬上限逻辑：先建图看实际 y 轴范围，超过 100 才收紧重画。
# -----------------------------------------------------------------------------
build_plot <- function(value_var, y_label, cap = Inf) {
  make_plot <- function(y_exp) {
    compare_plot(
      clean,
      value.var = value_var, group.by = "group_label",
      strategy = "auto", add_stat = "none",
      y_expand = y_exp,
      palette = "NPG", theme_use = theme_pub_base,
      title = "CFU Killing Assay", xlab = NULL, ylab = y_label
    ) +
      labs(subtitle = subtitle_for(value_var))
  }
  p <- make_plot(0.15)
  top <- ggplot_build(p)$layout$panel_params[[1]]$y.range[2]
  if (is.finite(cap) && top > cap) {
    p <- make_plot(0.02) +
      scale_y_continuous(limits = c(NA, cap), expand = expansion(mult = c(0.05, 0)))
  }
  p
}

p_killing  <- build_plot("killing_percent", "Killing Rate (%)", cap = 100)
p_survival <- build_plot("survival_percent", "Survival Rate (%)", cap = 100)
p_log10cfu <- build_plot("log10_cfu_mid", "log10(CFU/mL)", cap = Inf)
p_cfu <- build_plot("cfu_mid", "CFU/mL", cap = Inf) +
  scale_y_log10(labels = scales::label_scientific())

# -----------------------------------------------------------------------------
# 4. 导出：每张图 PDF + PNG(300dpi) + PPTX，另打包一份 all_plots.pptx
# -----------------------------------------------------------------------------
plots <- list(
  list(p = p_killing,  stub = "01_killing_percent",  w = 9, h = 9),
  list(p = p_survival, stub = "02_survival_percent", w = 9, h = 9),
  list(p = p_log10cfu, stub = "03_log10_cfu_mid",    w = 9, h = 9),
  list(p = p_cfu,      stub = "04_cfu_mid",          w = 9, h = 9)
)

for (plt in plots) {
  g <- fix_panel_size(plt$p, width = plt$w, height = plt$h)
  dw <- grid::convertWidth(sum(g$widths), "cm", valueOnly = TRUE)
  dh <- grid::convertHeight(sum(g$heights), "cm", valueOnly = TRUE)
  ggsave(file.path(FIG_DIR, paste0(plt$stub, ".pdf")), g, width = dw, height = dh, units = "cm")
  ggsave(file.path(FIG_DIR, paste0(plt$stub, ".png")), g, width = dw, height = dh, units = "cm", dpi = 300)
  save_plot_pptx(plt$p, file.path(FIG_DIR, paste0(plt$stub, ".pptx")), panel_width = plt$w, panel_height = plt$h)
}

named_plots <- setNames(
  list(p_killing, p_survival, p_log10cfu, p_cfu),
  c("Killing Rate", "Survival Rate", "log10 CFU/mL", "CFU/mL")
)
save_plots_pptx(named_plots, file.path(FIG_DIR, "all_plots.pptx"), panel_width = 9, panel_height = 9)

message("✓ 图表已生成 (", FIG_DIR, ")")
message("\n全部完成。输出目录: ", OUTPUT_DIR)
