# =============================================================================
# 共用设置：路径、库源码、浓度映射、配色、通用绘图/统计函数
# 由 01_prepare_data.R / 02_analyze_batch.R / 03_compare_batches.R 共同 source。
# 不单独运行。
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

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

BASE_DIR <- file.path(PROJ, "collaborators/Gilly/Growth_curve/SM")
RAW_DIR <- file.path(BASE_DIR, "raw")
DATA_DIR <- file.path(BASE_DIR, "data")
OUTPUT_DIR <- file.path(BASE_DIR, "output")
FIG_DIR <- file.path(OUTPUT_DIR, "figures")
STAT_DIR <- file.path(OUTPUT_DIR, "stats")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 浓度映射（A 为对照，B~F 按浓度从低到高排列；所有批次共用同一套材料浓度体系）
# 若浓度需要更正，改这里即可，下游脚本无需改动。
# -----------------------------------------------------------------------------
conc_map <- tibble::tribble(
  ~group , ~conc_ugml , ~label        ,
  "A"    ,   0        , "0 (Control)" ,
  "F"    ,  12.5      , "12.5"        ,
  "E"    ,  25        , "25"          ,
  "D"    ,  50        , "50"          ,
  "C"    , 100        , "100"         ,
  "B"    , 200        , "200"
)
write.csv(
  conc_map,
  file.path(DATA_DIR, "concentration_map.csv"),
  row.names = FALSE
)

control_label <- conc_map$label[conc_map$group == "A"]
dose_pairs <- purrr::map(setdiff(conc_map$label, control_label), function(l) {
  c(l, control_label)
})
n_groups <- length(conc_map$label)

group_colors <- get_colors("NPG", n = n_groups)
names(group_colors) <- conc_map$label

trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

save_fig <- function(p, dir, stub, w, h) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  g <- fix_panel_size(p, width = w, height = h)
  dw <- grid::convertWidth(sum(g$widths), "cm", valueOnly = TRUE)
  dh <- grid::convertHeight(sum(g$heights), "cm", valueOnly = TRUE)
  ggsave(
    file.path(dir, paste0(stub, ".pdf")),
    g,
    width = dw,
    height = dh,
    units = "cm"
  )
  ggsave(
    file.path(dir, paste0(stub, ".png")),
    g,
    width = dw,
    height = dh,
    units = "cm",
    dpi = 300
  )
  save_plot_pptx(
    p,
    file.path(dir, paste0(stub, ".pptx")),
    panel_width = w,
    panel_height = h
  )
}

build_stat_plot <- function(df, value_var, y_label, title) {
  compare_plot(
    df,
    value.var = value_var,
    group.by = "group",
    strategy = "auto",
    add_stat = "t.test",
    comparisons = dose_pairs,
    stat_label = "p.signif",
    hide_ns = FALSE,
    step_increase = 0.14,
    y_expand = 0.5,
    palette = "NPG",
    theme_use = theme_pub_base,
    title = title,
    xlab = "Nanozyme Concentration (ug/mL)",
    ylab = y_label
  ) +
    labs(subtitle = "Welch t-test vs Control, unadjusted p (n=3/group)") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

plot_growth_curve <- function(summary_data, mean_col, sd_col, y_label, title) {
  df <- summary_data %>%
    rename(y_mean = all_of(mean_col), y_sd = all_of(sd_col))
  ggplot(df, aes(x = time_h, y = y_mean, color = group, fill = group)) +
    geom_ribbon(
      aes(ymin = y_mean - y_sd, ymax = y_mean + y_sd),
      alpha = 0.15,
      color = NA
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.6) +
    scale_color_manual(values = group_colors, name = "Nanozyme (ug/mL)") +
    scale_fill_manual(values = group_colors, name = "Nanozyme (ug/mL)") +
    labs(
      title = title,
      subtitle = "Mean \u00b1 SD (n=3 wells/group)",
      x = "Time (h)",
      y = y_label
    ) +
    theme_pub_base()
}

run_ttest_table <- function(replicate_metrics, value_var) {
  purrr::map_dfr(dose_pairs, function(pr) {
    d1 <- replicate_metrics[[value_var]][replicate_metrics$group == pr[1]]
    d2 <- replicate_metrics[[value_var]][replicate_metrics$group == pr[2]]
    res <- suppressWarnings(t.test(d1, d2))
    tibble(
      metric = value_var,
      group1 = pr[1],
      group2 = pr[2],
      n1 = length(d1),
      n2 = length(d2),
      mean1 = mean(d1),
      mean2 = mean(d2),
      statistic = unname(res$statistic),
      df = unname(res$parameter),
      p_value = res$p.value
    )
  }) %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      significance = case_when(
        p_value < 0.0001 ~ "****",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}

# -----------------------------------------------------------------------------
# 解析要处理的 BATCH_ID：优先用脚本里手动设的 BATCH_ID 变量；
# 若通过 `Rscript xxx.R <batch_id>` 传了命令行参数，参数会覆盖它。
# -----------------------------------------------------------------------------
resolve_batch_id <- function(default_batch_id) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && nzchar(args[1])) args[1] else default_batch_id
}
