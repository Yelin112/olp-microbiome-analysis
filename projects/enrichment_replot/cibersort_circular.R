source("E:/打工人/Zeng/2023-3 OLP Microbiome/分析/环状热图_claude_2.R")

library(dplyr)
library(tidyr)
library(RColorBrewer)
library(circlize)

outdir <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/富集重画"

# ── 1. 读取数据 ───────────────────────────────────────────────────────────────
df <- read.csv(
  file.path(outdir, "Cibersort_data_long.csv"),
  row.names = 1,
  stringsAsFactors = FALSE
)
df$Composition <- as.numeric(df$Composition)

# ── 2. 样本排序：OLP 在前，HC 在后 ───────────────────────────────────────────
sample_anno <- df %>%
  distinct(sample, group) %>%
  mutate(group = factor(group, levels = c("OLP", "HC"))) %>%
  arrange(group, sample)

ordered_samples <- as.character(sample_anno$sample)
olp_samples <- as.character(sample_anno$sample[sample_anno$group == "OLP"])
hc_samples <- as.character(sample_anno$sample[sample_anno$group == "HC"])

# OLP | 空白分隔列 | HC
value_cols_sep <- c(olp_samples, "SEPARATOR", hc_samples)

# ── 3. 宽格式：行 = 细胞类型，列 = 样本 ──────────────────────────────────────
df_wide <- df %>%
  select(Celltype, sample, Composition) %>%
  pivot_wider(names_from = sample, values_from = Composition) %>%
  as.data.frame()

# ── 4. 细胞类型分组 ───────────────────────────────────────────────────────────
celltype_group_map <- c(
  "B.cells.naive" = "B cells",
  "B.cells.memory" = "B cells",
  "Plasma.cells" = "B cells",
  "T.cells.CD8" = "T cells",
  "T.cells.CD4.naive" = "T cells",
  "T.cells.CD4.memory.resting" = "T cells",
  "T.cells.CD4.memory.activated" = "T cells",
  "T.cells.follicular.helper" = "T cells",
  "T.cells.regulatory..Tregs." = "T cells",
  "T.cells.gamma.delta" = "T cells",
  "NK.cells.resting" = "NK cells",
  "NK.cells.activated" = "NK cells",
  "Monocytes" = "Myeloid",
  "Macrophages.M0" = "Myeloid",
  "Macrophages.M1" = "Myeloid",
  "Macrophages.M2" = "Myeloid",
  "Dendritic.cells.resting" = "Myeloid",
  "Dendritic.cells.activated" = "Myeloid",
  "Mast.cells.resting" = "Other",
  "Mast.cells.activated" = "Other",
  "Eosinophils" = "Other",
  "Neutrophils" = "Other"
)

df_wide$CellGroup <- celltype_group_map[df_wide$Celltype]

# ── 5. Wilcoxon 检验（OLP vs HC，每种细胞类型）───────────────────────────────
stat_res <- lapply(unique(df$Celltype), function(ct) {
  sub <- df[df$Celltype == ct, ]
  wt <- wilcox.test(Composition ~ group, data = sub, exact = FALSE)
  data.frame(
    Celltype = ct,
    p_value = wt$p.value,
    mean_OLP = mean(sub$Composition[sub$group == "OLP"]),
    mean_HC = mean(sub$Composition[sub$group == "HC"])
  )
}) %>%
  bind_rows()

# 有向 -log10(p)：正 = OLP > HC，负 = HC > OLP
stat_res$neg_log10_p <- -log10(stat_res$p_value + 1e-10)
stat_res$direction <- sign(stat_res$mean_OLP - stat_res$mean_HC)
stat_res$neg_log10_p_dir <- stat_res$neg_log10_p * stat_res$direction

# log2FC (OLP/HC)，小伪计数避免除零
stat_res$log2FC <- log2((stat_res$mean_OLP + 1e-4) / (stat_res$mean_HC + 1e-4))

# 稀疏度：该细胞类型在所有样本中 Composition = 0 的比例
sparsity_res <- df %>%
  group_by(Celltype) %>%
  summarise(sparsity = mean(Composition == 0), .groups = "drop")
stat_res <- left_join(stat_res, sparsity_res, by = "Celltype")

# 全样本最大 Composition（峰值丰度）
max_res <- df %>%
  group_by(Celltype) %>%
  summarise(max_comp = max(Composition), .groups = "drop")
stat_res <- left_join(stat_res, max_res, by = "Celltype")

# ── 6. 合并数据，按细胞大类和名称排序 ────────────────────────────────────────
df_plot <- df_wide %>%
  left_join(
    stat_res[, c(
      "Celltype",
      "neg_log10_p_dir",
      "neg_log10_p",
      "p_value",
      "log2FC",
      "sparsity",
      "max_comp"
    )],
    by = "Celltype"
  ) %>%
  mutate(
    CellGroup = factor(
      celltype_group_map[Celltype],
      levels = c("B cells", "T cells", "NK cells", "Myeloid", "Other")
    )
  ) %>%
  arrange(CellGroup, Celltype)

# 插入 OLP/HC 分隔列（NA_real_ 确保通过 is.numeric 校验）
df_plot$SEPARATOR <- NA_real_

df_plot <- df_plot %>%
  select(
    CellGroup,
    Celltype,
    all_of(value_cols_sep), # OLP列 + SEPARATOR + HC列
    neg_log10_p_dir,
    neg_log10_p,
    p_value,
    log2FC,
    sparsity,
    max_comp
  )

cat(
  "Data ready:",
  nrow(df_plot),
  "cell types,",
  length(ordered_samples),
  "samples\n"
)
cat(
  "OLP samples:",
  sum(sample_anno$group == "OLP"),
  "| HC samples:",
  sum(sample_anno$group == "HC"),
  "\n"
)

# ── 7. 颜色映射 ───────────────────────────────────────────────────────────────
# 热图：RdYlBu 翻转后跳过前3个深蓝色，从浅蓝(#ABD9E9)开始 → 黄 → 深红
col_fun <- colorRamp2(
  seq(0, 1, length.out = 8),
  rev(brewer.pal(11, "RdYlBu"))[4:11]
)

# 细胞大类颜色
group_pal <- c(
  "B cells" = "#4DBBD5",
  "T cells" = "#E64B35",
  "NK cells" = "#00A087",
  "Myeloid" = "#F39B7F",
  "Other" = "#8491B4"
)

# 显著性轨道颜色函数：红=OLP高，蓝=HC高，灰=ns
sig_thr <- -log10(0.05) # ≈ 1.301

dir_color_fn <- function(x) {
  ifelse(
    x >= sig_thr,
    "#E64B35", # OLP > HC，显著
    ifelse(
      x <= -sig_thr,
      "#4DBBD5", # HC > OLP，显著
      "#BBBBBB"
    )
  ) # 不显著
}

# log2FC 颜色（正=OLP高=红，负=HC高=蓝）
fc_color_fn <- function(x) {
  ifelse(x > 0, "#E64B35", "#4DBBD5")
}

# 稀疏度颜色（灰度：稀疏度越高越浅）
sparsity_color_fn <- function(x) {
  grays <- colorRampPalette(c("#555555", "#EEEEEE"))(100)
  grays[pmin(pmax(round(x * 99) + 1, 1), 100)]
}

# ── 8. 轨道配置 ───────────────────────────────────────────────────────────────
track_configs <- list(
  # 轨道1：有向显著性（-log10p * direction）
  list(
    column = "neg_log10_p_dir",
    type = "bars",
    height = 0.07,
    label = "-log10(p)",
    color = dir_color_fn,
    bar_width = 0.8,
    border = NA
  ),
  # 轨道2：log2 fold change（OLP/HC），展示效应量
  list(
    column = "log2FC",
    type = "bars",
    height = 0.07,
    label = "log2FC",
    color = fc_color_fn,
    bar_width = 0.8,
    border = NA
  ),
  # 轨道3：稀疏度（Composition=0 的样本比例），展示检出率
  list(
    column = "sparsity",
    type = "bars",
    height = 0.06,
    label = "Sparsity",
    color = sparsity_color_fn,
    bar_width = 0.8,
    border = NA
  )
)

# ── 9. 绘图函数（PDF 和 PNG 复用） ───────────────────────────────────────────
draw_plot <- function() {
  plot_advanced_circos(
    data = df_plot,
    group_col = "CellGroup",
    feature_col = "Celltype",
    value_cols = value_cols_sep, # OLP | SEPARATOR | HC
    group_palette = group_pal,
    heatmap_col = col_fun,
    track_configs = track_configs,
    show_rownames = TRUE,
    show_colnames = FALSE,
    gap_degree = 90,
    gap_between = 2,
    start_degree = 90,
    value_range = c(0, 1),
    plot_params = list(
      heatmap_rownames_cex = 0.68,
      heatmap_track_height = 0.3, # 增加热图宽度
      group_label_cex = 0.9,
      group_bg_alpha = 0.6
    ),
    group_label_params = list(
      track_height = 0.05
    ),
    show_legend = TRUE,
    legend_params = list(
      title = "Composition\n(scaled)",
      x = 0.82,
      y = 0.92
    )
  )

  # 显著性轨道图例
  legend(
    "bottomright",
    legend = c("OLP > HC (sig.)", "HC > OLP (sig.)", "Not significant"),
    fill = c("#E64B35", "#4DBBD5", "#BBBBBB"),
    border = NA,
    bty = "n",
    cex = 0.8,
    title = "Track1: -log10(p)"
  )

  # 样本分组说明
  legend(
    "bottomleft",
    legend = c(
      paste0("OLP n=", sum(sample_anno$group == "OLP"), " (inner)"),
      paste0("HC  n=", sum(sample_anno$group == "HC"), " (outer)")
    ),
    fill = c("#E64B35", "#4DBBD5"),
    border = NA,
    bty = "n",
    cex = 0.8,
    title = "Heatmap columns"
  )

  # log2FC 和 sparsity 图例
  legend(
    "topright",
    legend = c("OLP > HC", "HC > OLP"),
    fill = c("#E64B35", "#4DBBD5"),
    border = NA,
    bty = "n",
    cex = 0.8,
    title = "Track2: log2FC"
  )

  legend(
    "topleft",
    legend = c("0%  (always detected)", "100% (never detected)"),
    fill = c("#555555", "#EEEEEE"),
    border = "grey60",
    bty = "n",
    cex = 0.8,
    title = "Track3: Sparsity"
  )

  title(
    "CIBERSORT Immune Cell Composition: OLP vs HC",
    cex.main = 1.2,
    line = -1.5
  )
}

# ── 10. 保存 ─────────────────────────────────────────────────────────────────
pdf(
  file.path(outdir, "cibersort_circular_heatmap.pdf"),
  width = 12,
  height = 12
)
draw_plot()
dev.off()

png(
  file.path(outdir, "cibersort_circular_heatmap.png"),
  width = 12,
  height = 12,
  units = "in",
  res = 300
)
draw_plot()
dev.off()

cat("\nDone. Output:\n")
cat("  cibersort_circular_heatmap.pdf\n")
cat("  cibersort_circular_heatmap.png\n")

# ── 11. 显著性结果汇总 ────────────────────────────────────────────────────────
cat("\nSignificant cell types (raw p < 0.05):\n")
sig_cells <- stat_res[
  stat_res$p_value < 0.05,
  c("Celltype", "p_value", "mean_OLP", "mean_HC", "direction")
]
sig_cells$direction_label <- ifelse(
  sig_cells$direction > 0,
  "OLP > HC",
  "HC > OLP"
)
print(sig_cells[order(sig_cells$p_value), ])
