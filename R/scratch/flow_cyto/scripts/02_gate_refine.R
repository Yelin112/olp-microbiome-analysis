# 02_gate_refine.R
# 操作流程：
#   ① 逐段运行（Ctrl+Enter），不要整体 source()
#   ② 运行到"看图"标记处，观察散点图，记下目标群体所在的坐标范围
#   ③ 在"定 gate"段填写坐标（图上看到的值直接填，无需换算）
#   ④ 运行后续段：筛选事件 → 验证着色 → 导出 FCS

library(flowCore)
library(tidyverse)
# flowCore::filter 被 dplyr::filter 覆盖，全脚本统一用 flowCore::filter() 显式调用

source("R/scratch/flow_cyto/config/panel_T_cell.R")

DATA_DIR <- "R/scratch/flow_cyto/data/S1"
OUTPUT_DIR <- "R/scratch/flow_cyto/output/gated_fcs"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 有效荧光通道（排除仪器空通道）────────────────────────────────────────────
active_channels <- PANEL$channels |>
  keep(\(ch) !isTRUE(ch$empty)) |>
  map_chr("channel_area")

# ════════════════════════════════════════════════════════════════════════════
# STEP 1  读取 FCS + 荧光通道变换
#         ★ 修改下面这一行切换变换方式，gate 坐标和散点图坐标无需换算 ★
#
#   "logicle" → 值域 0~4.5，处理负值，FlowJo Biexponential 对应
#   "log10"   → 值域 0~5.4，轴标 10²/10³/10⁴，FlowJo Log 对应
# ════════════════════════════════════════════════════════════════════════════
DISP_TRANSFORM <- "log10" # ← 在此切换："log10" 或 "logicle"

TARGET_FILE <- "export_OXA-MN_Lymphocytes.fcs" # 修改这里以切换输入文件

fcs_raw <- read.FCS(
  file.path(DATA_DIR, TARGET_FILE),
  truncate_max_range = FALSE,
  emptyValue = FALSE
)

if (DISP_TRANSFORM == "logicle") {
  # Logicle（Biexponential）变换：w=0.5, t=262144, m=4.5, a=0（FlowJo 默认参数）
  lgcl_tl <- transformList(
    active_channels,
    logicleTransform(w = 0.5, t = 262144, m = 4.5, a = 0)
  )
  fcs_trans <- transform(fcs_raw, lgcl_tl)
  df_trans <- as.data.frame(exprs(fcs_trans))
} else {
  # Log10 变换：值≤0 夹到 1（log10(1)=0），避免 -Inf
  log10_safe <- function(x) log10(pmax(x, 1))
  fcs_trans <- fcs_raw
  exprs_mat <- exprs(fcs_trans)
  exprs_mat[, active_channels] <- log10_safe(exprs_mat[, active_channels])
  exprs(fcs_trans) <- exprs_mat
  df_trans <- as.data.frame(exprs(fcs_trans))
}

cat(sprintf(
  "已加载: %s  (%d 个事件)  变换: %s\n",
  TARGET_FILE,
  nrow(df_trans),
  DISP_TRANSFORM
))

# ════════════════════════════════════════════════════════════════════════════
# STEP 2  ★ 看图 ★
#         运行下面这段，观察散点图中各群体的位置
#         横轴 = CD4 (FL4-A)，纵轴 = CD8 (FL5-A)
#         logicle 变换后坐标通常在 0~4.5 之间
# ════════════════════════════════════════════════════════════════════════════
x_ch <- "FL1-A" # CD3
y_ch <- "FL8-A" # CD45

p_base <- function(df, x, y, title) {
  ggplot(df, aes(.data[[x]], .data[[y]])) +
    geom_hex(bins = 80) +
    scale_fill_gradientn(
      colors = c("#f0f0f0", "#6baed6", "#08306b"),
      name = "N"
    ) +
    labs(title = title, x = x, y = y) +
    theme_minimal(base_size = 11) +
    theme(aspect.ratio = 1, legend.position = "none")
}

p_raw <- ggplot(df_trans, aes(.data[[x_ch]], .data[[y_ch]])) +
  geom_hex(bins = 80) +
  scale_fill_gradientn(
    colors = c("#f0f0f0", "#6baed6", "#08306b"),
    name = "N"
  ) +
  labs(
    title = paste0(TARGET_FILE, "  |  CD4 vs CD8"),
    x = "CD4 (PC5.5)  [logicle]",
    y = "CD8 (PC7)  [logicle]",
    caption = "记下 CD4+ 群和 CD8+ 群的大致坐标，填入 STEP 3"
  ) +
  theme_minimal(base_size = 13) +
  theme(aspect.ratio = 1)

print(p_raw)

# ════════════════════════════════════════════════════════════════════════════
# STEP 2b  四种无标签算法对比
#          运行后输出 2×2 对比图，选满意的方法后保留在 METHOD_USE 里
# ════════════════════════════════════════════════════════════════════════════
library(patchwork)
library(MASS) # kde_modal 需要
source(file.path(Sys.getenv("R_TOOLKIT_ROOT"), "r_functions/dev/enrich_clusters/function.R"))

# softpull 的手动中心（仅 softpull 需要；其他方法自动估计）
enrich_centers <- data.frame(
  x = c(2.5, 3.0, 4.0, 4.5),
  y = c(2.5, 2.8, 3.8, 5.0)
)

# ── 公共参数 ──────────────────────────────────────────────────────────────
STRENGTH <- 0.5
# radius 留 NULL → 自动估算；或手动指定，如 0.8

p_hex <- function(df, title) {
  ggplot(df, aes(new_x, new_y)) +
    geom_hex(bins = 80) +
    scale_fill_gradientn(
      colors = c("#f0f0f0", "#6baed6", "#08306b"),
      name = "N"
    ) +
    labs(title = title, x = x_ch, y = y_ch) +
    theme_minimal(base_size = 11) +
    theme(aspect.ratio = 1, legend.position = "none")
}

res_softpull <- enrich_clusters_pro(
  df_trans,
  x_ch,
  y_ch,
  method = "softpull",
  centers = enrich_centers,
  coord_type = DISP_TRANSFORM,
  strength = STRENGTH
)

res_gmm <- enrich_clusters_pro(
  df_trans,
  x_ch,
  y_ch,
  method = "gmm",
  coord_type = DISP_TRANSFORM,
  strength = STRENGTH
)

res_kde <- enrich_clusters_pro(
  df_trans,
  x_ch,
  y_ch,
  method = "kde_modal",
  coord_type = DISP_TRANSFORM,
  strength = STRENGTH
)

res_ms <- enrich_clusters_pro(
  df_trans,
  x_ch,
  y_ch,
  method = "meanshift",
  coord_type = DISP_TRANSFORM,
  strength = STRENGTH
)

p_compare <- (p_base(df_trans, x_ch, y_ch, "原始") |
  p_hex(res_softpull, "softpull（手动中心）")) /
  (p_hex(res_gmm, "gmm（自动混合高斯）") |
    p_hex(res_kde, "kde_modal（密度峰）")) /
  (p_hex(res_ms, "meanshift（邻域均值）") |
    plot_spacer())
print(p_compare)
# ★ 确定方法后，修改下面 METHOD_USE，再运行 STEP 2c 导出

METHOD_USE <- "softpull" # ← 改为 gmm / kde_modal / meanshift

df_enriched <- switch(
  METHOD_USE,
  softpull = res_softpull,
  gmm = res_gmm,
  kde_modal = res_kde,
  meanshift = res_ms
)

# ════════════════════════════════════════════════════════════════════════════
# STEP 2c  反变换导出：变换坐标 → 原始荧光值 → FCS
#          自动根据 DISP_TRANSFORM 选择对应的逆变换
#          ⚠️ 文件名含 _visual_only，不得用于定量分析
# ════════════════════════════════════════════════════════════════════════════

orig_x <- exprs(fcs_raw)[, x_ch]
orig_y <- exprs(fcs_raw)[, y_ch]

if (DISP_TRANSFORM == "log10") {
  # log10 逆变换：10^x，精确无误差
  raw_x_new <- as.integer(pmin(pmax(round(10^df_enriched$new_x), 0), 262144))
  raw_y_new <- as.integer(pmin(pmax(round(10^df_enriched$new_y), 0), 262144))
  # log10 对原始值 ≤ 1 的事件无效（全被压缩到 0）
  # 这些点（补偿后负值、近零值）保留原始坐标，不做增强
  raw_x_new <- ifelse(orig_x > 1, raw_x_new, orig_x)
  raw_y_new <- ifelse(orig_y > 1, raw_y_new, orig_y)
} else {
  # logicle 逆变换：正向建查找表 → approxfun 线性插值（60k 点，误差 < 4 FU）
  lookup_raw <- seq(-5000, 262144, length.out = 60000)
  dummy_ff <- flowFrame(matrix(
    lookup_raw,
    ncol = 1,
    dimnames = list(NULL, x_ch)
  ))
  dummy_lgcl <- transformList(
    x_ch,
    logicleTransform(w = 0.5, t = 262144, m = 4.5, a = 0)
  )
  lookup_lgcl <- exprs(transform(dummy_ff, dummy_lgcl))[, x_ch]
  inv_lgcl <- approxfun(lookup_lgcl, lookup_raw, rule = 2)
  raw_x_new <- as.integer(pmin(
    pmax(round(inv_lgcl(df_enriched$new_x)), -5000),
    262144
  ))
  raw_y_new <- as.integer(pmin(
    pmax(round(inv_lgcl(df_enriched$new_y)), -5000),
    262144
  ))
}

fcs_visual <- fcs_raw
exprs(fcs_visual)[, x_ch] <- raw_x_new
exprs(fcs_visual)[, y_ch] <- raw_y_new

out_visual <- file.path(
  OUTPUT_DIR,
  paste0(tools::file_path_sans_ext(TARGET_FILE), "_visual_only.fcs")
)
write.FCS(fcs_visual, out_visual)
cat(sprintf(
  "\n已导出视觉增强 FCS（仅供展示）:\n  %s\n  替换通道: %s, %s\n",
  out_visual,
  x_ch,
  y_ch
))

# ════════════════════════════════════════════════════════════════════════════
# STEP 3  ★ 定 gate ★
#         根据上图观察，填写矩形框的四个边界（log10 变换后的坐标）
#         log10 坐标参考：log10(100)=2, log10(1000)=3, log10(10000)=4
#         CD4+ 群：CD4 高、CD8 低
#         CD8+ 群：CD4 低、CD8 高
# ════════════════════════════════════════════════════════════════════════════

# --- 在这里修改数值 ---
gate_cd4 <- rectangleGate(
  filterId = "CD4+",
  "FL4-A" = c(2.0, Inf), # CD4 下限（看图估计）
  "FL5-A" = c(-Inf, 2.0) # CD8 上限（看图估计）
)

gate_cd8 <- rectangleGate(
  filterId = "CD8+",
  "FL4-A" = c(-Inf, 2.0), # CD4 上限
  "FL5-A" = c(2.0, Inf) # CD8 下限
)
# ----------------------

# 在图上叠加 gate 框，确认位置是否准确
p_gated <- p_raw +
  annotate(
    "rect",
    xmin = gate_cd4@min[["FL4-A"]],
    xmax = Inf,
    ymin = -Inf,
    ymax = gate_cd4@max[["FL5-A"]],
    color = "#e41a1c",
    fill = NA,
    linewidth = 0.9,
    alpha = 0.8
  ) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = gate_cd8@max[["FL4-A"]],
    ymin = gate_cd8@min[["FL5-A"]],
    ymax = Inf,
    color = "#377eb8",
    fill = NA,
    linewidth = 0.9,
    alpha = 0.8
  ) +
  annotate(
    "label",
    x = 3.5,
    y = 0.5,
    label = "CD4+",
    color = "#e41a1c",
    size = 4
  ) +
  annotate(
    "label",
    x = 0.5,
    y = 3.5,
    label = "CD8+",
    color = "#377eb8",
    size = 4
  ) +
  labs(title = paste0(TARGET_FILE, "  |  Gate 预览（确认后再继续）"))

print(p_gated)
# ★ 如果框的位置不对，回到 STEP 3 调整数值，重新运行

# ════════════════════════════════════════════════════════════════════════════
# STEP 4  应用 gate，筛选事件
#         技巧：gate 作用于变换后的 fcs_trans，取出事件索引
#         用索引从原始 fcs_raw 里取对应行 → 导出的 FCS 保持原始荧光值
# ════════════════════════════════════════════════════════════════════════════
idx_cd4 <- flowCore::filter(fcs_trans, gate_cd4)@subSet
idx_cd8 <- flowCore::filter(fcs_trans, gate_cd8)@subSet

fcs_cd4_out <- fcs_raw[idx_cd4, ] # 原始值，CD4+ 事件
fcs_cd8_out <- fcs_raw[idx_cd8, ] # 原始值，CD8+ 事件

n_total <- nrow(exprs(fcs_raw))
cat(sprintf(
  "\n%s 总事件: %d\nCD4+ gate: %d 事件 (%.1f%%)\nCD8+ gate: %d 事件 (%.1f%%)\n",
  TARGET_FILE,
  n_total,
  sum(idx_cd4),
  100 * sum(idx_cd4) / n_total,
  sum(idx_cd8),
  100 * sum(idx_cd8) / n_total
))

# ════════════════════════════════════════════════════════════════════════════
# STEP 5  ★ 验证 ★  着色查看 gate 结果是否合理
# ════════════════════════════════════════════════════════════════════════════
df_verify <- df_trans |>
  mutate(
    pop = case_when(
      idx_cd4 ~ "CD4+",
      idx_cd8 ~ "CD8+",
      TRUE ~ "其他"
    )
  )

p_verify <- ggplot(df_verify, aes(.data[[x_ch]], .data[[y_ch]], color = pop)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(
    values = c("CD4+" = "#e41a1c", "CD8+" = "#377eb8", "其他" = "#d9d9d9"),
    name = "群体"
  ) +
  labs(
    title = "Gate 验证：CD4+ / CD8+ 着色",
    x = "CD4 (PC5.5)  [logicle]",
    y = "CD8 (PC7)  [logicle]"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_minimal(base_size = 13) +
  theme(aspect.ratio = 1)

print(p_verify)
# ★ 如果着色不对，回 STEP 3 调整 gate 坐标

# ════════════════════════════════════════════════════════════════════════════
# STEP 6  导出 FCS（FlowJo 可直接重新导入）
#         导出的是原始荧光值（未做 logicle 变换）
# ════════════════════════════════════════════════════════════════════════════
stem <- tools::file_path_sans_ext(TARGET_FILE)
out_cd4 <- file.path(OUTPUT_DIR, paste0(stem, "_CD4pos.fcs"))
out_cd8 <- file.path(OUTPUT_DIR, paste0(stem, "_CD8pos.fcs"))

write.FCS(fcs_cd4_out, out_cd4)
write.FCS(fcs_cd8_out, out_cd8)

cat(sprintf(
  "\n已导出（原始荧光值，FlowJo 兼容）:\n  %s\n  %s\n",
  out_cd4,
  out_cd8
))
