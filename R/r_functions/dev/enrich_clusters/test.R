# test.R — enrich_clusters_pro 测试用例

library(tidyverse)
library(MASS)
library(patchwork)

source("R/r_functions/dev/enrich_clusters/function.R")

# ── 测试数据 A：通用重叠正态分布（linear 坐标系）────────────────────────────
set.seed(123)
data_linear <- rbind(
  mvrnorm(300, mu = c(2, 2), Sigma = matrix(c(2, 1, 1, 2), 2)),
  mvrnorm(300, mu = c(5, 5), Sigma = matrix(c(2, -1, -1, 2), 2)),
  mvrnorm(300, mu = c(8, 2), Sigma = matrix(c(2, .5, .5, 2), 2))
) |>
  as.data.frame() |>
  rename(dim1 = V1, dim2 = V2) |>
  mutate(group = rep(c("Control", "OLP_Rm", "OLP_Other"), each = 300))

# ── 测试数据 B：模拟 logicle 变换后的流式数据（0~4.5 范围）──────────────────
set.seed(42)
data_flow <- rbind(
  mvrnorm(800, mu = c(0.5, 0.5), Sigma = matrix(c(.15, .02, .02, .15), 2)), # 死细胞/碎片
  mvrnorm(2000, mu = c(2.0, 1.0), Sigma = matrix(c(.20, .05, .05, .15), 2)), # CD4+
  mvrnorm(1500, mu = c(1.0, 2.5), Sigma = matrix(c(.15, .03, .03, .20), 2)), # CD8+
  mvrnorm(30, mu = c(3.0, 3.0), Sigma = matrix(c(.10, .01, .01, .10), 2)) # 稀有群体
) |>
  as.data.frame() |>
  rename(FL4_A = V1, FL5_A = V2) |>
  mutate(pop = rep(c("debris", "CD4+", "CD8+", "rare"), c(800, 2000, 1500, 30)))

p_base <- function(df, x, y, title, col = "group") {
  ggplot(df, aes(.data[[x]], .data[[y]], color = .data[[col]])) +
    geom_point(alpha = 0.4, size = 0.8) +
    theme_minimal(base_size = 11) +
    labs(title = title) +
    theme(legend.position = "bottom", aspect.ratio = 1)
}

# ════════════════════════════════════════════════════════════════════════════
# T1  coord_type 自动检测
# ════════════════════════════════════════════════════════════════════════════
cat("── T1: coord_type 自动检测 ──\n")
# linear 数据应检测为 "linear"
res_auto_lin <- enrich_clusters_pro(
  data_linear,
  "dim1",
  "dim2",
  method = "centroid",
  group_col = "group"
)
# flow 数据应检测为 "logicle"
res_auto_flow <- enrich_clusters_pro(
  data_flow,
  "FL4_A",
  "FL5_A",
  method = "centroid",
  group_col = "pop"
)
stopifnot(all(c("new_x", "new_y") %in% colnames(res_auto_flow)))
cat("T1 通过\n\n")

# ════════════════════════════════════════════════════════════════════════════
# T2  流式场景：coord_type="logicle"，radius 自动估算
# ════════════════════════════════════════════════════════════════════════════
cat("── T2: logicle 模式，radius 自动估算 ──\n")
res_flow_gau <- enrich_clusters_pro(
  data_flow,
  "FL4_A",
  "FL5_A",
  method = "gaussian",
  group_col = "pop",
  coord_type = "logicle"
)
cat("T2 通过\n\n")

# ════════════════════════════════════════════════════════════════════════════
# T3  min_events 稀有群体保护
#     rare 群体 30 个事件，min_events=50 时应跳过
# ════════════════════════════════════════════════════════════════════════════
cat("── T3: min_events 稀有群体保护 ──\n")
res_min <- enrich_clusters_pro(
  data_flow,
  "FL4_A",
  "FL5_A",
  method = "gaussian",
  group_col = "pop",
  coord_type = "logicle",
  min_events = 50
)
# rare 群体的 new_x/new_y 应与原始坐标相同（未移动）
rare_orig <- data_flow |> filter(pop == "rare")
rare_new <- res_min |> filter(pop == "rare")
stopifnot(all.equal(rare_orig$FL4_A, rare_new$new_x))
stopifnot(all.equal(rare_orig$FL5_A, rare_new$new_y))
cat("T3 通过：rare 群体（30 事件）坐标未发生移动\n\n")

# ════════════════════════════════════════════════════════════════════════════
# T4  umap/force 不建议用于流式：应打印 message 但不报错
# ════════════════════════════════════════════════════════════════════════════
# library(uwot)
# cat("── T4: umap 警告 ──\n")
# res_umap <- enrich_clusters_pro(data_linear, "dim1", "dim2",
#                                  method = "umap", n_neighbors = 15)
# cat("T4 通过\n\n")

# ════════════════════════════════════════════════════════════════════════════
# T5  手动指定参数（覆盖自动值）
# ════════════════════════════════════════════════════════════════════════════
cat("── T5: 手动 strength=0.6, radius=0.4 ──\n")
res_manual <- enrich_clusters_pro(
  data_flow,
  "FL4_A",
  "FL5_A",
  method = "gaussian",
  group_col = "pop",
  coord_type = "logicle",
  strength = 0.7,
  radius = 0.5
)
cat("T5 通过\n\n")

# ── 可视化对比 ────────────────────────────────────────────────────────────────
p0 <- p_base(data_flow, "FL4_A", "FL5_A", "原始流式数据", col = "pop")
p1 <- p_base(
  res_flow_gau,
  "new_x",
  "new_y",
  "gaussian (logicle, auto radius)",
  col = "pop"
)
p2 <- p_base(
  res_min,
  "new_x",
  "new_y",
  "min_events=50 (rare 不移动)",
  col = "pop"
)
p3 <- p_base(
  res_manual,
  "new_x",
  "new_y",
  "手动 strength=0.6 radius=0.4",
  col = "pop"
)

(p0 | p1) / (p2 | p3) + plot_layout(guides = "collect")

cat("所有测试通过\n")
