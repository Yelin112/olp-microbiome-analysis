library(tidyverse)
library(MASS)

set.seed(2026)

# 定义三个原始分布的中心和协方差矩阵
cluster1 <- mvrnorm(200, mu = c(2, 8), Sigma = matrix(c(0.5, 0.1, 0.1, 0.5), 2)) # 紧凑型
cluster2 <- mvrnorm(300, mu = c(8, 8), Sigma = matrix(c(2, -0.5, -0.5, 2), 2)) # 松散型
cluster3 <- mvrnorm(250, mu = c(5, 2), Sigma = matrix(c(3, 1.5, 1.5, 0.8), 2)) # 长条斜向

# 合并数据
df_raw <- rbind(cluster1, cluster2, cluster3) %>%
  as.data.frame() %>%
  rename(x = V1, y = V2) %>%
  mutate(original_group = factor(rep(c("A", "B", "C"), c(200, 300, 250))))

plot(df_raw)
# 定义你想要“引导”它们去的自定义引力中心
# 我们故意把引力中心设得比原始中心更远一点，观察拉伸效果
custom_centers <- data.frame(
  target_x = c(1, 9, 5),
  target_y = c(9, 9, 1)
)


# 增强函数
attract_logic <- function(df, centers_df, strength = 0.5, radius = 3) {
  df %>%
    rowwise() %>%
    mutate(
      # 寻找最近的自定义中心
      dists = list(sqrt(
        (centers_df$target_x - x)^2 + (centers_df$target_y - y)^2
      )),
      idx = which.min(dists),
      dist_min = dists[[idx]],

      # 方案 A: 固定比例收缩 (strength)
      # 方案 B: 距离衰减收缩 (pull)
      pull = exp(-(dist_min^2) / (2 * radius^2)) * strength,

      new_x = x + pull * (centers_df$target_x[idx] - x),
      new_y = y + pull * (centers_df$target_y[idx] - y)
    ) %>%
    ungroup()
}

# 应用算法：尝试不同的强度
df_enhanced <- attract_logic(df_raw, custom_centers, strength = 0.7, radius = 4)


library(patchwork)

# 图1：原始具有趋势的数据
p1 <- ggplot(df_enhanced, aes(x, y, color = original_group)) +
  geom_point(alpha = 0.4) +
  geom_point(
    data = custom_centers,
    aes(target_x, target_y),
    shape = 13,
    size = 5,
    color = "black"
  ) +
  theme_minimal() +
  labs(
    title = "原始数据 (具有初步聚类趋势)",
    subtitle = "十字为预定义的引力中心"
  ) +
  coord_fixed()

# 图2：增强后的数据
p2 <- ggplot(df_enhanced, aes(new_x, new_y, color = original_group)) +
  geom_point(alpha = 0.6) +
  geom_point(
    data = custom_centers,
    aes(target_x, target_y),
    shape = 13,
    size = 5,
    color = "black"
  ) +
  theme_minimal() +
  labs(title = "增强后的聚类", subtitle = "散点向中心发生非线性收缩") +
  coord_fixed()

p1 + p2 + plot_layout(guides = "collect")


library(tidyverse)
library(MASS)
library(patchwork)

#' 增强聚类视觉效果函数
#' @param data 原始数据框
#' @param x_col X轴列名
#' @param y_col Y轴列名
#' @param mode 算法模式: "linear" (固定比例), "gaussian" (距离衰减), "centroid" (向组均值收缩)
#' @param centers 自定义中心坐标框 (target_x, target_y)。若 mode="centroid" 则无需提供
#' @param group_col 组名列（仅在 mode="centroid" 时必须）
#' @param strength 强度系数 (0-1)
#' @param radius 高斯引力的半径（仅在 mode="gaussian" 时有效）

enrich_clusters <- function(
  data,
  x_col,
  y_col,
  mode = "linear",
  centers = NULL,
  group_col = NULL,
  strength = 0.5,
  radius = 2
) {
  # 内部计算：为每个点分配目标中心
  df_proc <- data %>% rename(tmp_x = !!sym(x_col), tmp_y = !!sym(y_col))

  # 模式 1: 向各组自身的均值中心收缩 (最常用)
  if (mode == "centroid") {
    group_centers <- df_proc %>%
      group_by(!!sym(group_col)) %>%
      summarise(
        target_x = mean(tmp_x),
        target_y = mean(tmp_y),
        .groups = "drop"
      )
    df_proc <- df_proc %>% left_join(group_centers, by = group_col)
  } else {
    # 模式 2 & 3: 向自定义的中心点收缩
    df_proc <- df_proc %>%
      rowwise() %>%
      mutate(
        dists = list(sqrt(
          (centers$target_x - tmp_x)^2 + (centers$target_y - tmp_y)^2
        )),
        idx = which.min(dists),
        target_x = centers$target_x[idx],
        target_y = centers$target_y[idx]
      ) %>%
      ungroup()
  }

  # 执行算法逻辑
  df_result <- df_proc %>%
    mutate(
      dist_to_target = sqrt((target_x - tmp_x)^2 + (target_y - tmp_y)^2),
      # 计算位移系数
      alpha = case_when(
        mode == "gaussian" ~ exp(-(dist_to_target^2) / (2 * radius^2)) *
          strength,
        TRUE ~ strength
      ),
      # 坐标变换
      new_x = tmp_x + alpha * (target_x - tmp_x),
      new_y = tmp_y + alpha * (target_y - tmp_y)
    ) %>%
    # 还原列名并清理多余列
    rename(!!sym(x_col) := tmp_x, !!sym(y_col) := tmp_y) %>%
    select(-target_x, -target_y, -dist_to_target, -alpha)

  return(df_result)
}

set.seed(888)
sim_data <- rbind(
  mvrnorm(500, mu = c(2, 2), Sigma = matrix(c(1.2, 0.8, 0.8, 1.2), 2)),
  mvrnorm(500, mu = c(5, 7), Sigma = matrix(c(2, -0.3, -0.3, 1), 2)),
  mvrnorm(500, mu = c(8, 3), Sigma = matrix(c(1, 0.5, 0.5, 1.5), 2))
) %>%
  as.data.frame() %>%
  rename(PC1 = V1, PC2 = V2) %>%
  mutate(group = rep(c("Control", "OLP_Rm", "OLP_Other"), each = 500))

# 定义三个预设引力中心
my_centers <- data.frame(target_x = c(1.5, 5, 8.5), target_y = c(1.5, 7.5, 2.5))


# 1. 线性模式 (Linear): 全局均匀收缩
res_lin <- enrich_clusters(
  sim_data,
  "PC1",
  "PC2",
  mode = "linear",
  centers = my_centers,
  strength = 0.4
)

# 2. 高斯模式 (Gaussian): 保护边缘，中心收缩 (调优参数 radius=3)
res_gau <- enrich_clusters(
  sim_data,
  "PC1",
  "PC2",
  mode = "gaussian",
  centers = my_centers,
  strength = 0.8,
  radius = 3
)

# 3. 质心模式 (Centroid): 自动向组内重心靠拢
res_cen <- enrich_clusters(
  sim_data,
  "PC1",
  "PC2",
  mode = "centroid",
  group_col = "group",
  strength = 0.5
)


library(tidyverse)
library(MASS)
library(uwot)
library(igraph)
library(patchwork)

# ================= 策略库 (Strategy Library) =================

# 1. 线性收缩
algo_linear <- function(p_mat, target_mat, strength, ...) {
  p_mat + strength * (target_mat - p_mat)
}

# 2. 高斯衰减引力
algo_gaussian <- function(p_mat, target_mat, strength, radius, ...) {
  # 计算每个点到其对应中心的距离
  dists <- sqrt(rowSums((target_mat - p_mat)^2))
  pull <- exp(-(dists^2) / (2 * radius^2)) * strength
  p_mat + pull * (target_mat - p_mat)
}

# 3. UMAP 流形变换
algo_umap <- function(p_mat, ...) {
  # UMAP 是全局算法，忽略 target_mat
  res <- uwot::umap(p_mat, n_components = 2, ...)
  return(as.matrix(res))
}

# 4. Fruchterman-Reingold 力导向布局
algo_force <- function(p_mat, strength = 1, ...) {
  dist_m <- as.matrix(dist(p_mat))
  # 建立近邻图：只连接距离最近的前 5% 的点
  adj_m <- dist_m < quantile(dist_m, 0.05)
  g <- graph_from_adjacency_matrix(adj_m, mode = "undirected", diag = FALSE)
  layout <- layout_with_fr(g, ...)
  # 将布局缩放到原始数据的大致尺度
  layout <- layout * (max(p_mat) / max(layout)) * strength
  return(as.matrix(layout))
}

# ================= 主调度函数 (Main Dispatcher) =================

enrich_clusters_pro <- function(
  data,
  x_col,
  y_col,
  method = "gaussian",
  group_col = NULL,
  centers = NULL,
  ...
) {
  # 1. 准备数据
  df_input <- data %>% rename(curr_x = !!sym(x_col), curr_y = !!sym(y_col))
  p_mat <- as.matrix(df_input[, c("curr_x", "curr_y")])

  # 2. 中心点识别逻辑 (针对线性、高斯、质心模式)
  target_mat <- NULL
  if (method %in% c("linear", "gaussian", "centroid")) {
    if (method == "centroid" || (is.null(centers) && !is.null(group_col))) {
      # 自动计算各组均值作为中心
      grp_centers <- df_input %>%
        group_by(!!sym(group_col)) %>%
        summarise(t_x = mean(curr_x), t_y = mean(curr_y), .groups = "drop")
      temp_df <- df_input %>% left_join(grp_centers, by = group_col)
      target_mat <- as.matrix(temp_df[, c("t_x", "t_y")])
      # 质心模式本质上是特定 centers 下的线性收缩
      if (method == "centroid") method <- "linear"
    } else if (!is.null(centers)) {
      # 寻找最近的自定义中心
      dists_to_centers <- as.matrix(dist(rbind(p_mat, as.matrix(centers))))
      # 提取点到中心的子矩阵
      point_to_center_dists <- dists_to_centers[
        1:nrow(p_mat),
        (nrow(p_mat) + 1):ncol(dists_to_centers)
      ]
      nearest_idx <- apply(point_to_center_dists, 1, which.min)
      target_mat <- as.matrix(centers[nearest_idx, ])
    }
  }

  # 3. 动态调用算法
  algo_func <- match.fun(paste0("algo_", method))

  # 执行计算
  if (method %in% c("umap", "force")) {
    new_coords <- algo_func(p_mat, ...)
  } else {
    new_coords <- algo_func(p_mat, target_mat, ...)
  }

  # 4. 封装结果
  df_output <- df_input %>%
    mutate(new_x = new_coords[, 1], new_y = new_coords[, 2]) %>%
    rename(!!sym(x_col) := curr_x, !!sym(y_col) := curr_y)

  return(df_output)
}


set.seed(123)
# 模拟三个正态分布，增加方差使其重叠
data_sim <- rbind(
  mvrnorm(300, mu = c(2, 2), Sigma = matrix(c(2, 1, 1, 2), 2)),
  mvrnorm(300, mu = c(5, 5), Sigma = matrix(c(2, -1, -1, 2), 2)),
  mvrnorm(300, mu = c(8, 2), Sigma = matrix(c(2, 0.5, 0.5, 2), 2))
) %>%
  as.data.frame() %>%
  rename(dim1 = V1, dim2 = V2) %>%
  mutate(group = rep(c("Control", "OLP_Rm", "OLP_Other"), each = 300))


# 1. 质心收缩 (最为稳健)
res_cen <- enrich_clusters_pro(
  data_sim,
  "dim1",
  "dim2",
  method = "centroid",
  group_col = "group",
  strength = 0.5
)

# 2. 高斯增强 (保护离群点)
res_gau <- enrich_clusters_pro(
  data_sim,
  "dim1",
  "dim2",
  method = "gaussian",
  group_col = "group",
  strength = 0.8,
  radius = 2
)

# 3. UMAP 变换 (重塑拓扑)
res_umap <- enrich_clusters_pro(
  data_sim,
  "dim1",
  "dim2",
  method = "umap",
  n_neighbors = 20,
  min_dist = 0.1
)

# 4. 力导向布局 (增强网络感)
res_force <- enrich_clusters_pro(
  data_sim,
  "dim1",
  "dim2",
  method = "force",
  niter = 500
)


# 绘图函数
p_theme <- function(df, title) {
  ggplot(df, aes(new_x, new_y, color = group)) +
    geom_point(alpha = 0.5, size = 1) +
    theme_minimal() +
    labs(title = title) +
    theme(legend.position = "none")
}

p0 <- ggplot(data_sim, aes(dim1, dim2, color = group)) +
  geom_point(alpha = 0.5, size = 1) +
  theme_minimal() +
  labs(title = "Original")

(p0 + p_theme(res_cen, "Centroid (Linear)")) /
  (p_theme(res_gau, "Gaussian (Decay)") + p_theme(res_umap, "UMAP Manifold")) /
  (p_theme(res_force, "Force-Directed"))


library(tidyverse)
library(MASS)
library(uwot)
library(igraph)
library(patchwork)

# ================= 策略库 (Strategy Library) =================

# 1. 线性收缩
algo_linear <- function(p_mat, target_mat, strength, ...) {
  p_mat + strength * (target_mat - p_mat)
}

# 2. 高斯衰减引力
algo_gaussian <- function(p_mat, target_mat, strength, radius, ...) {
  dists <- sqrt(rowSums((target_mat - p_mat)^2))
  pull <- exp(-(dists^2) / (2 * radius^2)) * strength
  p_mat + pull * (target_mat - p_mat)
}

# 3. UMAP 流形变换
algo_umap <- function(p_mat, ...) {
  res <- uwot::umap(p_mat, n_components = 2, ...)
  return(as.matrix(res))
}

# 4. 力导向布局
algo_force <- function(p_mat, strength = 1, ...) {
  dist_m <- as.matrix(dist(p_mat))
  adj_m <- dist_m < quantile(dist_m, 0.05)
  g <- graph_from_adjacency_matrix(adj_m, mode = "undirected", diag = FALSE)
  layout <- layout_with_fr(g, ...)
  layout <- layout * (max(p_mat) / max(layout)) * strength
  return(as.matrix(layout))
}

# ================= 主调度函数 (Main Dispatcher) =================

enrich_clusters_pro <- function(
  data,
  x_col,
  y_col,
  method = "gaussian",
  group_col = NULL,
  centers = NULL,
  ...
) {
  # 1. 数据准备
  df_input <- data %>% rename(curr_x = !!sym(x_col), curr_y = !!sym(y_col))
  p_mat <- as.matrix(df_input[, c("curr_x", "curr_y")])

  # 2. 中心点识别逻辑 (支持按组匹配)
  target_mat <- NULL

  if (method %in% c("linear", "gaussian", "centroid")) {
    # A. 自动计算重心模式 (Centroid Mode)
    if (method == "centroid" || (is.null(centers) && !is.null(group_col))) {
      grp_centers <- df_input %>%
        group_by(!!sym(group_col)) %>%
        summarise(t_x = mean(curr_x), t_y = mean(curr_y), .groups = "drop")
      temp_df <- df_input %>% left_join(grp_centers, by = group_col)
      target_mat <- as.matrix(temp_df[, c("t_x", "t_y")])
      if (method == "centroid") method <- "linear"
    } else if (!is.null(centers)) {
      # B. 按组精准匹配模式 (Group-specific Centers)
      # 如果 centers 包含分组列，则进行 Join 匹配
      if (!is.null(group_col) && group_col %in% colnames(centers)) {
        # 提取坐标并重命名以防冲突
        centers_mapped <- centers %>%
          dplyr::select(all_of(group_col), t_x = 1, t_y = 2)

        temp_df <- df_input %>%
          left_join(centers_mapped, by = group_col)

        target_mat <- as.matrix(temp_df[, c("t_x", "t_y")])
      } else {
        # C. 最近邻匹配模式 (Nearest Neighbor)
        # 仅取前两列作为坐标
        c_mat <- as.matrix(centers[, 1:2])
        dists_to_centers <- as.matrix(dist(rbind(p_mat, c_mat)))
        point_to_center_dists <- dists_to_centers[
          1:nrow(p_mat),
          (nrow(p_mat) + 1):ncol(dists_to_centers)
        ]
        nearest_idx <- apply(as.matrix(point_to_center_dists), 1, which.min)
        target_mat <- c_mat[nearest_idx, ]
      }
    }
  }

  # 3. 动态调度策略
  algo_func <- match.fun(paste0("algo_", method))

  if (method %in% c("umap", "force")) {
    new_coords <- algo_func(p_mat, ...)
  } else {
    new_coords <- algo_func(p_mat, target_mat, ...)
  }

  # 4. 结果还原
  df_output <- df_input %>%
    mutate(new_x = new_coords[, 1], new_y = new_coords[, 2]) %>%
    rename(!!sym(x_col) := curr_x, !!sym(y_col) := curr_y)

  return(df_output)
}


# --- 准备数据 ---
set.seed(123)
data_sim <- rbind(
  mvrnorm(300, mu = c(2, 2), Sigma = matrix(c(2, 1, 1, 2), 2)),
  mvrnorm(300, mu = c(5, 5), Sigma = matrix(c(2, -1, -1, 2), 2)),
  mvrnorm(300, mu = c(8, 2), Sigma = matrix(c(2, 0.5, 0.5, 2), 2))
) %>%
  as.data.frame() %>%
  rename(dim1 = V1, dim2 = V2) %>%
  mutate(group = rep(c("Control", "OLP_Rm", "OLP_Other"), each = 300))

# --- 定义按组指定的质心 (Group-specific Centers) ---
# 我们故意把 Control 往左拉，OLP_Rm 往上拉，OLP_Other 往右拉
my_custom_centers <- data.frame(
  x = c(0, 5, 10),
  y = c(0, 10, 0),
  group = c("Control", "OLP_Rm", "OLP_Other") # 对应分组
)

# --- 运行验证 ---
# 1. 原始质心模式 (自动计算重心)
res_auto <- enrich_clusters_pro(
  data_sim,
  "dim1",
  "dim2",
  method = "centroid",
  group_col = "group",
  strength = 0.5
)

# 2. 自定义分组质心模式 (精准匹配我们设定的 0, 5, 10 坐标)
res_custom <- enrich_clusters_pro(
  data_sim,
  "dim1",
  "dim2",
  method = "linear",
  group_col = "group",
  centers = my_custom_centers,
  strength = 0.5
)

# --- 可视化对比 ---
p1 <- ggplot(data_sim, aes(dim1, dim2, color = group)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Original Overlap Data")

p2 <- ggplot(res_auto, aes(new_x, new_y, color = group)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Method: Centroid (Auto Mean)")

p3 <- ggplot(res_custom, aes(new_x, new_y, color = group)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Method: Linear (Custom Group Centers)")

(p1 | p2 | p3) + plot_layout(guides = "collect")
