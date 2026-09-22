# enrich_clusters_pro
# 对已知分组的 2D 散点图做视觉增强：将各群体的点向群体中心收缩，
# 使重叠群体在视觉上更清晰分离，不改变数据的统计性质。
#
# 依赖：tidyverse, uwot（method="umap"时）, igraph（method="force"时）

# ── 内部：策略函数 ────────────────────────────────────────────────────────────

.algo_linear <- function(p_mat, target_mat, strength, ...) {
  p_mat + strength * (target_mat - p_mat)
}

.algo_gaussian <- function(p_mat, target_mat, strength, radius, ...) {
  dists <- sqrt(rowSums((target_mat - p_mat)^2))
  pull  <- exp(-(dists^2) / (2 * radius^2)) * strength
  p_mat + pull * (target_mat - p_mat)
}

.algo_umap <- function(p_mat, ...) {
  as.matrix(uwot::umap(p_mat, n_components = 2, ...))
}

.algo_force <- function(p_mat, strength = 1, ...) {
  dist_m <- as.matrix(dist(p_mat))
  adj_m  <- dist_m < quantile(dist_m, 0.05)
  g      <- igraph::graph_from_adjacency_matrix(adj_m, mode = "undirected", diag = FALSE)
  layout <- igraph::layout_with_fr(g, ...)
  layout * (max(p_mat) / max(layout)) * strength
}

# GMM 自动中心：拟合高斯混合模型，提取成分均值作为 softpull 中心
.algo_gmm <- function(p_mat, strength, radius, k = NULL, ...) {
  if (!requireNamespace("mclust", quietly = TRUE))
    stop("[enrich_clusters_pro] gmm 模式需要安装 mclust：install.packages('mclust')")
  suppressPackageStartupMessages(library(mclust))   # 必须 attach 才能解析内部函数
  fit   <- Mclust(p_mat, G = if (is.null(k)) 1:9 else k, verbose = FALSE)
  means <- fit$parameters$mean
  c_mat <- if (is.matrix(means)) t(means) else matrix(means, 1, 2)
  message(sprintf("[enrich_clusters_pro] GMM 检测到 %d 个成分", nrow(c_mat)))
  .algo_softpull(p_mat, c_mat, strength, radius)
}

# KDE 密度峰：在 2D 核密度估计网格上找局部极大值，作为 softpull 中心
.algo_kde_modal <- function(p_mat, strength, radius,
                             grid_n = 100, min_rel_height = 0.15, ...) {
  kde  <- MASS::kde2d(p_mat[, 1], p_mat[, 2], h = radius * 2,
                      n = grid_n,
                      lims = c(range(p_mat[, 1]), range(p_mat[, 2])))
  z    <- kde$z
  thr  <- max(z) * min_rel_height
  modes <- NULL
  for (i in 2:(grid_n - 1)) {
    for (j in 2:(grid_n - 1)) {
      patch <- z[(i - 1):(i + 1), (j - 1):(j + 1)]
      if (z[i, j] >= thr && z[i, j] == max(patch))
        modes <- rbind(modes, c(kde$x[i], kde$y[j]))
    }
  }
  if (is.null(modes)) {
    warning("[enrich_clusters_pro] kde_modal 未检测到局部峰值，返回原始数据")
    return(p_mat)
  }
  message(sprintf("[enrich_clusters_pro] KDE 检测到 %d 个密度峰", nrow(modes)))
  .algo_softpull(p_mat, modes, strength, radius)
}

# Mean Shift（近似版）：向邻域加权均值迭代收敛
# 用 n_ref 个采样点作参照，避免 O(n²) 全量计算
.algo_meanshift <- function(p_mat, strength, radius,
                             n_ref = 2000, n_iter = 1, ...) {
  n   <- nrow(p_mat)
  ref <- if (n > n_ref) p_mat[sample(n, n_ref), ] else p_mat
  b2  <- rowSums(ref^2)
  cur <- p_mat
  for (iter in seq_len(n_iter)) {
    a2  <- rowSums(cur^2)
    d2  <- outer(a2, b2, "+") - 2 * tcrossprod(cur, ref)  # n × n_ref
    w   <- exp(-d2 / (2 * radius^2))
    ws  <- rowSums(w)
    cur <- cbind(drop(w %*% ref[, 1]) / ws,
                 drop(w %*% ref[, 2]) / ws)
  }
  p_mat + strength * (cur - p_mat)
}

# 软吸引：每个点同时被所有中心吸引，引力按 Gaussian 距离加权平均
# 不做硬分配 → 无 Voronoi 边界，适合无标签连续分布
.algo_softpull <- function(p_mat, c_mat, strength, radius, ...) {
  n_pts  <- nrow(p_mat)
  n_ctr  <- nrow(c_mat)
  # 计算每个点到每个中心的平方距离（向量化，n_pts × n_ctr）
  d2_mat <- matrix(0, n_pts, n_ctr)
  for (k in seq_len(n_ctr)) {
    d2_mat[, k] <- rowSums((p_mat - matrix(c_mat[k, ], n_pts, 2, byrow = TRUE))^2)
  }
  # Gaussian 权重，按行归一化（softmax-like）
  w_mat      <- exp(-d2_mat / (2 * radius^2))
  w_mat      <- w_mat / rowSums(w_mat)
  # 加权平均中心作为每个点的目标
  target_mat <- w_mat %*% c_mat
  p_mat + strength * (target_mat - p_mat)
}

# ── 内部：根据群体重心间距自动估算 radius ─────────────────────────────────────
.estimate_radius <- function(p_mat, group_vec) {
  unique_grps <- unique(group_vec)
  if (length(unique_grps) < 2) {
    return(diff(range(p_mat)) / 8)
  }
  centroids <- t(vapply(unique_grps, function(g) {
    colMeans(p_mat[group_vec == g, , drop = FALSE])
  }, numeric(2)))
  mean(dist(centroids)) / 4
}

# ── 内部：检测坐标系类型 ──────────────────────────────────────────────────────
# logicle 变换后：值域约 0~4.5，无极端大值
.detect_coord_type <- function(p_mat) {
  if (max(p_mat) < 6 && min(p_mat) > -1.5) "logicle" else "linear"
}

# ── 主函数 ───────────────────────────────────────────────────────────────────

#' 散点图聚类视觉增强（流式数据优化版）
#'
#' @param data       data.frame，含坐标列和分组列
#' @param x_col      X 轴列名（字符串）
#' @param y_col      Y 轴列名（字符串）
#' @param method     算法模式：
#'   "gaussian"   距离衰减引力，推荐用于流式（默认）
#'   "softpull"   软吸引：每点被所有中心按 Gaussian 加权平均吸引，无 Voronoi 边界
#'   "gmm"        GMM 自动中心：mclust 拟合混合高斯，成分均值作为 softpull 中心
#'   "kde_modal"  KDE 密度峰：核密度局部极大值作为 softpull 中心
#'   "meanshift"  近似 Mean Shift：向采样邻域加权均值迭代收敛
#'   "centroid"   向各组均值均匀收缩
#'   "linear"     向指定中心均匀收缩
#'   "umap"       UMAP 流形变换（⚠️ 破坏坐标轴生物学含义，不建议流式使用）
#'   "force"      力导向布局（⚠️ 同上）
#' @param group_col  分组列名；centroid 模式必须，其他模式有 group_col 时可按组处理
#' @param centers    自定义中心坐标（见下）：
#'   NULL（默认）：centroid 自动计算，其他走最近邻
#'   含 group_col 列的 data.frame：按组精准匹配
#'   只含坐标的 data.frame：最近邻匹配
#' @param coord_type 坐标系类型，影响 strength/radius 默认值：
#'   "logicle"（流式 logicle 变换后，0~4.5 范围）→ strength=0.35, radius=自动
#'   "linear"（线性/PCA 等）→ strength=0.5, radius=自动
#'   "auto"（默认）→ 根据数据范围自动判断
#' @param strength   收缩强度 0~1；NULL 时根据 coord_type 自动设置
#' @param radius     高斯引力半径；NULL 时根据群体中心间距自动估算
#' @param min_events 事件数低于此值的群体跳过增强（保护稀有群体）；默认 50
#' @param ...        传递给 umap/force 的额外参数（如 n_neighbors, niter）
#' @return 原始 data.frame 加上 new_x, new_y 两列（变换后坐标）
enrich_clusters_pro <- function(
  data,
  x_col,
  y_col,
  method     = "gaussian",
  group_col  = NULL,
  centers    = NULL,
  coord_type = "auto",
  strength   = NULL,
  radius     = NULL,
  min_events = 50,
  ...
) {
  # ── 0. 参数预检 ─────────────────────────────────────────────────────────────
  if (method %in% c("umap", "force")) {
    message(sprintf(
      "[enrich_clusters_pro] method='%s' 会破坏坐标轴的生物学含义，",
      method
    ), "不建议用于流式数据。如需继续，请忽略此提示。")
  }

  # ── 1. 准备坐标矩阵 ─────────────────────────────────────────────────────────
  df_input <- data |>
    dplyr::rename(curr_x = !!rlang::sym(x_col), curr_y = !!rlang::sym(y_col))
  p_mat <- as.matrix(df_input[, c("curr_x", "curr_y")])

  # ── 2. 自动检测坐标系，设置 strength/radius 默认值 ─────────────────────────
  if (coord_type == "auto") {
    coord_type <- .detect_coord_type(p_mat)
    message(sprintf("[enrich_clusters_pro] coord_type 自动检测为: %s", coord_type))
  }

  if (is.null(strength)) {
    strength <- if (coord_type == "logicle") 0.35 else 0.5
  }

  # ── 3. 中心点解析 ──────────────────────────────────────────────────────────
  target_mat  <- NULL
  c_mat       <- NULL   # softpull 专用：中心坐标矩阵
  group_vec   <- if (!is.null(group_col)) df_input[[group_col]] else NULL

  # softpull：直接提取中心矩阵，不做逐点最近邻分配
  if (method == "softpull") {
    if (is.null(centers)) stop("[enrich_clusters_pro] softpull 模式需要提供 centers 参数")
    coord_cols <- if (!is.null(group_col) && group_col %in% colnames(centers)) {
      setdiff(colnames(centers), group_col)[1:2]
    } else {
      colnames(centers)[1:2]
    }
    c_mat <- as.matrix(centers[, coord_cols])
  }

  if (method %in% c("linear", "gaussian", "centroid")) {

    if (method == "centroid" || (is.null(centers) && !is.null(group_col))) {
      # A. 自动重心模式
      grp_centers <- df_input |>
        dplyr::group_by(!!rlang::sym(group_col)) |>
        dplyr::summarise(t_x = mean(curr_x), t_y = mean(curr_y), .groups = "drop")
      temp_df    <- df_input |> dplyr::left_join(grp_centers, by = group_col)
      target_mat <- as.matrix(temp_df[, c("t_x", "t_y")])
      if (method == "centroid") method <- "linear"

    } else if (!is.null(centers)) {

      if (!is.null(group_col) && group_col %in% colnames(centers)) {
        # B. 按组精准匹配
        coord_cols     <- setdiff(colnames(centers), group_col)[1:2]
        centers_mapped <- centers |>
          dplyr::select(dplyr::all_of(c(group_col, coord_cols))) |>
          dplyr::rename(t_x = dplyr::all_of(coord_cols[1]),
                        t_y = dplyr::all_of(coord_cols[2]))
        temp_df    <- df_input |> dplyr::left_join(centers_mapped, by = group_col)
        target_mat <- as.matrix(temp_df[, c("t_x", "t_y")])

      } else {
        # C. 最近邻匹配
        c_mat               <- as.matrix(centers[, 1:2])
        all_dists           <- as.matrix(dist(rbind(p_mat, c_mat)))
        p2c                 <- all_dists[1:nrow(p_mat),
                                         (nrow(p_mat) + 1):ncol(all_dists)]
        nearest_idx         <- apply(as.matrix(p2c), 1, which.min)
        target_mat          <- c_mat[nearest_idx, ]
      }
    }
  }

  # ── 4. min_events：稀有群体保护 ────────────────────────────────────────────
  if (!is.null(group_vec) && !is.null(target_mat) && min_events > 0) {
    grp_counts   <- table(group_vec)
    small_groups <- names(grp_counts)[grp_counts < min_events]

    if (length(small_groups) > 0) {
      skip_mask <- group_vec %in% small_groups
      # 稀有群体的目标 = 自身位置（位移为 0，不移动）
      target_mat[skip_mask, ] <- p_mat[skip_mask, ]
      message(sprintf(
        "[enrich_clusters_pro] %d 个稀有群体（共 %d 个事件）事件数 < %d，已跳过增强: %s",
        length(small_groups), sum(skip_mask), min_events,
        paste(small_groups, collapse = ", ")
      ))
    }
  }

  # ── 5. radius 自动估算 ────────────────────────────────────────────────────
  if (method %in% c("gaussian", "softpull", "meanshift", "kde_modal", "gmm") &&
      is.null(radius)) {
    radius <- if (!is.null(group_vec)) {
      .estimate_radius(p_mat, group_vec)
    } else if (coord_type == "logicle") {
      0.5
    } else {
      diff(range(p_mat)) / 8
    }
    message(sprintf("[enrich_clusters_pro] radius 自动估算: %.3f", radius))
  }

  # ── 6. 调度算法 ────────────────────────────────────────────────────────────
  algo_func  <- match.fun(paste0(".algo_", method))
  new_coords <- if (method %in% c("umap", "force")) {
    algo_func(p_mat, strength = strength, ...)
  } else if (method == "softpull") {
    algo_func(p_mat, c_mat, strength = strength, radius = radius, ...)
  } else if (method %in% c("gmm", "kde_modal", "meanshift")) {
    algo_func(p_mat, strength = strength, radius = radius, ...)
  } else {
    algo_func(p_mat, target_mat, strength = strength, radius = radius, ...)
  }

  # ── 7. 结果封装 ────────────────────────────────────────────────────────────
  df_input |>
    dplyr::mutate(new_x = new_coords[, 1], new_y = new_coords[, 2]) |>
    dplyr::rename(!!rlang::sym(x_col) := curr_x,
                  !!rlang::sym(y_col) := curr_y)
}
