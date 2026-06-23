#' ============================================================================
#' filter_abundance: 微生物物种丰度过滤函数
#' ============================================================================
#'
#' 基于流行率、丰度占比、最小计数、方差等指标过滤低丰度/不重要的物种。
#' 支持多条件组合过滤，默认参数基于微生物组学领域的常用阈值。
#'
#' 输入: 行=物种(features)，列=样本(samples) 的矩阵或 data.frame。
#'
#' @param abund_table matrix 或 data.frame，行=物种, 列=样本。
#'   可以是原始计数 (counts) 或相对丰度 (relative abundance)。
#' @param method 过滤策略，默认 "combined"。可选:
#'   \describe{
#'     \item{"prevalence"}{仅按流行率过滤}
#'     \item{"abundance"}{仅按丰度过滤}
#'     \item{"variance"}{仅按方差过滤}
#'     \item{"combined"}{同时满足流行率 + 丰度条件 (推荐)}
#'     \item{"liberal"}{宽松模式: 流行率 OR 丰度满足其一即保留}
#'     \item{"custom"}{仅使用用户指定的条件}
#'   }
#' @param prev_threshold 流行率阈值，默认 0.1 (10%)。
#'   物种在至少该比例的样本中需检出 (计数 > 0)。
#'   常用范围: 0.05-0.2，文献中 10% 是使用最多的阈值。
#' @param abund_threshold 平均相对丰度阈值，默认 1e-5 (0.001%)。
#'   物种跨样本的平均相对丰度需高于该值。
#'   注意: 该阈值始终基于相对丰度 (比例)，函数内部会自动转换。
#' @param min_count 最小计数阈值，默认 2。
#'   物种在至少 min_count_n 个样本中的计数需 >= min_count。
#'   仅在输入为计数数据时有意义。
#' @param min_count_n 满足最小计数的样本数阈值，默认 NULL。
#'   NULL 时自动设为 prev_threshold * 样本数 (与流行率一致)。
#'   也可直接指定一个整数。
#' @param var_quantile 方差分位数阈值，默认 0.1 (移除方差最低的 10%)。
#'   基于物种在所有样本中的方差，过滤方差极低的物种。
#'   仅在 method = "variance" 或 "custom" 且 filter_by_variance = TRUE 时生效。
#' @param max_one_group_zero 最大允许的"某组全零"比例，默认 NULL (不启用)。
#'   设为如 0.5 时，若一个物种在超过 50% 的分组中全部为零则移除。
#'   需要提供 group 参数。
#' @param group 分组向量，默认 NULL。仅在使用分组相关过滤条件时需要。
#' @param filter_unclassified 是否移除未分类物种，默认 TRUE。
#'   移除行名中包含 "unclassified", "Unclassified", "uncultured",
#'   "Unknown", "unknown", "norank", "unidentified" 等的物种。
#' @param unclassified_patterns 自定义未分类物种的匹配模式，默认 NULL 使用内置模式。
#' @param filter_by_variance 是否启用方差过滤，默认 FALSE。
#'   设为 TRUE 时在 combined/custom 模式中额外应用方差过滤。
#' @param verbose 是否打印过滤信息，默认 TRUE。
#'
#' @return list 包含:
#'   \describe{
#'     \item{filtered}{过滤后的丰度矩阵 (data.frame)}
#'     \item{removed}{被移除的物种名向量}
#'     \item{stats}{过滤统计信息 (data.frame)，含每个物种的流行率、丰度等指标}
#'     \item{summary}{过滤结果摘要}
#'   }
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(c(rpois(100, 50), rpois(100, 2), rep(0, 100)),
#'               nrow = 30, ncol = 10,
#'               dimnames = list(paste0("Sp", 1:30), paste0("S", 1:10)))
#'
#' # 推荐默认参数 (combined: 流行率 10% + 相对丰度 0.001%)
#' res <- filter_abundance(mat)
#'
#' # 仅按流行率过滤
#' res <- filter_abundance(mat, method = "prevalence", prev_threshold = 0.2)
#'
#' # 宽松过滤
#' res <- filter_abundance(mat, method = "liberal")
#'
#' @export

filter_abundance <- function(
  abund_table,
  method = c(
    "combined",
    "prevalence",
    "abundance",
    "variance",
    "liberal",
    "custom"
  ),
  prev_threshold = 0.1,
  abund_threshold = 1e-5,
  min_count = 2,
  min_count_n = NULL,
  var_quantile = 0.1,
  max_one_group_zero = NULL,
  group = NULL,
  filter_unclassified = TRUE,
  unclassified_patterns = NULL,
  filter_by_variance = FALSE,
  verbose = TRUE
) {
  method <- match.arg(method)
  abund_table <- as.data.frame(abund_table)
  n_features <- nrow(abund_table)
  n_samples <- ncol(abund_table)

  if (is.null(rownames(abund_table))) {
    rownames(abund_table) <- paste0("Feature_", seq_len(n_features))
  }

  # ======================== 判断数据类型 ========================
  max_val <- max(abund_table, na.rm = TRUE)
  is_count <- all(abund_table == round(abund_table), na.rm = TRUE) &&
    max_val > 1
  data_type <- if (is_count) "counts" else "relative abundance"

  # 计算相对丰度 (用于统一评估丰度阈值)
  if (is_count) {
    lib_sizes <- colSums(abund_table)
    lib_sizes[lib_sizes == 0] <- 1
    rel_table <- sweep(abund_table, 2, lib_sizes, "/")
  } else {
    rel_table <- abund_table
  }

  # min_count_n 默认值
  if (is.null(min_count_n)) {
    min_count_n <- max(1, ceiling(prev_threshold * n_samples))
  }

  # ======================== 计算各指标 ========================
  stats <- data.frame(
    Feature = rownames(abund_table),
    # 流行率: 检出比例 (count > 0 的样本比例)
    Prevalence = rowMeans(abund_table > 0, na.rm = TRUE),
    # 检出样本数
    N_detected = rowSums(abund_table > 0, na.rm = TRUE),
    # 平均相对丰度
    Mean_rel_abund = rowMeans(rel_table, na.rm = TRUE),
    # 最大相对丰度
    Max_rel_abund = apply(rel_table, 1, max, na.rm = TRUE),
    # 平均计数 (原始数据)
    Mean_count = rowMeans(abund_table, na.rm = TRUE),
    # 总计数
    Total_count = rowSums(abund_table, na.rm = TRUE),
    # 方差 (相对丰度)
    Variance = apply(rel_table, 1, var, na.rm = TRUE),
    # 变异系数
    CV = apply(rel_table, 1, function(x) {
      m <- mean(x, na.rm = TRUE)
      if (m == 0) {
        return(0)
      }
      sd(x, na.rm = TRUE) / m
    }),
    stringsAsFactors = FALSE
  )

  # 满足最小计数条件的样本数
  if (is_count) {
    stats$N_above_min_count <- rowSums(abund_table >= min_count, na.rm = TRUE)
  } else {
    stats$N_above_min_count <- stats$N_detected
  }

  # ======================== 各过滤条件 ========================

  # 1. 流行率过滤
  pass_prev <- stats$Prevalence >= prev_threshold

  # 2. 丰度过滤
  pass_abund <- stats$Mean_rel_abund >= abund_threshold

  # 3. 最小计数过滤
  pass_count <- stats$N_above_min_count >= min_count_n

  # 4. 方差过滤
  var_cutoff <- quantile(stats$Variance, probs = var_quantile, na.rm = TRUE)
  pass_var <- stats$Variance >= var_cutoff

  # 5. 分组零值过滤
  pass_group <- rep(TRUE, n_features)
  if (!is.null(max_one_group_zero) && !is.null(group)) {
    group <- as.factor(group)
    group_names <- levels(group)
    n_groups <- length(group_names)

    n_allzero_groups <- sapply(seq_len(n_features), function(i) {
      sum(sapply(group_names, function(g) {
        idx <- which(group == g)
        all(abund_table[i, idx] == 0)
      }))
    })
    pass_group <- (n_allzero_groups / n_groups) <= max_one_group_zero
  }

  # 6. 未分类物种过滤
  pass_classified <- rep(TRUE, n_features)
  if (filter_unclassified) {
    if (is.null(unclassified_patterns)) {
      unclassified_patterns <- c(
        "unclassified",
        "Unclassified",
        "UNCLASSIFIED",
        "uncultured",
        "Uncultured",
        "[Uu]nknown",
        "[Uu]nidentified",
        "norank",
        "no_rank",
        "Norank",
        "unresolved",
        "Unresolved",
        "incertae_sedis",
        "Incertae_[Ss]edis",
        "ambiguous_taxa",
        "metagenome",
        "^[a-z]__$", # 空的分类标签如 "g__"
        "_sp\\.$",
        "_sp$", # 不确定的种名如 "Genus_sp."
        "bacterium$" # 通用名如 "uncultured bacterium"
      )
    }
    pattern <- paste(unclassified_patterns, collapse = "|")
    pass_classified <- !grepl(pattern, rownames(abund_table))
  }

  # ======================== 组合过滤逻辑 ========================
  keep <- switch(
    method,
    "combined" = pass_prev & pass_abund & pass_classified,
    "prevalence" = pass_prev & pass_classified,
    "abundance" = pass_abund & pass_classified,
    "variance" = pass_var & pass_classified,
    "liberal" = (pass_prev | pass_abund) & pass_classified,
    "custom" = pass_classified # 仅未分类过滤, 其余由用户控制
  )

  # combined 模式额外应用最小计数过滤 (仅计数数据)
  if (method == "combined" && is_count) {
    keep <- keep & pass_count
  }

  # 可选: 额外叠加方差过滤
  if (filter_by_variance && method %in% c("combined", "custom")) {
    keep <- keep & pass_var
  }

  # 可选: 分组零值过滤
  if (!is.null(max_one_group_zero) && !is.null(group)) {
    keep <- keep & pass_group
  }

  # ======================== 执行过滤 ========================
  filtered <- abund_table[keep, , drop = FALSE]
  removed <- rownames(abund_table)[!keep]

  # 过滤统计
  stats$Pass_prevalence <- pass_prev
  stats$Pass_abundance <- pass_abund
  stats$Pass_count <- pass_count
  stats$Pass_variance <- pass_var
  stats$Pass_classified <- pass_classified
  stats$Kept <- keep

  # 分类移除原因
  stats$Remove_reason <- ""
  stats$Remove_reason[!pass_classified] <- "unclassified"
  stats$Remove_reason[
    pass_classified & !pass_prev & !pass_abund
  ] <- "low_prev+low_abund"
  stats$Remove_reason[
    pass_classified & !pass_prev & pass_abund
  ] <- "low_prevalence"
  stats$Remove_reason[
    pass_classified & pass_prev & !pass_abund
  ] <- "low_abundance"
  stats$Remove_reason[keep] <- "kept"

  # ======================== 汇总信息 ========================
  n_kept <- sum(keep)
  n_removed <- sum(!keep)
  n_removed_uncl <- sum(!pass_classified)
  n_removed_prev <- sum(pass_classified & !pass_prev)
  n_removed_abund <- sum(pass_classified & !pass_abund)

  summary_info <- list(
    input_features = n_features,
    input_samples = n_samples,
    data_type = data_type,
    method = method,
    kept_features = n_kept,
    removed_features = n_removed,
    removed_unclassified = n_removed_uncl,
    prev_threshold = prev_threshold,
    abund_threshold = abund_threshold,
    min_count = if (is_count) min_count else NA,
    min_count_n = if (is_count) min_count_n else NA,
    pct_features_kept = round(n_kept / n_features * 100, 1),
    pct_reads_kept = round(sum(abund_table[keep, ]) / sum(abund_table) * 100, 2)
  )

  if (verbose) {
    cat("\n========================================\n")
    cat("      物种过滤结果\n")
    cat("========================================\n")
    cat("输入数据:   ", n_features, "个物种 ×", n_samples, "个样本\n")
    cat("数据类型:   ", data_type, "\n")
    cat("过滤策略:   ", method, "\n")
    cat("----------------------------------------\n")
    cat(
      "流行率阈值:  >=",
      prev_threshold * 100,
      "% 样本 (>=",
      ceiling(prev_threshold * n_samples),
      "个样本)\n"
    )
    cat(
      "丰度阈值:    >=",
      format(abund_threshold, scientific = TRUE),
      "(平均相对丰度)\n"
    )
    if (is_count) {
      cat(
        "最小计数:    >=",
        min_count,
        "reads 在 >=",
        min_count_n,
        "个样本中\n"
      )
    }
    if (filter_by_variance) {
      cat("方差阈值:    移除最低", var_quantile * 100, "% 方差\n")
    }
    cat("未分类过滤: ", ifelse(filter_unclassified, "是", "否"), "\n")
    cat("----------------------------------------\n")
    cat(
      "保留物种:   ",
      n_kept,
      "/",
      n_features,
      " (",
      summary_info$pct_features_kept,
      "%)\n"
    )
    cat("移除物种:   ", n_removed, "\n")
    if (n_removed_uncl > 0) {
      cat("  其中未分类:", n_removed_uncl, "\n")
    }
    cat("保留 reads:  ", summary_info$pct_reads_kept, "%\n")
    cat("========================================\n\n")
  }

  list(
    filtered = filtered,
    removed = removed,
    stats = stats,
    summary = summary_info
  )
}


# ============================================================================
# 便捷封装: 常用预设方案
# ============================================================================

#' 严格过滤 (适合发表级分析)
#' 流行率 >= 20%, 平均相对丰度 >= 0.01%, 最小计数 >= 5 in >= 20% 样本
#' @export
filter_strict <- function(abund_table, ...) {
  filter_abundance(
    abund_table,
    method = "combined",
    prev_threshold = 0.2,
    abund_threshold = 1e-4,
    min_count = 5,
    filter_by_variance = TRUE,
    var_quantile = 0.1,
    ...
  )
}

#' 宽松过滤 (适合探索性分析)
#' 流行率 >= 5%, 平均相对丰度 >= 0.0001%
#' @export
filter_lenient <- function(abund_table, ...) {
  filter_abundance(
    abund_table,
    method = "combined",
    prev_threshold = 0.05,
    abund_threshold = 1e-6,
    min_count = 1,
    ...
  )
}

#' 差异分析专用过滤
#' 流行率 >= 10%, 最小计数 >= 2 in >= 10% 样本, 不过滤未分类物种
#' @export
filter_for_difftest <- function(abund_table, ...) {
  filter_abundance(
    abund_table,
    method = "combined",
    prev_threshold = 0.1,
    abund_threshold = 1e-5,
    min_count = 2,
    filter_unclassified = FALSE,
    ...
  )
}

#' 网络分析专用过滤 (较严格, 减少稀疏性)
#' 流行率 >= 30%, 方差过滤移除最低 20%
#' @export
filter_for_network <- function(abund_table, ...) {
  filter_abundance(
    abund_table,
    method = "combined",
    prev_threshold = 0.3,
    abund_threshold = 1e-4,
    filter_by_variance = TRUE,
    var_quantile = 0.2,
    ...
  )
}


# ============================================================================
# 可视化: 过滤阈值探索
# ============================================================================

#' 绘制流行率-丰度分布图, 帮助选择过滤阈值
#'
#' @param abund_table 丰度矩阵 (行=物种, 列=样本)
#' @param prev_threshold 流行率阈值线 (红色虚线)
#' @param abund_threshold 丰度阈值线 (蓝色虚线)
#' @return ggplot
#' @export
plot_filter_threshold <- function(
  abund_table,
  prev_threshold = 0.1,
  abund_threshold = 1e-5
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("需要 ggplot2: install.packages('ggplot2')")
  }

  abund_table <- as.data.frame(abund_table)

  # 计算相对丰度
  max_val <- max(abund_table, na.rm = TRUE)
  if (max_val > 1) {
    lib_sizes <- colSums(abund_table)
    lib_sizes[lib_sizes == 0] <- 1
    rel_table <- sweep(abund_table, 2, lib_sizes, "/")
  } else {
    rel_table <- abund_table
  }

  df <- data.frame(
    Feature = rownames(abund_table),
    Prevalence = rowMeans(abund_table > 0, na.rm = TRUE),
    Mean_abund = rowMeans(rel_table, na.rm = TRUE)
  )
  df$Status <- ifelse(
    df$Prevalence >= prev_threshold & df$Mean_abund >= abund_threshold,
    "Kept",
    "Removed"
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = Prevalence, y = Mean_abund, color = Status)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_color_manual(
      values = c("Kept" = "#2ecc71", "Removed" = "#e74c3c")
    ) +
    ggplot2::geom_vline(
      xintercept = prev_threshold,
      linetype = "dashed",
      color = "#e74c3c",
      linewidth = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = abund_threshold,
      linetype = "dashed",
      color = "#3498db",
      linewidth = 0.5
    ) +
    ggplot2::annotate(
      "text",
      x = prev_threshold + 0.02,
      y = max(df$Mean_abund) * 0.5,
      label = paste0("Prev = ", prev_threshold * 100, "%"),
      color = "#e74c3c",
      size = 3,
      hjust = 0
    ) +
    ggplot2::annotate(
      "text",
      x = 0.8,
      y = abund_threshold * 3,
      label = paste0("Abund = ", format(abund_threshold, scientific = TRUE)),
      color = "#3498db",
      size = 3
    ) +
    ggplot2::labs(
      x = "Prevalence (fraction of samples detected)",
      y = "Mean Relative Abundance (log scale)",
      title = "Feature Filtering Threshold",
      subtitle = paste0(
        "Kept: ",
        sum(df$Status == "Kept"),
        " | Removed: ",
        sum(df$Status == "Removed"),
        " / ",
        nrow(df),
        " total features"
      )
    ) +
    ggplot2::theme_bw()
  p
}


#' 绘制不同过滤阈值下保留的物种数和 reads 占比
#'
#' @param abund_table 丰度矩阵
#' @return ggplot
#' @export
plot_filter_sensitivity <- function(abund_table) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("需要 ggplot2")
  }

  abund_table <- as.data.frame(abund_table)
  n_features <- nrow(abund_table)
  total_reads <- sum(abund_table)

  # 计算相对丰度
  max_val <- max(abund_table, na.rm = TRUE)
  if (max_val > 1) {
    lib_sizes <- colSums(abund_table)
    lib_sizes[lib_sizes == 0] <- 1
    rel_table <- sweep(abund_table, 2, lib_sizes, "/")
  } else {
    rel_table <- abund_table
  }

  prev_values <- rowMeans(abund_table > 0, na.rm = TRUE)
  abund_values <- rowMeans(rel_table, na.rm = TRUE)
  read_sums <- rowSums(abund_table)

  # 流行率灵敏度
  prev_thresholds <- seq(0, 0.5, by = 0.02)
  prev_sens <- data.frame(
    Threshold = prev_thresholds,
    N_kept = sapply(prev_thresholds, function(t) sum(prev_values >= t)),
    Reads_pct = sapply(prev_thresholds, function(t) {
      sum(read_sums[prev_values >= t]) / total_reads * 100
    }),
    Filter_type = "Prevalence"
  )

  # 丰度灵敏度
  abund_thresholds <- 10^seq(-7, -2, by = 0.25)
  abund_sens <- data.frame(
    Threshold = abund_thresholds,
    N_kept = sapply(abund_thresholds, function(t) sum(abund_values >= t)),
    Reads_pct = sapply(abund_thresholds, function(t) {
      sum(read_sums[abund_values >= t]) / total_reads * 100
    }),
    Filter_type = "Abundance"
  )

  # 流行率图
  p1 <- ggplot2::ggplot(prev_sens, ggplot2::aes(x = Threshold * 100)) +
    ggplot2::geom_line(
      ggplot2::aes(y = N_kept),
      color = "#2ecc71",
      linewidth = 1
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = Reads_pct / 100 * n_features),
      color = "#3498db",
      linewidth = 1,
      linetype = "dashed"
    ) +
    ggplot2::scale_y_continuous(
      name = "Features Kept",
      sec.axis = ggplot2::sec_axis(
        ~ . / n_features * 100,
        name = "Reads Retained (%)"
      )
    ) +
    ggplot2::geom_vline(
      xintercept = 10,
      linetype = "dotted",
      color = "gray40"
    ) +
    ggplot2::labs(
      x = "Prevalence Threshold (%)",
      title = "Prevalence Filter Sensitivity",
      subtitle = "Green = features kept, Blue dashed = reads retained"
    ) +
    ggplot2::theme_bw()

  # 丰度图
  p2 <- ggplot2::ggplot(abund_sens, ggplot2::aes(x = Threshold)) +
    ggplot2::geom_line(
      ggplot2::aes(y = N_kept),
      color = "#2ecc71",
      linewidth = 1
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = Reads_pct / 100 * n_features),
      color = "#3498db",
      linewidth = 1,
      linetype = "dashed"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_continuous(
      name = "Features Kept",
      sec.axis = ggplot2::sec_axis(
        ~ . / n_features * 100,
        name = "Reads Retained (%)"
      )
    ) +
    ggplot2::geom_vline(
      xintercept = 1e-5,
      linetype = "dotted",
      color = "gray40"
    ) +
    ggplot2::labs(
      x = "Mean Relative Abundance Threshold (log scale)",
      title = "Abundance Filter Sensitivity",
      subtitle = "Green = features kept, Blue dashed = reads retained"
    ) +
    ggplot2::theme_bw()

  list(prevalence_plot = p1, abundance_plot = p2)
}
