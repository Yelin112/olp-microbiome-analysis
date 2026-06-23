#' 检测最佳可视化策略
#'
#' 根据每组样本量自动选择最合适的可视化策略
#'
#' @param data 数据框
#' @param group_col 分组列名（字符串）
#' @param threshold_small 小样本量阈值，默认5（n < 此值使用纯散点）
#' @param threshold_medium 中等样本量阈值，默认10（n < 此值使用柱状图+散点）
#' @param threshold_large 大样本量阈值，默认20（n < 此值使用箱线图+散点）
#'
#' @return 包含策略名称和样本量统计的列表：
#'   \item{strategy}{策略名称: "pure_scatter", "bar_points", "boxplot_points", "pure_boxplot"}
#'   \item{n_per_group}{各组样本量的命名向量}
#'   \item{max_n}{最大组样本量}
#'   \item{median_n}{样本量中位数}
#'   \item{min_n}{最小组样本量}
#'
#' @importFrom dplyr group_by summarise pull
#' @importFrom rlang sym
#'
#' @examples
#' data <- data.frame(group = rep(c("A", "B"), each = 5), value = rnorm(10))
#' detect_strategy(data, "group")
#'
#' @export
detect_strategy <- function(data,
                             group_col,
                             threshold_small = 5,
                             threshold_medium = 10,
                             threshold_large = 20) {

  # 参数验证
  if (!is.data.frame(data)) {
    stop("data 必须是 data.frame")
  }
  if (!group_col %in% names(data)) {
    stop("列 '", group_col, "' 在数据框中不存在")
  }
  if (threshold_small >= threshold_medium) {
    stop("threshold_small (", threshold_small, ") 必须小于 threshold_medium (", threshold_medium, ")")
  }
  if (threshold_medium >= threshold_large) {
    stop("threshold_medium (", threshold_medium, ") 必须小于 threshold_large (", threshold_large, ")")
  }

  # 计算各组样本量
  n_summary <- data %>%
    dplyr::group_by(!!rlang::sym(group_col)) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  n_per_group <- setNames(n_summary$n, as.character(n_summary[[group_col]]))

  # 统计量
  max_n    <- max(n_per_group)
  median_n <- median(n_per_group)
  min_n    <- min(n_per_group)

  # 策略判断（以 max_n 为准）
  if (max_n < threshold_small) {
    strategy <- "pure_scatter"
    message("策略: 纯散点图 (max n=", max_n, " < ", threshold_small, ")")
  } else if (max_n < threshold_medium) {
    strategy <- "bar_points"
    message("策略: 柱状图+散点 (max n=", max_n, " < ", threshold_medium, ")")
  } else if (max_n < threshold_large) {
    strategy <- "boxplot_points"
    message("策略: 箱线图+散点 (max n=", max_n, " < ", threshold_large, ")")
  } else {
    strategy <- "pure_boxplot"
    message("策略: 纯箱线图 (max n=", max_n, " >= ", threshold_large, ")")
  }

  return(list(
    strategy     = strategy,
    n_per_group  = n_per_group,
    max_n        = max_n,
    median_n     = median_n,
    min_n        = min_n
  ))
}
