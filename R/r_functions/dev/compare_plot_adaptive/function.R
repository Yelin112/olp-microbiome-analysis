# 加载策略模块
source(file.path(dirname(sys.frame(1)$ofile), "strategies/detect_strategy.R"))
source(file.path(dirname(sys.frame(1)$ofile), "strategies/pure_scatter.R"))
source(file.path(dirname(sys.frame(1)$ofile), "strategies/bar_points.R"))
source(file.path(dirname(sys.frame(1)$ofile), "strategies/boxplot_points.R"))
source(file.path(dirname(sys.frame(1)$ofile), "strategies/pure_boxplot.R"))

#' 自适应组间比较可视化
#'
#' 根据每组样本量自动选择最佳可视化策略，生成发表级图形
#'
#' @param data 数据框，包含分组变量和数值变量
#' @param group_col 分组列名（字符串）
#' @param value_col 数值列名（字符串）
#' @param strategy 策略选择。"auto"（默认）自动检测；或手动指定：
#'   "pure_scatter", "bar_points", "boxplot_points", "pure_boxplot"
#' @param threshold_small 样本量阈值，n < 此值使用纯散点，默认 5
#' @param threshold_medium 样本量阈值，n < 此值使用柱状图+散点，默认 10
#' @param threshold_large 样本量阈值，n < 此值使用箱线图+散点，默认 20
#' @param use_violin 在 boxplot_points/pure_boxplot 策略中是否使用小提琴图，默认 FALSE
#' @param show_points 在 pure_boxplot 策略中是否叠加散点，默认 FALSE
#' @param palette 配色方案名称或颜色向量。内置: "NPG", "AAAS", "NEJM", "Lancet", "JCO", "JAMA", "D3"
#' @param xlab X轴标签，默认使用 group_col
#' @param ylab Y轴标签，默认使用 value_col
#' @param title 图标题，默认为空
#' @param verbose 是否打印策略信息，默认 TRUE
#'
#' @return 包含以下元素的列表：
#'   \item{plot}{ggplot 对象}
#'   \item{strategy}{使用的策略名称}
#'   \item{n_per_group}{各组样本量命名向量}
#'   \item{metadata}{包含 max_n, median_n, min_n, thresholds 的列表}
#'
#' @importFrom ggplot2 labs
#' @importFrom dplyr filter
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#'
#' @examples
#' # 示例1：小样本（自动使用纯散点）
#' data_small <- data.frame(
#'   group = rep(c("Control", "TreatA", "TreatB"), each = 3),
#'   value = rnorm(9, mean = rep(c(30, 35, 28), each = 3), sd = 2)
#' )
#' result <- compare_plot_adaptive(data_small, "group", "value")
#' print(result$plot)
#'
#' # 示例2：大样本（自动使用箱线图）
#' data_large <- data.frame(
#'   group = rep(c("Control", "TreatA", "TreatB"), each = 25),
#'   value = rnorm(75, mean = rep(c(30, 35, 28), each = 25), sd = 5)
#' )
#' result <- compare_plot_adaptive(data_large, "group", "value")
#' print(result$plot)
#'
#' # 示例3：手动指定策略
#' result <- compare_plot_adaptive(data_small, "group", "value", strategy = "bar_points")
#'
#' @export
compare_plot_adaptive <- function(data,
                                   group_col,
                                   value_col,
                                   strategy          = "auto",
                                   threshold_small   = 5,
                                   threshold_medium  = 10,
                                   threshold_large   = 20,
                                   use_violin        = FALSE,
                                   show_points       = FALSE,
                                   palette           = "NPG",
                                   xlab              = NULL,
                                   ylab              = NULL,
                                   title             = NULL,
                                   verbose           = TRUE) {

  # ============================================================================
  # 参数验证
  # ============================================================================

  if (!is.data.frame(data)) {
    stop("data 必须是 data.frame")
  }
  if (!group_col %in% names(data)) {
    stop("列 '", group_col, "' 在数据框中不存在")
  }
  if (!value_col %in% names(data)) {
    stop("列 '", value_col, "' 在数据框中不存在")
  }

  valid_strategies <- c("auto", "pure_scatter", "bar_points", "boxplot_points", "pure_boxplot")
  if (!strategy %in% valid_strategies) {
    stop("无效的策略 '", strategy, "'。可选: ", paste(valid_strategies, collapse = ", "))
  }

  # ============================================================================
  # 数据预处理
  # ============================================================================

  # 移除缺失值
  n_before <- nrow(data)
  data <- data[!is.na(data[[group_col]]) & !is.na(data[[value_col]]), ]
  n_removed <- n_before - nrow(data)
  if (n_removed > 0 && verbose) {
    message("已移除 ", n_removed, " 行缺失值")
  }
  if (nrow(data) == 0) {
    stop("移除缺失值后数据为空，请检查数据")
  }

  # 转为因子（保持原始顺序）
  if (!is.factor(data[[group_col]])) {
    data[[group_col]] <- factor(data[[group_col]], levels = unique(data[[group_col]]))
  }

  # ============================================================================
  # 策略检测
  # ============================================================================

  if (strategy == "auto") {
    detection <- detect_strategy(
      data, group_col,
      threshold_small, threshold_medium, threshold_large
    )
    active_strategy <- detection$strategy
    n_info          <- detection
  } else {
    if (verbose) message("手动指定策略: ", strategy)
    # 仍然获取样本量信息
    n_info <- detect_strategy(data, group_col,
                               threshold_small, threshold_medium, threshold_large)
    n_info$strategy <- strategy
    active_strategy <- strategy
  }

  # ============================================================================
  # 颜色配置
  # ============================================================================

  n_groups <- nlevels(data[[group_col]])
  colors   <- .get_colors(palette, n_groups)

  # ============================================================================
  # 绘图
  # ============================================================================

  p <- switch(active_strategy,
    pure_scatter   = plot_pure_scatter(data, group_col, value_col, colors),
    bar_points     = plot_bar_points(data, group_col, value_col, colors),
    boxplot_points = plot_boxplot_points(data, group_col, value_col, colors,
                                          use_violin = use_violin),
    pure_boxplot   = plot_pure_boxplot(data, group_col, value_col, colors,
                                        use_violin  = use_violin,
                                        show_points = show_points)
  )

  # ============================================================================
  # 标签
  # ============================================================================

  p <- p + ggplot2::labs(
    title = title,
    x     = xlab %||% group_col,
    y     = ylab %||% value_col
  )

  # ============================================================================
  # 返回结果
  # ============================================================================

  return(list(
    plot        = p,
    strategy    = active_strategy,
    n_per_group = n_info$n_per_group,
    metadata    = list(
      max_n      = n_info$max_n,
      median_n   = n_info$median_n,
      min_n      = n_info$min_n,
      thresholds = c(
        small  = threshold_small,
        medium = threshold_medium,
        large  = threshold_large
      )
    )
  ))
}


# ==============================================================================
# 内部辅助函数
# ==============================================================================

#' 获取调色板（内部）
#' @keywords internal
.get_colors <- function(palette, n_groups) {
  sci_palettes <- list(
    "NPG"    = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F",
                 "#8491B4","#91D1C2","#DC0000","#7E6148"),
    "AAAS"   = c("#3B4992","#EE0000","#008B45","#631879","#008280",
                 "#BB0021","#5F559B","#A20056","#808180"),
    "NEJM"   = c("#BC3C29","#0072B5","#E18727","#20854E","#7876B1",
                 "#6F99AD","#FFDC91","#EE4C97"),
    "Lancet" = c("#00468B","#ED0000","#42B540","#0099B4","#925E9F",
                 "#FDAF91","#AD002A","#ADB6B6"),
    "JCO"    = c("#0073C2","#EFC000","#868686","#CD534C","#7AA6DC",
                 "#003C67","#8F7700","#3B3B3B"),
    "JAMA"   = c("#374E55","#DF8F44","#00A1D5","#B24745","#79AF97",
                 "#6A6599","#80796B"),
    "D3"     = c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
                 "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")
  )

  if (length(palette) > 1) {
    base_cols <- palette
  } else if (palette %in% names(sci_palettes)) {
    base_cols <- sci_palettes[[palette]]
  } else if (palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    max_n     <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    base_cols <- RColorBrewer::brewer.pal(max(3, min(n_groups, max_n)), palette)
  } else {
    warning("配色方案 '", palette, "' 未找到，使用默认 NPG 配色")
    base_cols <- sci_palettes[["NPG"]]
  }

  if (n_groups > length(base_cols)) {
    colorRampPalette(base_cols)(n_groups)
  } else {
    base_cols[seq_len(n_groups)]
  }
}

#' 空值合并操作符（内部）
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a
