# ============================================================================
# compare_plot_adaptive — 自适应组间比较可视化
# 版本: v3 (管线式架构)
# 更新日期: 2025-04-05 → 2026-06-26 (v3 重构)
# ============================================================================
#
# 架构说明:
#   策略函数只返回 geom 层列表，统计检验 / 分面 / 主题等共享功能由统一管线叠加。
#   新增策略只需写一个 layer builder 文件并在 switch() 中注册即可继承全部共享功能。

library(ggplot2)
library(ggpubr)
library(dplyr)

# ---- 加载工具模块（如调用方已加载则跳过，避免路径解析问题）----
if (!exists("%||%")) {
  UTILS_DIR <- file.path(
    dirname(sys.frame(1)$ofile), "..", "..", "..", "utils"
  )
  source(file.path(UTILS_DIR, "helpers.R"))         # %||%, %ni%
}
if (!exists("get_colors")) {
  if (!exists("UTILS_DIR")) {
    UTILS_DIR <- file.path(
      dirname(sys.frame(1)$ofile), "..", "..", "..", "utils"
    )
  }
  source(file.path(UTILS_DIR, "palette_system.R"))  # get_colors(), scale_fill_pub_d()
}

# ---- 加载策略模块 ----

STRAT_DIR <- file.path(dirname(sys.frame(1)$ofile), "strategies")
source(file.path(STRAT_DIR, "detect_strategy.R"))
source(file.path(STRAT_DIR, "layers_pure_scatter.R"))
source(file.path(STRAT_DIR, "layers_bar_points.R"))
source(file.path(STRAT_DIR, "layers_boxplot_points.R"))
source(file.path(STRAT_DIR, "layers_pure_boxplot.R"))


# ============================================================================
# 主函数
# ============================================================================

#' 自适应组间比较可视化（v3 管线式架构）
#'
#' @description
#' 根据每组样本量自动选择最佳可视化策略，输出发表级图形。
#' 支持统计检验、分面、斑马纹背景、自定义主题等共享功能 —
#' 这些功能对所有策略生效，无需在每个策略文件中重复实现。
#'
#' @param data 数据框，包含分组变量和数值变量。
#' @param group_col 分组列名（字符串）。
#' @param value_col 数值列名（字符串）。
#' @param strategy 策略选择。\code{"auto"}（默认）自动检测样本量后选择；
#'   或手动指定: \code{"pure_scatter"}, \code{"bar_points"},
#'   \code{"boxplot_points"}, \code{"pure_boxplot"}。
#' @param threshold_small 纯散点策略上限（n < 此值），默认 5。
#' @param threshold_medium 柱状图+散点策略上限（n < 此值），默认 10。
#' @param threshold_large 箱线图+散点策略上限（n < 此值），默认 20。
#' @param use_violin 在 boxplot 策略中是否使用小提琴图替代箱线图，默认 FALSE。
#' @param show_points 在 pure_boxplot 策略中是否叠加散点，默认 FALSE。
#' @param palette 配色方案。预设名（\code{"NPG"}, \code{"AAAS"},
#'   \code{"NEJM"}, \code{"Lancet"}, \code{"JCO"}, \code{"JAMA"}, \code{"D3"}）
#'   或自定义颜色向量。默认 \code{"NPG"}。
#' @param add_bg 是否添加斑马纹背景。默认 FALSE。
#' @param bg_color 斑马纹背景颜色。默认 \code{"#0000000D"}。
#' @param add_stat 统计检验方法: \code{"none"}, \code{"t.test"},
#'   \code{"wilcox.test"}, \code{"anova"}, \code{"kruskal.test"}。
#' @param stat_label 统计标签格式: \code{"p.signif"}（星号）,
#'   \code{"p.format"}（"p = 0.001"）。
#' @param comparisons 指定比较对列表，如 \code{list(c("A","B"), c("A","C"))}。
#'   若为 NULL 且方法为两两比较且组数>2，自动生成全部组合（限前10组）。
#' @param hide_ns 是否隐藏非显著比较。默认 FALSE。
#' @param step_increase 显著性标注垂直间距增量。默认 0.12。
#' @param y_expand Y轴上方扩展比例，为显著性标注留空间。默认 0.15。
#' @param split.by 分面列名（字符串），用于 \code{facet_wrap(~ split.by)}。
#' @param theme_use 主题函数。默认 \code{ggplot2::theme_classic}。
#' @param xlab X轴标签，默认使用 group_col。
#' @param ylab Y轴标签，默认使用 value_col。
#' @param title 图标题。
#' @param verbose 是否打印策略信息。默认 TRUE。
#' @param ... 传递给 \code{theme_use} 的额外参数（如 \code{base_size}）。
#'
#' @return 包含以下元素的列表:
#'   \item{plot}{ggplot 对象}
#'   \item{strategy}{使用的策略名称}
#'   \item{n_per_group}{各组样本量命名向量}
#'   \item{metadata}{包含 max_n, median_n, min_n, thresholds 的列表}
#'
#' @importFrom ggplot2 ggplot aes labs scale_fill_manual
#' @importFrom ggplot2 theme_classic facet_wrap scale_y_continuous expansion
#' @importFrom ggplot2 geom_rect element_blank theme
#' @importFrom dplyr filter
#' @importFrom ggpubr stat_compare_means
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' \dontrun{
#' # 小样本自动使用纯散点
#' data_small <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 3),
#'   value = rnorm(9)
#' )
#' result <- compare_plot_adaptive(data_small, "group", "value")
#' print(result$plot)
#'
#' # 添加统计检验
#' result <- compare_plot_adaptive(
#'   data_small, "group", "value",
#'   add_stat = "t.test",
#'   comparisons = list(c("A", "B"), c("A", "C"))
#' )
#' }
#' @export
compare_plot_adaptive <- function(
    data,
    group_col,
    value_col,
    strategy          = c("auto", "pure_scatter", "bar_points",
                           "boxplot_points", "pure_boxplot"),
    threshold_small   = 5,
    threshold_medium  = 10,
    threshold_large   = 20,
    use_violin        = FALSE,
    show_points       = FALSE,
    palette           = "NPG",
    add_bg            = FALSE,
    bg_color          = "#0000000D",
    add_stat          = c("none", "t.test", "wilcox.test",
                           "anova", "kruskal.test"),
    stat_label        = c("p.signif", "p.format"),
    comparisons       = NULL,
    hide_ns           = FALSE,
    step_increase     = 0.12,
    y_expand          = 0.15,
    split.by          = NULL,
    theme_use         = ggplot2::theme_classic,
    xlab              = NULL,
    ylab              = NULL,
    title             = NULL,
    verbose           = TRUE,
    ...
) {
  # ============================================================================
  # 管线阶段 1: 参数校验
  # ============================================================================

  strategy   <- match.arg(strategy)
  add_stat   <- match.arg(add_stat)
  stat_label <- match.arg(stat_label)

  if (!is.data.frame(data)) stop("data 必须是 data.frame")
  if (!group_col %in% names(data))
    stop("列 '", group_col, "' 在数据框中不存在")
  if (!value_col %in% names(data))
    stop("列 '", value_col, "' 在数据框中不存在")
  if (!is.null(split.by) && !split.by %in% names(data))
    stop("分面列 '", split.by, "' 在数据框中不存在")

  # ============================================================================
  # 管线阶段 2: 数据预处理
  # ============================================================================

  n_before <- nrow(data)
  data <- data[!is.na(data[[group_col]]) & !is.na(data[[value_col]]), ]
  n_removed <- n_before - nrow(data)
  if (n_removed > 0 && verbose) {
    message("已移除 ", n_removed, " 行缺失值")
  }
  if (nrow(data) == 0) stop("移除缺失值后数据为空，请检查数据")

  if (!is.null(split.by)) {
    data <- data[!is.na(data[[split.by]]), ]
  }

  # 转为因子（保持原始顺序）
  if (!is.factor(data[[group_col]])) {
    data[[group_col]] <- factor(data[[group_col]],
                                 levels = unique(data[[group_col]]))
  }

  # ============================================================================
  # 管线阶段 3: 策略检测
  # ============================================================================

  if (strategy == "auto") {
    detection <- detect_strategy(data, group_col,
                                  threshold_small, threshold_medium,
                                  threshold_large)
    active_strategy <- detection$strategy
    if (verbose) {
      message("策略: ", active_strategy,
              " (max n per group = ", detection$max_n, ")")
    }
  } else {
    if (verbose) message("手动指定策略: ", strategy)
    detection <- detect_strategy(data, group_col,
                                  threshold_small, threshold_medium,
                                  threshold_large)
    detection$strategy <- strategy
    active_strategy <- strategy
  }

  # ============================================================================
  # 管线阶段 4: 颜色配置
  # ============================================================================

  n_groups <- nlevels(data[[group_col]])
  colors <- get_colors(palette, n = n_groups, type = "discrete")

  # ============================================================================
  # 管线阶段 5: 初始化 ggplot — 建立统一 aes 映射
  # ============================================================================

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x    = .data[[group_col]],
      y    = .data[[value_col]],
      fill = .data[[group_col]]
    )
  )

  # ============================================================================
  # 管线阶段 6: 背景层
  # ============================================================================

  if (add_bg) {
    p <- p + .add_zebra_background(data[[group_col]], bg_color)
  }

  # ============================================================================
  # 管线阶段 7: 策略图层 — 每个策略返回 geom 列表，统一用 + 叠加
  # ============================================================================

  p <- p + switch(
    active_strategy,
    pure_scatter   = build_layers_pure_scatter(),
    bar_points     = build_layers_bar_points(),
    boxplot_points = build_layers_boxplot_points(use_violin = use_violin),
    pure_boxplot   = build_layers_pure_boxplot(
      use_violin  = use_violin,
      show_points = show_points
    )
  )

  # ============================================================================
  # 管线阶段 8: 统计检验层
  # ============================================================================

  if (add_stat != "none") {
    p <- p + .add_stat_test(
      data          = data,
      group_col     = group_col,
      method        = add_stat,
      label         = stat_label,
      comparisons   = comparisons,
      hide_ns       = hide_ns,
      step_increase = step_increase
    )
  }

  # ============================================================================
  # 管线阶段 9: Y轴扩展（为显著性标注留空间）
  # ============================================================================

  if (add_stat != "none") {
    p <- p + ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, y_expand))
    )
  }

  # ============================================================================
  # 管线阶段 10: 分面
  # ============================================================================

  if (!is.null(split.by)) {
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste("~", split.by)),
      scales = "free"
    )
  }

  # ============================================================================
  # 管线阶段 11: 主题与标签（统一收官）
  # ============================================================================

  p <- p +
    ggplot2::scale_fill_manual(values = colors, name = group_col) +
    ggplot2::labs(
      title = title,
      x     = xlab %||% group_col,
      y     = ylab %||% value_col
    )

  if (is.function(theme_use)) {
    p <- p + theme_use(...)
  }

  if (add_bg) {
    p <- p + ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
  }

  # ============================================================================
  # 管线阶段 12: 组装返回列表
  # ============================================================================

  return(list(
    plot         = p,
    strategy     = active_strategy,
    n_per_group  = detection$n_per_group,
    metadata     = list(
      max_n      = detection$max_n,
      median_n   = detection$median_n,
      min_n      = detection$min_n,
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


#' 添加斑马纹背景（内部）
#' @keywords internal
.add_zebra_background <- function(group_var, bg_color) {
  n_x_groups <- nlevels(as.factor(group_var))
  bg_df <- data.frame(
    xmin = seq(1, n_x_groups, 2) - 0.5,
    xmax = seq(1, n_x_groups, 2) + 0.5
  )
  ggplot2::geom_rect(
    data        = bg_df,
    mapping     = ggplot2::aes(xmin = xmin, xmax = xmax,
                               ymin = -Inf, ymax = Inf),
    fill        = bg_color,
    inherit.aes = FALSE,
    show.legend = FALSE
  )
}


#' 添加统计检验（内部）
#' @keywords internal
.add_stat_test <- function(data, group_col, method, label,
                            comparisons, hide_ns, step_increase) {
  n_groups <- nlevels(data[[group_col]])

  # 自动生成两两比较（仅 t.test / wilcox.test，组数>2 且未手动指定）
  if (is.null(comparisons) && n_groups > 2 &&
      method %in% c("t.test", "wilcox.test")) {
    group_levels <- levels(data[[group_col]])
    comparisons <- utils::combn(group_levels, 2, simplify = FALSE)
    if (length(comparisons) > 10) {
      warning("组数较多（>5组），仅显示前 10 个比较。建议手动指定 comparisons。")
      comparisons <- comparisons[1:10]
    }
  }

  # 整体检验（ANOVA/Kruskal）无 pairwise 时用 p.format 显示具体 p 值
  label_format <- if (is.null(comparisons) && n_groups > 2) "p.format" else label

  ggpubr::stat_compare_means(
    comparisons   = comparisons,
    method        = method,
    label         = label_format,
    hide.ns       = hide_ns,
    step.increase = step_increase,
    size          = 3.5,
    bracket.size  = 0.5,
    tip.length    = 0.02
  )
}



