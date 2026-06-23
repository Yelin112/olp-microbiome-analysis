library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)

#' 组间比较可视化函数 (compare_plot)
#'
#' @description
#' 高度集成的 ggplot2 包装器，用于快速绘制高质量的组间比较图。
#' 支持箱线图、小提琴图、柱状图等，内置统计检验、斑马纹背景、趋势线及智能配色。
#'
#' @param data 数据框 (data.frame)。
#' @param value.var Y轴数值列名 (字符串)。
#' @param group.by X轴分组列名 (字符串)。
#' @param fill.by 填充颜色列名 (字符串)。默认与 group.by 相同。若指定不同列，则开启"组内簇状比较"模式。
#' @param split.by 分面列名 (字符串)。用于 facet_wrap 分面。
#' @param add_bg 是否添加斑马纹背景 (逻辑值)。默认 FALSE。
#' @param bg_color 背景颜色 (字符串)。默认为 "#0000000D" (5%透明度的黑色)。
#' @param plot_type 主图类型: "violin", "box", "bar", "dot"。
#' @param add_box 是否叠加箱线图 (逻辑值)。仅在 plot_type = "violin" 时有效。
#' @param add_point 是否叠加散点 (逻辑值)。
#' @param point_shape 散点形状 (整数)。默认 21 (可填充圆点)。常用: 16(实心圆), 21(可填充圆), 22(可填充方), 23(可填充菱形)。
#' @param point_size 散点大小 (数值)。默认 2。
#' @param point_alpha 散点透明度 (0-1)。默认 0.6。
#' @param point_fill 散点填充颜色。默认 "same" (与箱体相同但更透明)，可设为 "white"、"black" 或任意颜色。
#' @param point_color 散点边框颜色 (字符串)。默认 "black"。
#' @param point_stroke 散点边框粗细 (数值)。默认 0.3。
#' @param point_jitter 散点抖动幅度 (数值)。默认 0.2，设为 0 则无抖动。
#' @param add_trend 是否叠加趋势线 (逻辑值)。
#' @param add_stat 统计方法: "none", "mean", "median", "t.test", "wilcox.test", "anova", "kruskal.test"。
#' @param stat_label 统计标签格式: "p.signif" (星号), "p.format" (p=0.001)。默认 "p.signif"。
#' @param comparisons 指定比较组 (List)。例如 list(c("A","B"), c("A","C"))。
#' @param hide_ns 是否隐藏非显著结果 (逻辑值)。默认 TRUE。
#' @param step_increase 显著性标注的垂直间距增量。默认 0.12，可调整避免重叠。
#' @param y_expand Y轴上方扩展比例，为显著性标注留出空间。默认 0.15 (15%)。
#' @param palette 配色方案。可以是预设名 ("NPG", "AAAS", "Paired"...) 或自定义颜色向量。
#' @param theme_use 使用的主题函数。例如 theme_classic, theme_bw, theme_minimal。
#' @param title,xlab,ylab 标题与轴标签。
#' @param strategy 自适应策略。"auto"（默认）根据样本量自动选择绘图类型；
#'   "none" 禁用自动选择（使用 plot_type 手动指定）；
#'   或直接指定: "pure_scatter"（n<5，散点+均值线）, "bar_points"（5≤n<10，柱+散点）,
#'   "boxplot_points"（10≤n<20，箱线图+散点）, "pure_boxplot"（n≥20，箱线图）。
#' @param threshold_small 纯散点策略的上限（n < 此值），默认 5。
#' @param threshold_medium 柱状图+散点策略的上限（n < 此值），默认 10。
#' @param threshold_large 箱线图+散点策略的上限（n < 此值），默认 20。
#' @param ... 传递给 theme_use 的其他参数 (如 base_size)。
#'
#' @return ggplot 对象
#' @export
#'
#' @examples
#' # 基础用法
#' compare_plot(iris, value.var = "Sepal.Length", group.by = "Species")
#'
#' # 添加统计检验
#' compare_plot(iris, value.var = "Sepal.Length", group.by = "Species",
#'              add_stat = "anova", comparisons = list(c("setosa", "versicolor")))
#'
#' # 组内比较（簇状图）
#' compare_plot(data, value.var = "value", group.by = "time", fill.by = "treatment",
#'              add_stat = "t.test", plot_type = "box")
#'
compare_plot <- function(
  data,
  value.var,
  group.by,
  fill.by = NULL,
  split.by = NULL,
  add_bg = FALSE,
  bg_color = "#0000000D",
  plot_type = c("violin", "box", "bar", "dot", "scatter"),
  add_box = FALSE,
  add_point = FALSE,
  point_shape = 21,
  point_size = 2,
  point_alpha = 1,
  point_fill = "same",
  point_color = "black",
  point_stroke = 0.7,
  point_jitter = 0.2,
  add_trend = FALSE,
  add_stat = c(
    "none",
    "mean",
    "median",
    "t.test",
    "wilcox.test",
    "anova",
    "kruskal.test"
  ),
  stat_label = c("p.signif", "p.format"),
  comparisons = NULL,
  hide_ns = FALSE,
  step_increase = 0.12,
  y_expand = 0.15,
  palette = "NPG",
  xlab = NULL,
  ylab = NULL,
  title = NULL,
  theme_use = theme_classic,
  strategy = c(
    "auto",
    "none",
    "pure_scatter",
    "bar_points",
    "boxplot_points",
    "pure_boxplot"
  ),
  threshold_small = 5,
  threshold_medium = 10,
  threshold_large = 20,
  box_args = list(),
  violin_args = list(),
  bar_args = list(),
  errorbar_args = list(),
  crossbar_args = list(),
  ...
) {
  # ============================================================================
  # 参数校验与预处理
  # ============================================================================

  plot_type <- match.arg(plot_type)
  add_stat <- match.arg(add_stat)
  stat_label <- match.arg(stat_label)
  strategy <- match.arg(strategy)

  # 设置默认标签
  xlab <- xlab %||% group.by
  ylab <- ylab %||% value.var

  # 确定填充变量
  fill_var <- fill.by %||% group.by

  # 检查必需列是否存在
  cols_needed <- unique(c(value.var, group.by, fill_var, split.by))
  cols_needed <- cols_needed[!is.null(cols_needed)]

  if (!all(cols_needed %in% colnames(data))) {
    missing_cols <- setdiff(cols_needed, colnames(data))
    stop("以下列名在数据框中不存在: ", paste(missing_cols, collapse = ", "))
  }

  # 数据预处理：移除缺失值
  data <- data %>%
    filter(!is.na(.data[[value.var]]) & !is.na(.data[[group.by]]))

  if (nrow(data) == 0) {
    stop("移除缺失值后数据为空，请检查数据。")
  }

  # 转换为因子（保持原始顺序）
  for (col in unique(c(group.by, fill_var, split.by))) {
    if (!is.null(col) && !is.factor(data[[col]])) {
      data[[col]] <- factor(data[[col]], levels = unique(data[[col]]))
    }
  }

  # ============================================================================
  # 自适应策略解析（根据样本量自动选择绘图类型）
  # ============================================================================

  if (strategy != "none") {
    # 以 fill_var 分组计算样本量（组内比较模式下更精确）
    n_per_group <- tapply(data[[value.var]], data[[fill_var]], length)
    max_n <- max(n_per_group)

    # 自动检测时根据阈值确定策略
    if (strategy == "auto") {
      if (max_n < threshold_small) {
        strategy <- "pure_scatter"
      } else if (max_n < threshold_medium) {
        strategy <- "bar_points"
      } else if (max_n < threshold_large) {
        strategy <- "boxplot_points"
      } else {
        strategy <- "pure_boxplot"
      }
      message("自适应策略: ", strategy, " (max n per group = ", max_n, ")")
    }

    # 将策略映射到图形参数
    if (strategy == "pure_scatter") {
      plot_type <- "scatter"
      add_point <- TRUE
      point_size <- 3
      point_stroke <- 0.5
    } else if (strategy == "bar_points") {
      plot_type <- "bar"
      add_point <- TRUE
      point_fill <- "white"
      point_jitter <- 0.2
      point_size <- 3
    } else if (strategy == "boxplot_points") {
      plot_type <- "box"
      add_point <- TRUE
      point_size <- 3
    } else if (strategy == "pure_boxplot") {
      plot_type <- "box"
      # add_point 保持用户设置（默认 FALSE）
    }
  }

  # ============================================================================
  # 颜色配置
  # ============================================================================

  final_colors <- get_color_palette(palette, nlevels(data[[fill_var]]))

  # ============================================================================
  # 绘图初始化
  # ============================================================================

  p <- ggplot(
    data,
    aes(x = .data[[group.by]], y = .data[[value.var]], fill = .data[[fill_var]])
  )

  is_grouped <- (fill_var != group.by)
  dodge_width <- 0.8
  pos <- if (is_grouped) position_dodge(width = dodge_width) else "identity"

  # ============================================================================
  # 背景层
  # ============================================================================

  if (add_bg) {
    p <- p + add_zebra_background(data[[group.by]], bg_color)
  }

  # ============================================================================
  # 主图层
  # ============================================================================

  p <- p +
    add_main_geom(
      plot_type = plot_type,
      pos = pos,
      is_grouped = is_grouped,
      fill_var = fill_var,
      box_args = box_args,
      violin_args = violin_args,
      bar_args = bar_args,
      errorbar_args = errorbar_args,
      crossbar_args = crossbar_args
    )

  # ============================================================================
  # 叠加层
  # ============================================================================

  # 趋势线
  if (add_trend) {
    p <- p + add_trend_line(is_grouped, fill_var, pos)
  }

  # 叠加箱线图（仅 violin 图）
  if (add_box && plot_type == "violin") {
    p <- p + add_overlay_box(pos, fill_var)
  }

  # 散点
  if (add_point) {
    p <- p +
      add_jitter_points(
        is_grouped = is_grouped,
        dodge_width = dodge_width,
        fill_var = fill_var,
        point_shape = point_shape,
        point_size = point_size,
        point_alpha = point_alpha,
        point_fill = point_fill,
        point_color = point_color,
        point_stroke = point_stroke,
        point_jitter = point_jitter,
        final_colors = final_colors
      )
  }

  # ============================================================================
  # 统计层
  # ============================================================================

  # 统计点（均值/中位数）
  if (add_stat %in% c("mean", "median")) {
    p <- p + add_stat_point(add_stat, pos, fill_var)
  }

  # 统计检验
  if (add_stat %in% c("t.test", "wilcox.test", "anova", "kruskal.test")) {
    p <- p +
      add_stat_test(
        method = add_stat,
        label = stat_label,
        comparisons = comparisons,
        is_grouped = is_grouped,
        fill_var = fill_var,
        group.by = group.by,
        data = data,
        hide_ns = hide_ns,
        step_increase = step_increase
      )
  }

  # ============================================================================
  # Y轴扩展（为显著性标注留空间）
  # ============================================================================

  if (add_stat %in% c("t.test", "wilcox.test", "anova", "kruskal.test")) {
    y_range <- range(data[[value.var]], na.rm = TRUE)
    y_margin <- diff(y_range) * y_expand
    p <- p +
      scale_y_continuous(
        expand = expansion(mult = c(0.05, y_expand))
      )
  }

  # ============================================================================
  # 分面
  # ============================================================================

  if (!is.null(split.by)) {
    p <- p + facet_wrap(as.formula(paste("~", split.by)), scales = "free")
  }

  # ============================================================================
  # 主题与样式
  # ============================================================================

  p <- p +
    scale_fill_manual(
      values = final_colors,
      name = fill_var,
      na.value = "transparent"
    ) +
    labs(title = title, x = xlab, y = ylab)

  # 应用主题
  if (is.function(theme_use)) {
    p <- p + theme_use(...)
  }

  # 移除网格线（如有背景）
  if (add_bg) {
    p <- p +
      theme(
        panel.grid = element_blank()
      )
  }

  # 趋势线颜色（统一为黑色，避免图例混淆）
  if (is_grouped && add_trend) {
    p <- p +
      scale_color_manual(
        values = setNames(
          rep("black", nlevels(data[[fill_var]])),
          levels(data[[fill_var]])
        ),
        guide = "none"
      )
  }

  return(p)
}


# ==============================================================================
# 辅助函数
# ==============================================================================

#' 获取调色板
#' @keywords internal
get_color_palette <- function(palette, n_groups) {
  # 预设科研配色
  sci_palettes <- list(
    "NPG" = c(
      "#E64B35",
      "#4DBBD5",
      "#00A087",
      "#3C5488",
      "#F39B7F",
      "#8491B4",
      "#91D1C2",
      "#DC0000",
      "#7E6148"
    ),
    "AAAS" = c(
      "#3B4992",
      "#EE0000",
      "#008B45",
      "#631879",
      "#008280",
      "#BB0021",
      "#5F559B",
      "#A20056",
      "#808180"
    ),
    "NEJM" = c(
      "#BC3C29",
      "#0072B5",
      "#E18727",
      "#20854E",
      "#7876B1",
      "#6F99AD",
      "#FFDC91",
      "#EE4C97"
    ),
    "Lancet" = c(
      "#00468B",
      "#ED0000",
      "#42B540",
      "#0099B4",
      "#925E9F",
      "#FDAF91",
      "#AD002A",
      "#ADB6B6"
    ),
    "JCO" = c(
      "#0073C2",
      "#EFC000",
      "#868686",
      "#CD534C",
      "#7AA6DC",
      "#003C67",
      "#8F7700",
      "#3B3B3B"
    ),
    "JAMA" = c(
      "#374E55",
      "#DF8F44",
      "#00A1D5",
      "#B24745",
      "#79AF97",
      "#6A6599",
      "#80796B"
    ),
    "D3" = c(
      "#1F77B4",
      "#FF7F0E",
      "#2CA02C",
      "#D62728",
      "#9467BD",
      "#8C564B",
      "#E377C2",
      "#7F7F7F",
      "#BCBD22",
      "#17BECF"
    )
  )

  # 解析配色方案
  if (length(palette) == 1 && is.character(palette)) {
    if (palette %in% names(sci_palettes)) {
      base_cols <- sci_palettes[[palette]]
    } else if (palette %in% rownames(RColorBrewer::brewer.pal.info)) {
      max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
      base_cols <- RColorBrewer::brewer.pal(
        max(3, min(n_groups, max_n)),
        palette
      )
    } else {
      warning("配色方案 '", palette, "' 未找到，使用默认 NPG 配色。")
      base_cols <- sci_palettes[["NPG"]]
    }
  } else {
    base_cols <- palette
  }

  # 扩展颜色以匹配分组数
  if (n_groups > length(base_cols)) {
    final_colors <- colorRampPalette(base_cols)(n_groups)
  } else {
    final_colors <- base_cols[1:n_groups]
  }

  return(final_colors)
}


#' 添加斑马纹背景
#' @keywords internal
add_zebra_background <- function(group_var, bg_color) {
  n_x_groups <- nlevels(as.factor(group_var))
  bg_df <- data.frame(x_pos = seq(1, n_x_groups, 2)) %>%
    mutate(xmin = x_pos - 0.5, xmax = x_pos + 0.5)

  geom_rect(
    data = bg_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = bg_color,
    inherit.aes = FALSE,
    show.legend = FALSE
  )
}


#' 添加主图层
#' @keywords internal
add_main_geom <- function(
  plot_type,
  pos,
  is_grouped,
  fill_var,
  box_args,
  violin_args,
  bar_args,
  errorbar_args,
  crossbar_args
) {
  geom_list <- list()

  if (plot_type == "violin") {
    vln <- modifyList(
      list(
        scale = "width",
        trim = FALSE,
        alpha = 0.6,
        color = "black",
        linewidth = 0.5,
        position = pos
      ),
      violin_args
    )
    geom_list <- list(do.call(geom_violin, vln))
  } else if (plot_type == "box") {
    box <- modifyList(
      list(
        width = if (is_grouped) 0.6 else 0.7,
        alpha = 0.6,
        outlier.shape = NA,
        color = "black",
        linewidth = 0.8,
        position = pos
      ),
      box_args
    )
    geom_list <- list(
      do.call(geom_boxplot, box),
      stat_summary(
        aes(group = .data[[fill_var]]),
        fun = mean,
        geom = "point",
        shape = 21,
        size = 2.5,
        fill = "white",
        color = "black",
        stroke = 0.8,
        position = pos,
        show.legend = FALSE
      )
    )
  } else if (plot_type == "bar") {
    bar <- modifyList(
      list(
        fun = mean,
        geom = "bar",
        width = 0.7,
        alpha = 0.6,
        color = NA,
        linewidth = 0.5,
        position = pos
      ),
      bar_args
    )
    err <- modifyList(
      list(
        fun.data = mean_se,
        geom = "errorbar",
        width = 0.25,
        linewidth = 0.6,
        position = pos
      ),
      errorbar_args
    )
    geom_list <- list(do.call(stat_summary, bar), do.call(stat_summary, err))
  } else if (plot_type == "dot") {
    err <- modifyList(
      list(
        fun.data = mean_se,
        geom = "errorbar",
        width = 0.25,
        linewidth = 0.6,
        position = pos
      ),
      errorbar_args
    )
    geom_list <- list(
      stat_summary(
        fun = mean,
        geom = "point",
        size = 4,
        shape = 21,
        color = "black",
        stroke = 0.8,
        position = pos
      ),
      do.call(stat_summary, err)
    )
  } else if (plot_type == "scatter") {
    cbr <- modifyList(
      list(
        fun = mean,
        geom = "crossbar",
        width = 0.3,
        linewidth = 1.2,
        color = "black",
        fatten = 1,
        position = pos
      ),
      crossbar_args
    )
    err <- modifyList(
      list(
        fun.data = mean_se,
        geom = "errorbar",
        width = 0.18,
        linewidth = 0.8,
        color = "black",
        position = pos
      ),
      errorbar_args
    )
    geom_list <- list(do.call(stat_summary, cbr), do.call(stat_summary, err))
  }

  return(geom_list)
}


#' 添加趋势线
#' @keywords internal
add_trend_line <- function(is_grouped, fill_var, pos) {
  stat_summary(
    fun = mean,
    geom = "line",
    aes(
      group = if (is_grouped) .data[[fill_var]] else 1,
      color = if (is_grouped) .data[[fill_var]] else NULL
    ),
    linewidth = 1,
    linetype = "solid",
    position = pos,
    show.legend = FALSE
  )
}


#' 添加叠加箱线图
#' @keywords internal
add_overlay_box <- function(pos, fill_var) {
  list(
    geom_boxplot(
      width = 0.15,
      fill = "white",
      outlier.shape = NA,
      color = "black",
      linewidth = 0.5,
      alpha = 0.9,
      position = pos
    ),
    stat_summary(
      aes(group = .data[[fill_var]]),
      fun = mean,
      geom = "point",
      shape = 21,
      size = 2,
      fill = "white",
      color = "black",
      stroke = 0.7,
      position = pos,
      show.legend = FALSE
    )
  )
}


#' 添加散点
#' @keywords internal
add_jitter_points <- function(
  is_grouped,
  dodge_width,
  fill_var,
  point_shape,
  point_size,
  point_alpha,
  point_fill,
  point_color,
  point_stroke,
  point_jitter,
  final_colors
) {
  # 处理填充颜色
  if (point_fill == "same") {
    # fill 继承自分组的 aes 映射（与主图颜色一致）
    use_fill_aes <- TRUE
    fill_value <- NULL
  } else {
    # 使用指定的固定颜色
    use_fill_aes <- FALSE
    fill_value <- point_fill
  }

  # 根据形状决定是否需要填充
  # 形状 21-25 支持填充和边框，16-20 为实心
  needs_fill <- point_shape %in% 21:25

  if (is_grouped) {
    # 组内比较模式
    if (needs_fill && use_fill_aes) {
      # 使用继承的 fill，点颜色与箱体一致但更透明
      geom_point(
        aes(fill = .data[[fill_var]]), # 继承填充颜色
        position = position_jitterdodge(
          jitter.width = point_jitter,
          dodge.width = dodge_width
        ),
        shape = point_shape,
        size = point_size,
        alpha = point_alpha,
        color = point_color,
        stroke = point_stroke
      )
    } else if (needs_fill && !use_fill_aes) {
      # 使用固定填充颜色
      geom_point(
        position = position_jitterdodge(
          jitter.width = point_jitter,
          dodge.width = dodge_width
        ),
        shape = point_shape,
        size = point_size,
        alpha = point_alpha,
        fill = fill_value,
        color = point_color,
        stroke = point_stroke
      )
    } else {
      # 实心形状，不需要 fill
      geom_point(
        position = position_jitterdodge(
          jitter.width = point_jitter,
          dodge.width = dodge_width
        ),
        shape = point_shape,
        size = point_size,
        alpha = point_alpha,
        color = point_color
      )
    }
  } else {
    # 单组模式
    if (needs_fill && use_fill_aes) {
      geom_jitter(
        aes(fill = .data[[fill_var]]),
        width = point_jitter,
        shape = point_shape,
        size = point_size,
        alpha = point_alpha,
        color = point_color,
        stroke = point_stroke
      )
    } else if (needs_fill && !use_fill_aes) {
      geom_jitter(
        width = point_jitter,
        shape = point_shape,
        size = point_size,
        alpha = point_alpha,
        fill = fill_value,
        color = point_color,
        stroke = point_stroke
      )
    } else {
      geom_jitter(
        width = point_jitter,
        shape = point_shape,
        size = point_size,
        alpha = point_alpha,
        color = point_color
      )
    }
  }
}


#' 添加统计点
#' @keywords internal
add_stat_point <- function(stat_fun, pos, fill_var) {
  stat_summary(
    aes(group = .data[[fill_var]]),
    fun = stat_fun,
    geom = "point",
    shape = 23,
    size = 3.5,
    fill = "yellow",
    color = "black",
    stroke = 0.8,
    position = pos,
    show.legend = FALSE
  )
}


#' 添加统计检验
#' @keywords internal
add_stat_test <- function(
  method,
  label,
  comparisons,
  is_grouped,
  fill_var,
  group.by,
  data,
  hide_ns,
  step_increase
) {
  if (is_grouped) {
    # 组内比较模式
    stat_compare_means(
      aes(group = .data[[fill_var]]),
      method = method,
      label = label,
      hide.ns = hide_ns,
      step.increase = step_increase,
      size = 3.5,
      bracket.size = 0.5,
      tip.length = 0.02
    )
  } else {
    # 组间比较模式
    n_groups <- nlevels(data[[group.by]])

    # 自动生成比较对（如果未指定）
    if (
      is.null(comparisons) &&
        n_groups > 2 &&
        method %in% c("t.test", "wilcox.test")
    ) {
      # 生成所有两两比较
      group_levels <- levels(data[[group.by]])
      comparisons <- combn(group_levels, 2, simplify = FALSE)
      # 限制比较数量（避免过于拥挤）
      if (length(comparisons) > 10) {
        warning("组数较多，仅显示前10个比较。建议手动指定 comparisons 参数。")
        comparisons <- comparisons[1:10]
      }
    }

    # 整体检验标签（ANOVA/Kruskal）
    if (is.null(comparisons) && n_groups > 2) {
      label_format <- "p.format"
    } else {
      label_format <- label
    }

    stat_compare_means(
      comparisons = comparisons,
      method = method,
      label = label_format,
      hide.ns = hide_ns,
      step.increase = step_increase,
      size = 3.5,
      bracket.size = 0.5,
      tip.length = 0.02
    )
  }
}


#' 空值合并操作符
#' @keywords internal
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
