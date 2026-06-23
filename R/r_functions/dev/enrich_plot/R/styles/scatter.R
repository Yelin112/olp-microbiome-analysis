# ==============================================================================
# plot_enrich_scatter(): 双向散点图
# 创建: 2026-04-24
#
# 布局：
#   x 轴 → Count × direction（正 = 第一组，负 = 第二组），标签显示绝对值
#   y 轴 → -log10(p.adjust)，刻度线绘制在 x = 0 右侧（自定义浮动 y 轴）
#   点   → shape 21，fill = 分组，size = Count
#   标签 → geom_text_repel，colour = 分组
#
# 依赖
# - core/std_format.R: prep_enrich_dual()
# - ggrepel
# ==============================================================================


# ------------------------------------------------------------------------------
# plot_enrich_scatter()
# ------------------------------------------------------------------------------

#' 双向散点图（Count × direction vs -log10 p.adjust）
#'
#' @param data 长度为 2 的 named list（如 list(Up = ..., Down = ...)），
#'   每个元素为 enrichResult / EnrichStd / data.frame
#' @param col_map 列名映射
#' @param top_n 每组取前 N 个 term，默认 15
#' @param filter_padj p_adjust 过滤阈值，NULL 表示不过滤
#' @param group_colors 两组点的填充色，named vector 或长度为 2 的 vector
#' @param group_labels 两组的显示名称（图例标签），NULL 时使用 list 名称
#' @param x_limits x 轴范围，NULL 时自动（±max(Count) 向上取整 × 1.2）
#' @param y_tick_step y 轴刻度间隔，NULL 时自动
#' @param size_range 点大小范围，默认 c(3, 12)
#' @param label_size 文字标签字号，默认 3
#' @param max_overlaps geom_text_repel 最大重叠数，默认 50
#' @param point_alpha 点透明度，默认 0.9
#' @param bg_color 面板背景色，默认 "#fffef5"
#' @param legend_title 分组图例标题，默认 "Group"
#' @param title 图标题，默认 NULL
#' @param base_size 基础字号，默认 11
#'
#' @return ggplot 对象
#'
#' @examples
#' \dontrun{
#' plot_enrich_scatter(
#'   list(Up = kegg_up, Down = kegg_down),
#'   top_n = 15,
#'   group_colors = c(Up = "#3369e7", Down = "#ff6c5f")
#' )
#' }
#'
#' @export
library(ggplot2)


plot_enrich_scatter <- function(
  data,
  col_map       = NULL,
  top_n         = 15,
  filter_padj   = NULL,
  group_colors  = c("#3369e7", "#ff6c5f"),
  group_labels  = NULL,
  x_var         = "gene_count",   # "gene_count" | "gene_ratio" | "rich_factor"
  y_var         = "p_adjust",     # "p_adjust" | "p_value"
  size_var      = "gene_count",   # "gene_count" | "gene_ratio" | "rich_factor" | NULL
  point_shape   = 21,
  point_size    = 4,              # size_var = NULL 时的固定点大小
  stroke_color  = "white",        # 点描边色（filled shape 时有效）
  stroke_width  = 0.3,            # 点描边宽度
  x_limits      = NULL,
  y_tick_step   = NULL,
  size_range    = c(3, 12),
  label_size    = 3,
  max_overlaps  = 50,
  point_alpha   = 0.9,
  bg_color      = "#fffef5",
  legend_title  = "Group",
  title         = NULL,
  base_size     = 11
) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("plot_enrich_scatter() 需要 ggrepel: install.packages('ggrepel')")
  }

  # 仅支持双组
  is_dual <- is.list(data) &&
    !inherits(data, c("data.frame", "enrichResult", "EnrichStd")) &&
    length(data) == 2
  if (!is_dual) {
    stop("plot_enrich_scatter() 需要长度为 2 的 named list，如 list(Up = ..., Down = ...)")
  }

  # 1. 数据准备 ----
  plot_data <- prep_enrich_dual(
    data_list   = data,
    top_n       = top_n,
    filter_padj = filter_padj,
    col_map     = col_map
  )

  group_names <- levels(plot_data$group)
  if (is.null(names(group_colors))) {
    group_colors <- setNames(group_colors[seq_along(group_names)], group_names)
  }
  if (is.null(group_labels)) {
    group_labels <- setNames(group_names, group_names)
  }

  # 2. x / y / size 列计算 ----
  # x：原始量级（不做 -log10），乘以 direction 产生正负方向
  .get_raw_col <- function(df, var, fallback = "gene_count") {
    if (var == fallback) return(df[[fallback]])
    if (var %in% colnames(df) && !all(is.na(df[[var]]))) return(df[[var]])
    warning(var, " 不可用，退回 ", fallback)
    df[[fallback]]
  }

  x_raw         <- .get_raw_col(plot_data, x_var, "gene_count")
  plot_data$x_pos <- x_raw * plot_data$direction

  # y：取 -log10
  y_col         <- if (y_var == "p_value") plot_data$p_value else plot_data$p_adjust
  plot_data$y_pos <- -log10(y_col)

  # size
  if (!is.null(size_var)) {
    plot_data$.size <- .get_raw_col(plot_data, size_var, "gene_count")
  }

  # 轴标签
  x_label <- switch(x_var,
    gene_ratio  = "GeneRatio",
    rich_factor = "Rich Factor",
    "Count"
  )
  y_label <- if (y_var == "p_value") "-log10(p.value)" else "-log10(p.adjust)"
  size_label <- if (is.null(size_var)) NULL else switch(size_var,
    gene_ratio  = "GeneRatio",
    rich_factor = "Rich Factor",
    "Count"
  )

  # 3. 坐标轴范围 ----
  max_x <- max(abs(x_raw), na.rm = TRUE)
  if (is.null(x_limits)) {
    x_lim <- if (x_var == "gene_count") {
      ceiling(max_x * 1.2 / 5) * 5
    } else {
      ceiling(max_x * 1.2 * 1000) / 1000
    }
    x_limits <- c(-x_lim, x_lim)
  }

  max_y <- ceiling(-log10(min(y_col, na.rm = TRUE)))
  if (is.null(y_tick_step)) {
    y_tick_step <- if (max_y <= 10) 2 else if (max_y <= 30) 5 else 10
  }
  y_breaks   <- seq(0, max_y, by = y_tick_step)
  tick_width <- diff(x_limits) * 0.02

  # 4. 自定义 y 轴数据（浮动在 x = 0 处）----
  axis_data <- data.frame(
    x     = 0,
    y     = y_breaks,
    xend  = tick_width,
    yend  = y_breaks,
    label = as.character(y_breaks),
    stringsAsFactors = FALSE
  )
  axis_data$label[axis_data$y == 0] <- ""

  # 5. 判断点类型：filled shape（21-25）用 fill 映射，solid 用 colour 映射 ----
  is_filled <- point_shape %in% 21:25

  if (is_filled) {
    if (!is.null(size_var)) {
      pt_layer <- geom_point(
        aes(fill = group, size = .size),
        shape = point_shape, alpha = point_alpha,
        color = stroke_color, stroke = stroke_width
      )
    } else {
      pt_layer <- geom_point(
        aes(fill = group),
        shape = point_shape, size = point_size, alpha = point_alpha,
        color = stroke_color, stroke = stroke_width
      )
    }
    fill_scale   <- scale_fill_manual(name = legend_title, values = group_colors,
                                      labels = group_labels)
    colour_scale <- scale_colour_manual(values = group_colors, guide = "none")
  } else {
    if (!is.null(size_var)) {
      pt_layer <- geom_point(
        aes(colour = group, size = .size),
        shape = point_shape, alpha = point_alpha
      )
    } else {
      pt_layer <- geom_point(
        aes(colour = group),
        shape = point_shape, size = point_size, alpha = point_alpha
      )
    }
    fill_scale   <- NULL
    colour_scale <- scale_colour_manual(name = legend_title, values = group_colors,
                                        labels = group_labels)
  }

  size_scale <- if (!is.null(size_var)) {
    scale_size_continuous(name = size_label, range = size_range)
  } else {
    NULL
  }

  # 6. 绘图 ----
  p <- ggplot(plot_data, aes(x = x_pos, y = y_pos)) +

    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +

    geom_segment(
      aes(x = x, y = y, xend = xend, yend = yend),
      data = axis_data, inherit.aes = FALSE, linewidth = 0.4
    ) +
    geom_text(
      aes(x = x, y = y, label = label),
      data = axis_data, inherit.aes = FALSE,
      hjust = 1.3, size = base_size / 3.8
    ) +

    pt_layer +

    ggrepel::geom_text_repel(
      aes(label = term_name, colour = group),
      max.overlaps = max_overlaps,
      size         = label_size,
      show.legend  = FALSE
    ) +

    labs(x = x_label, y = y_label, title = title) +

    fill_scale +
    colour_scale +
    size_scale +

    scale_x_continuous(limits = x_limits, labels = abs) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +

    theme_bw(base_size = base_size) +
    theme(
      axis.line.y      = element_blank(),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      panel.border     = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = bg_color),
      plot.title       = element_text(size = base_size + 3, face = "bold")
    )

  p
}


# 注册
if (exists("register_enrich_plot", mode = "function")) {
  register_enrich_plot(
    type        = "scatter",
    fn          = plot_enrich_scatter,
    description = "双向散点图（x = Count×direction，y = -log10 p，ggrepel 标签）",
    overwrite   = TRUE
  )
}
