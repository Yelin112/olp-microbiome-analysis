# ==============================================================================
# plot_enrich_lollipop(): 富集分析棒棒糖图
# 创建: 2026-04-24
# ==============================================================================

# 依赖
# - core/std_format.R: as_enrich_std(), prep_enrich_data(), prep_enrich_dual()
# - ggplot2


# ------------------------------------------------------------------------------
# 内部核心：构建发散棒棒糖图层（单/双组通用）
# ------------------------------------------------------------------------------

library(ggplot2)


.build_lollipop <- function(
  plot_data,          # 已处理好的 data.frame，含 x, term_name, group, direction
  group_colors,       # 长度为 1（单组）或 2（双组）的颜色向量
  group_labels,       # 图例标签，NULL 时用 group_colors 的 names
  size_var,           # 映射到点大小的列名，NULL 时等大
  size_range,         # 点大小范围
  label_offset,       # 标签与零线的水平偏移量（自动计算，可覆盖）
  point_size,         # size_var = NULL 时的固定点大小
  line_width,
  legend_title,
  x_label,
  is_dual             # 是否双向发散模式
) {
  # 标签位置：在零线一侧偏移
  if (is_dual) {
    plot_data$label_x <- ifelse(
      plot_data$direction == 1, -label_offset, label_offset
    )
    plot_data$label_hjust <- ifelse(plot_data$direction == 1, 1, 0)
  }

  # 基础映射
  if (is_dual) {
    aes_base <- aes(x = x, y = term_name, colour = group)
  } else {
    aes_base <- aes(x = x, y = term_name, colour = x)
  }

  p <- ggplot(plot_data, aes_base) +
    geom_segment(
      aes(xend = 0, yend = term_name),
      linewidth  = line_width,
      show.legend = FALSE
    )

  # 点：按 size_var 映射或固定大小
  if (!is.null(size_var) && size_var %in% colnames(plot_data)) {
    p <- p + geom_point(aes(size = .data[[size_var]]))
  } else {
    p <- p + geom_point(size = point_size)
  }

  # 双向模式：y 轴标签改用 geom_text
  if (is_dual) {
    p <- p +
      geom_text(
        aes(x = label_x, label = as.character(term_name), hjust = label_hjust),
        show.legend = FALSE,
        colour = "black",
        size = 3.5
      ) +
      geom_hline(yintercept = 0, linewidth = 0.5, color = "grey40")
  }

  # 配色
  if (is_dual) {
    color_vals <- setNames(group_colors, levels(plot_data$group))
    p <- p + scale_colour_manual(
      name   = legend_title,
      values = color_vals,
      labels = if (!is.null(group_labels)) group_labels else waiver()
    )
  } else {
    p <- p + scale_colour_gradientn(
      colours = colorRampPalette(group_colors)(50),
      name    = "-log10(p.adjust)"
    )
  }

  # 大小图例
  if (!is.null(size_var) && size_var %in% colnames(plot_data)) {
    size_name <- switch(size_var,
      gene_count  = "Gene Count",
      gene_ratio  = "Gene Ratio",
      rich_factor = "Rich Factor",
      size_var
    )
    p <- p + scale_size_continuous(range = size_range, name = size_name)
  }

  # x 轴：双向时显示绝对值
  if (is_dual) {
    x_lim <- max(abs(plot_data$x), na.rm = TRUE) * 1.15
    p <- p + scale_x_continuous(
      labels = abs,
      limits = c(-x_lim, x_lim)
    )
  }

  p
}


# ------------------------------------------------------------------------------
# plot_enrich_lollipop()
# ------------------------------------------------------------------------------

#' 富集分析棒棒糖图
#'
#' 支持两种模式：
#' - **单组模式**：`data` 为单个富集结果，x 轴为 -log10(p.adjust)，颜色渐变
#' - **双向模式**：`data` 为长度 2 的命名 list，两组对称展示，
#'   x 轴正负分别代表两组（参考模板样式）
#'
#' @param data 富集分析结果，接受：
#'   - EnrichStd / enrichResult / data.frame（单组模式）
#'   - 长度为 2 的命名 list（双向模式），如 `list(Up = kegg_up, Down = kegg_down)`
#' @param col_map  列名映射，见 as_enrich_std()
#' @param top_n    每组显示前 N 个 term（双向模式为每组各 top_n），默认 10
#' @param filter_padj p_adjust 过滤阈值，默认 NULL
#' @param terms    手动指定 term_name，仅单组模式有效
#' @param group_colors 颜色向量。
#'   - 单组：长度 2，作为渐变色的两端，默认 c('#ffdbdb', '#f29c98')
#'   - 双向：长度 2，两组各一色，默认 c('#009dd3', '#f29c98')
#' @param group_labels 双向模式图例标签，NULL 时使用 data_list 的名称
#' @param size_var  映射到点大小的列，可选 "gene_count"、"gene_ratio"、
#'   "rich_factor"，NULL 时所有点等大
#' @param size_range 点大小范围，默认 c(2, 8)
#' @param point_size size_var = NULL 时的固定点大小，默认 4
#' @param line_width 线段宽度，默认 1.2
#' @param legend_title 双向模式的图例标题，默认 "Group"
#' @param x_label x 轴标签，默认 "-log10(p.adjust)"
#' @param title 图标题，默认 NULL
#' @param base_size 基础字号，默认 11
#'
#' @return ggplot 对象
#'
#' @examples
#' \dontrun{
#' # 单组
#' plot_enrich_lollipop(kegg_result, top_n = 15)
#'
#' # 双向（两组对比）
#' plot_enrich_lollipop(
#'   list(Up = kegg_up, Down = kegg_down),
#'   group_colors = c('#009dd3', '#f29c98'),
#'   size_var = "gene_count"
#' )
#' }
#'
#' @export
plot_enrich_lollipop <- function(
  data,
  col_map       = NULL,
  top_n         = 10,
  filter_padj   = NULL,
  terms         = NULL,
  group_colors  = NULL,
  group_labels  = NULL,
  size_var      = "gene_count",
  size_range    = c(2, 8),
  point_size    = 4,
  line_width    = 1.2,
  legend_title  = "Group",
  x_var         = "p_adjust",
  x_label       = NULL,
  title         = NULL,
  base_size     = 11
) {
  if (is.null(x_label)) {
    x_label <- switch(x_var,
      gene_ratio  = "GeneRatio",
      gene_count  = "Gene Count",
      "-log10(p.adjust)"
    )
  }
  is_dual <- .is_dual_input(data)

  # 默认配色
  if (is.null(group_colors)) {
    group_colors <- if (is_dual) c('#009dd3', '#f29c98') else c('#ffdbdb', '#f29c98')
  }

  # 数据准备
  if (is_dual) {
    plot_data <- prep_enrich_dual(data, top_n = top_n,
                                  filter_padj = filter_padj, col_map = col_map,
                                  x_var = x_var)
    label_offset <- max(abs(plot_data$x), na.rm = TRUE) * 0.03
  } else {
    std <- as_enrich_std(data, col_map = col_map)
    std <- prep_enrich_data(std, top_n = top_n, filter_padj = filter_padj,
                            terms = terms)
    std$group     <- factor("single")
    std$direction <- 1L
    std$x         <- .compute_x_col(std, x_var)
    std <- std[order(std$x), ]
    std$term_name <- factor(std$term_name, levels = std$term_name)
    plot_data     <- std
    label_offset  <- NULL
  }

  # 构建图形
  p <- .build_lollipop(
    plot_data    = plot_data,
    group_colors = group_colors,
    group_labels = group_labels,
    size_var     = size_var,
    size_range   = size_range,
    label_offset = if (is_dual) label_offset else NULL,
    point_size   = point_size,
    line_width   = line_width,
    legend_title = legend_title,
    x_label      = x_label,
    is_dual      = is_dual
  )

  # 主题
  base_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid      = ggplot2::element_blank(),
      panel.border    = ggplot2::element_blank(),
      axis.line.x     = ggplot2::element_line(color = "black"),
      axis.text       = ggplot2::element_text(size = base_size),
      axis.title      = ggplot2::element_text(size = base_size + 2),
      plot.title      = ggplot2::element_text(size = base_size + 4, face = "bold"),
      plot.margin     = ggplot2::margin(t = 10, r = 10, b = 20, l = 20),
      legend.text     = ggplot2::element_text(size = base_size - 1)
    )

  if (is_dual) {
    base_theme <- base_theme +
      ggplot2::theme(
        axis.text.y  = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  }

  p <- p +
    ggplot2::labs(x = x_label, y = NULL, title = title) +
    base_theme

  p
}


# 注册
if (exists("register_enrich_plot", mode = "function")) {
  register_enrich_plot(
    type        = "lollipop",
    fn          = plot_enrich_lollipop,
    description = "棒棒糖图（支持单组渐变色 / 双向发散两组对比）",
    overwrite   = TRUE
  )
}
