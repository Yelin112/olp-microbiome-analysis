# ==============================================================================
# plot_enrich_bar(): 富集分析条形图
# 创建: 2026-04-24
# ==============================================================================

# 依赖
# - core/std_format.R: as_enrich_std(), prep_enrich_data(), prep_enrich_dual()
# - ggplot2


# ------------------------------------------------------------------------------
# plot_enrich_bar()
# ------------------------------------------------------------------------------

#' 富集分析条形图
#'
#' 支持两种模式：
#' - **单组模式**：水平条形图，x 轴为 -log10(p.adjust)，颜色渐变
#' - **双向模式**：`data` 为长度 2 的命名 list，两组对称展示（geom_col 版本）
#'
#' @param data 富集分析结果，接受：
#'   - EnrichStd / enrichResult / data.frame（单组模式）
#'   - 长度为 2 的命名 list（双向模式），如 `list(Up = kegg_up, Down = kegg_down)`
#' @param col_map 列名映射，见 as_enrich_std()
#' @param top_n 每组前 N 个 term，默认 10
#' @param filter_padj p_adjust 过滤阈值，默认 NULL
#' @param terms 手动指定 term_name，仅单组模式有效
#' @param group_colors 颜色向量。
#'   - 单组：长度 2，作为渐变色两端，默认 c('#ffdbdb', '#f29c98')
#'   - 双向：长度 2，两组各一色，默认 c('#009dd3', '#f29c98')
#' @param group_labels 双向模式图例标签，NULL 时使用 data_list 的名称
#' @param legend_title 双向模式图例标题，默认 "Group"
#' @param x_label x 轴标签，默认 "-log10(p.adjust)"
#' @param title 图标题，默认 NULL
#' @param bar_width 条形宽度，默认 0.7
#' @param base_size 基础字号，默认 11
#'
#' @return ggplot 对象
#'
#' @examples
#' \dontrun{
#' # 单组
#' plot_enrich_bar(kegg_result, top_n = 15)
#'
#' # 双向
#' plot_enrich_bar(
#'   list(Up = kegg_up, Down = kegg_down),
#'   group_colors = c('#009dd3', '#f29c98')
#' )
#' }
#'
#' @export
library(ggplot2)


plot_enrich_bar <- function(
  data,
  col_map       = NULL,
  top_n         = 10,
  filter_padj   = NULL,
  terms         = NULL,
  group_colors  = NULL,
  group_labels  = NULL,
  legend_title  = "Group",
  x_var         = "p_adjust",
  x_label       = NULL,
  color_var     = NULL,   # 单组模式下渐变色映射列，NULL = 同 x_var；可选 "p_adjust","gene_ratio","gene_count"
  color_label   = NULL,   # 渐变图例标题，NULL 时自动推断
  title         = NULL,
  bar_width     = 0.7,
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

  if (is.null(group_colors)) {
    group_colors <- if (is_dual) c('#009dd3', '#f29c98') else c('#ffdbdb', '#f29c98')
  }

  # 数据准备
  if (is_dual) {
    plot_data    <- prep_enrich_dual(data, top_n = top_n,
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

  # 双向标签位置
  if (is_dual) {
    plot_data$label_x     <- ifelse(plot_data$direction == 1, -label_offset, label_offset)
    plot_data$label_hjust <- ifelse(plot_data$direction == 1, 1, 0)
  }

  # 绘图
  if (is_dual) {
    color_vals <- setNames(group_colors, levels(plot_data$group))
    p <- ggplot(plot_data, aes(x = x, y = term_name, fill = group)) +
      geom_col(width = bar_width, show.legend = TRUE) +
      geom_text(
        aes(x = label_x, label = as.character(term_name), hjust = label_hjust),
        show.legend = FALSE, colour = "black", size = 3.5
      ) +
      geom_hline(yintercept = 0, linewidth = 0.5, color = "grey40") +
      scale_fill_manual(
        name   = legend_title,
        values = color_vals,
        labels = if (!is.null(group_labels)) group_labels else waiver()
      ) +
      scale_x_continuous(
        labels = abs,
        limits = c(-1, 1) * max(abs(plot_data$x), na.rm = TRUE) * 1.15
      )
  } else {
    # 颜色映射列：默认与 x 相同；color_var 指定时独立计算
    if (is.null(color_var)) {
      plot_data$.fill_val <- plot_data$x
      c_label <- x_label
    } else {
      plot_data$.fill_val <- .compute_x_col(plot_data, color_var)
      c_label <- if (!is.null(color_label)) color_label else switch(color_var,
        gene_ratio  = "GeneRatio",
        gene_count  = "Gene Count",
        "-log10(p.adjust)"
      )
    }
    p <- ggplot(plot_data, aes(x = x, y = term_name, fill = .fill_val)) +
      geom_col(width = bar_width) +
      scale_fill_gradientn(
        colours = colorRampPalette(group_colors)(50),
        name    = c_label
      )
  }

  # 主题
  base_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid   = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line.x  = ggplot2::element_line(color = "black"),
      axis.text    = ggplot2::element_text(size = base_size),
      axis.title   = ggplot2::element_text(size = base_size + 2),
      plot.title   = ggplot2::element_text(size = base_size + 4, face = "bold"),
      plot.margin  = ggplot2::margin(t = 10, r = 10, b = 20, l = 20),
      legend.text  = ggplot2::element_text(size = base_size - 1)
    )

  if (is_dual) {
    base_theme <- base_theme +
      ggplot2::theme(
        axis.text.y  = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  }

  p + ggplot2::labs(x = x_label, y = NULL, title = title) + base_theme
}


# 注册
if (exists("register_enrich_plot", mode = "function")) {
  register_enrich_plot(
    type        = "bar",
    fn          = plot_enrich_bar,
    description = "条形图（支持单组渐变色 / 双向发散两组对比）",
    overwrite   = TRUE
  )
}
