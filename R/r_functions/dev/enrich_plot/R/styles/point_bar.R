# ==============================================================================
# plot_enrich_point_bar(): 分组 Bar + 折线点图（patchwork 拼图版）
# 创建: 2026-04-24
#
# 布局：
#   每个分类 = 左面板（折线+点，GeneRatio，x 轴反转）
#             + 右面板（水平条形，-log10 p.adjust，含通路名称标签）
#   多分类 → patchwork::wrap_plots() 拼为网格
#
# 依赖
# - core/std_format.R: as_enrich_std(), prep_enrich_combined()
# - patchwork
# ==============================================================================


# ------------------------------------------------------------------------------
# plot_enrich_point_bar()
# ------------------------------------------------------------------------------

#' 分组 Bar + 折线点图（patchwork 拼图版）
#'
#' 将多分类富集结果拼图展示，每个分类由两个面板组成：
#' - 左面板：GeneRatio 折线 + 点（x 轴反转，最右 = 0）
#' - 右面板：-log10(p.adjust) 水平条形 + 通路名称标签
#'
#' @param data 输入数据，支持两种形式：
#'   - **named list**：元素为 enrichResult / EnrichStd / data.frame，
#'     名称作为分类标签
#'   - **data.frame**：已合并数据，需同时指定 `ontology_col`
#' @param ontology_col 分类列名，data.frame 输入时使用
#' @param col_map 列名映射，见 as_enrich_std()
#' @param top_n 每个分类显示前 N 个 term，默认 5
#' @param category_order 分类排列顺序，NULL 时使用输入顺序
#' @param category_colors 分类配色，named vector，NULL 时自动生成
#' @param x_var_left 左侧点/线使用的 x 变量，默认 "gene_ratio"，
#'   不可用时自动降级为 "rich_factor" → "gene_count"
#' @param x_left_scale 左侧值缩放倍数。gene_ratio 默认×100（转百分比）
#' @param x_left_label 左面板 x 轴标签，NULL 时根据 x_var_left 自动生成
#' @param x_buffer 右侧 x 轴额外留白（expand mult），默认 0.1
#' @param label_size 通路名称字号，默认 3（右面板 y 轴文字）
#' @param point_size 折线点大小，默认 3
#' @param line_color 折线颜色，默认 "grey50"
#' @param line_width 折线宽度，默认 0.8
#' @param bar_alpha 条形透明度，默认 0.6
#' @param bar_width 条形宽度，默认 0.7
#' @param panel_widths 左右面板宽度比，默认 c(1, 2)
#' @param facet_ncol 列数，NULL 时自动（sqrt 取整）
#' @param title 图标题，默认 NULL
#' @param base_size 基础字号，默认 11
#'
#' @return patchwork 对象（可直接 print / ggsave）
#'
#' @examples
#' \dontrun{
#' plot_enrich_point_bar(
#'   list(BP = go_bp, CC = go_cc, MF = go_mf, KEGG = kegg),
#'   top_n = 5
#' )
#' }
#'
#' @export
library(ggplot2)


plot_enrich_point_bar <- function(
  data,
  ontology_col    = NULL,
  col_map         = NULL,
  top_n           = 5,
  category_order  = NULL,
  category_colors = NULL,
  x_var_left      = "gene_ratio",
  x_left_scale    = NULL,
  x_left_label    = NULL,
  x_buffer        = 0.1,
  label_size      = 3,
  point_size      = 3,
  line_color      = "grey50",
  line_width      = 0.8,
  bar_alpha       = 0.6,
  bar_width       = 0.7,
  panel_widths    = c(1, 2),
  facet_ncol      = NULL,
  title           = NULL,
  base_size       = 11
) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("plot_enrich_point_bar() 需要 patchwork: install.packages('patchwork')")
  }
  library(patchwork)

  # 1. 数据准备 ----
  plot_data <- prep_enrich_combined(
    data           = data,
    ontology_col   = ontology_col,
    col_map        = col_map,
    top_n          = top_n,
    category_order = category_order
  )

  # 2. 组内索引（最显著 = within_id 最大 = 图中最上方）----
  plot_data <- do.call(rbind, lapply(
    split(plot_data, plot_data$category),
    function(grp) {
      grp <- grp[order(grp$p_adjust, decreasing = TRUE), , drop = FALSE]
      grp$within_id <- seq_len(nrow(grp))
      grp
    }
  ))
  rownames(plot_data) <- NULL

  # 3. 左侧 x 变量 ----
  x_var_left <- .resolve_left_xvar(plot_data, x_var_left)

  # 4. 缩放 ----
  if (is.null(x_left_scale)) {
    x_left_scale <- if (x_var_left == "gene_ratio") 100 else 1
  }
  plot_data$.x_left <- plot_data[[x_var_left]] * x_left_scale

  if (is.null(x_left_label)) {
    x_left_label <- switch(x_var_left,
      gene_ratio  = "GeneRatio (%)",
      rich_factor = "RichFactor",
      gene_count  = "GeneCount"
    )
  }

  # 5. 分类配色 ----
  cats <- levels(plot_data$category)
  if (is.null(category_colors)) {
    default_pal <- c(
      '#1cc7d0', '#2dde98', '#ff6c5f', '#3369e7',
      '#fd8d3c', '#9e9ac8', '#41ab5d', '#80b1d3'
    )
    category_colors <- setNames(default_pal[seq_along(cats)], cats)
  }

  # 6. 公共主题 ----
  base_theme <- theme_bw(base_size = base_size) +
    theme(
      panel.grid   = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = base_size - 1)
    )

  # 7. 每个分类生成配对面板 ----
  cat_plots <- lapply(cats, function(cat) {
    d        <- plot_data[plot_data$category == cat, ]
    d        <- d[order(d$within_id), ]                # geom_line 按数据行顺序连接
    y_levels <- d$term_name
    d$yf     <- factor(d$term_name, levels = y_levels)
    col      <- category_colors[[cat]]

    # 右面板：水平条形 + 通路名称叠在条形上
    p_right <- ggplot(d, aes(y = yf, x = -log10(p_adjust))) +
      geom_col(
        fill  = col,
        alpha = bar_alpha,
        width = bar_width
      ) +
      geom_text(
        aes(x = 0, label = term_name),
        hjust      = 0,
        nudge_x    = 0.05,
        size       = label_size,
        color      = "black",
        lineheight = 0.85
      ) +
      scale_x_continuous(
        expand = expansion(mult = c(0, x_buffer))
      ) +
      labs(
        title = cat,
        x     = "-log10(p.adjust)",
        y     = NULL
      ) +
      base_theme +
      theme(
        plot.title  = element_text(hjust = 0.5, face = "bold", size = base_size + 1),
        axis.text.y = element_blank()
      )

    # 左面板：折线 + 点，x 轴反转（0 靠近右侧 y 轴）
    p_left <- ggplot(d, aes(y = yf, x = .x_left)) +
      geom_line(
        aes(group = 1),
        color     = line_color,
        linewidth = line_width,
        orientation = "y"
      ) +
      geom_point(
        color = col,
        size  = point_size
      ) +
      scale_x_reverse(
        expand = expansion(mult = c(x_buffer, 0))
      ) +
      labs(x = x_left_label, y = NULL) +
      base_theme +
      theme(axis.text.y = element_blank())

    p_left + p_right + plot_layout(widths = panel_widths)
  })

  # 8. 组合所有分类 ----
  ncol_use <- if (is.null(facet_ncol)) ceiling(sqrt(length(cat_plots))) else facet_ncol
  result   <- wrap_plots(cat_plots, ncol = ncol_use)

  if (!is.null(title)) {
    result <- result + plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(size = base_size + 3, face = "bold", hjust = 0.5)
      )
    )
  }

  result
}


# ------------------------------------------------------------------------------
# 内部工具：左侧 x 变量可用性检查与降级
# ------------------------------------------------------------------------------

.resolve_left_xvar <- function(data, x_var) {
  candidates <- c("gene_ratio", "rich_factor", "gene_count")
  candidates <- c(x_var, setdiff(candidates, x_var))

  for (v in candidates) {
    if (v %in% colnames(data) && any(!is.na(data[[v]]))) {
      if (v != x_var) {
        warning("x_var_left '", x_var, "' 全为 NA，自动切换为 '", v, "'")
      }
      return(v)
    }
  }
  stop("gene_ratio / rich_factor / gene_count 均不可用，无法绘制左侧折线点图")
}


# 注册
if (exists("register_enrich_plot", mode = "function")) {
  register_enrich_plot(
    type        = "point_bar",
    fn          = plot_enrich_point_bar,
    description = "拼图 Bar+折线点（左=GeneRatio, 右=-log10 p, patchwork）",
    overwrite   = TRUE
  )
}
