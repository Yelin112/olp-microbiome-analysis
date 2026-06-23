# ==============================================================================
# plot_enrich_radial(): 环形树状 + 极坐标条形图
# 创建: 2026-04-24
#
# 布局：
#   外圈 → 极坐标条形图：Count 高度，-log10(pvalue) 透明度，分类着色
#          + 通路名称文字（沿圆弧排列）
#   内圈 → 圆形树状图（ggraph dendrogram）：根 → 分类 → 通路
#          两圈通过 patchwork::inset_element 叠加
#
# 依赖
# - core/std_format.R: prep_enrich_combined()
# - ggraph, tidygraph, patchwork
# - ggnewscale (可选), stringr (可选)
# ==============================================================================


# ------------------------------------------------------------------------------
# plot_enrich_radial()
# ------------------------------------------------------------------------------

#' 环形树状 + 极坐标条形图
#'
#' @param data 输入数据（named list 或 data.frame），见 prep_enrich_combined()
#' @param ontology_col 分类列名（data.frame 输入时使用）
#' @param col_map 列名映射
#' @param top_n 每个分类显示前 N 个 term，默认 8
#' @param category_order 分类排列顺序
#' @param category_colors 分类配色（深色）named vector，NULL 时使用内置配色
#' @param category_colors_light 边的配色（浅色），NULL 时使用内置浅色配色
#' @param bar_label_size 外圈通路名称字号，默认 2
#' @param node_size 内圈分类节点大小，默认 3
#' @param node_label_size 内圈分类节点标签字号，默认 3
#' @param node_label_width 内圈标签自动折行宽度（字符数），默认 12
#' @param inner_radius 极坐标内圆半径（负值越大洞越大），NULL 时自动计算
#' @param inset_left   内嵌图形左边界，默认 0.2
#' @param inset_right  内嵌图形右边界，默认 0.8
#' @param inset_bottom 内嵌图形下边界，默认 0.2
#' @param inset_top    内嵌图形上边界，默认 0.8
#' @param title 图标题，默认 NULL
#' @param base_size 基础字号，默认 11
#'
#' @return patchwork 对象
#'
#' @examples
#' \dontrun{
#' plot_enrich_radial(
#'   list(BP = go_bp, CC = go_cc, MF = go_mf, KEGG = kegg),
#'   top_n = 8
#' )
#' }
#'
#' @export
library(ggplot2)


plot_enrich_radial <- function(
  data,
  ontology_col          = NULL,
  col_map               = NULL,
  top_n                 = 8,
  category_order        = NULL,
  category_colors       = NULL,
  category_colors_light = NULL,
  bar_label_size        = 2,
  node_size             = 3,
  node_label_size       = 3,
  node_label_width      = 12,
  inner_radius          = NULL,
  inset_left            = 0.2,
  inset_right           = 0.8,
  inset_bottom          = 0.2,
  inset_top             = 0.8,
  title                 = NULL,
  base_size             = 11
) {
  for (pkg in c("ggraph", "tidygraph", "patchwork")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("plot_enrich_radial() 需要 ", pkg, "，请先 install.packages('", pkg, "')")
    }
  }
  library(ggraph)
  library(tidygraph)
  library(patchwork)

  # 1. 数据准备 ----
  plot_data <- prep_enrich_combined(
    data           = data,
    ontology_col   = ontology_col,
    col_map        = col_map,
    top_n          = top_n,
    category_order = category_order
  )

  cats <- levels(plot_data$category)

  # 标准列 → from/to/Count/pvalue 格式
  edge_data <- data.frame(
    from   = as.character(plot_data$category),
    to     = as.character(plot_data$term_name),
    Count  = plot_data$gene_count,
    pvalue = plot_data$p_value,
    stringsAsFactors = FALSE
  )

  # 2. 分类配色 ----
  default_colors <- c(
    '#437f64', '#d88c51', '#466277', '#be5960', '#3b4368',
    '#b0897c', '#4992c2', '#a975be', '#4996a2'
  )
  default_colors_light <- c(
    '#84ae9b', '#e0b28e', '#5f707d', '#d69293', '#5e6fa3',
    '#cba7a2', '#82aecb', '#bfa0cc', '#71c2cc'
  )

  n_cats <- length(cats)
  make_pal <- function(base_pal, n) {
    if (n <= length(base_pal)) base_pal[seq_len(n)]
    else colorRampPalette(base_pal)(n)
  }

  if (is.null(category_colors)) {
    cols <- setNames(make_pal(default_colors, n_cats), cats)
  } else {
    cols <- category_colors
  }

  if (is.null(category_colors_light)) {
    cols_light <- setNames(make_pal(default_colors_light, n_cats), cats)
  } else {
    cols_light <- category_colors_light
  }

  # 3. 构建图结构（根 → 分类 → 通路）----
  edges_graph <- rbind(
    edge_data[, c("from", "to", "Count", "pvalue")],
    data.frame(
      from = "pathway", to = unique(edge_data$from),
      Count = 0, pvalue = 0,
      stringsAsFactors = FALSE
    )
  )
  edges_graph$colour <- ifelse(
    edges_graph$from == "pathway", edges_graph$to, edges_graph$from
  )

  all_nodes <- unique(c(edge_data$from, edge_data$to))
  vertices <- data.frame(
    name  = c("pathway", all_nodes),
    level = c(0L,
              ifelse(all_nodes %in% unique(edge_data$from), 1L, 2L)),
    stringsAsFactors = FALSE
  )

  # 标签折行（分类节点）
  vertices$label_wrapped <- if (requireNamespace("stringr", quietly = TRUE)) {
    stringr::str_wrap(vertices$name, node_label_width)
  } else {
    vertices$name
  }

  gr <- tbl_graph(nodes = vertices, edges = edges_graph, directed = TRUE)

  # 4. leaf 数据（极坐标条形图）----
  leaf <- vertices[vertices$level == 2, , drop = FALSE]
  leaf <- merge(leaf, edge_data, by.x = "name", by.y = "to", all.x = FALSE)
  leaf <- leaf[order(match(leaf$from, cats), leaf$name), ]   # 同类聚在一起
  leaf$name  <- factor(leaf$name, levels = leaf$name)
  leaf$id    <- seq_len(nrow(leaf))
  n_leaf     <- nrow(leaf)
  leaf$angle <- 90 - 360 * leaf$id / n_leaf
  leaf$hjust <- ifelse(leaf$angle < -90, 1, 0)
  leaf$angle <- ifelse(leaf$angle < -90, leaf$angle + 180, leaf$angle)

  max_count <- max(leaf$Count, na.rm = TRUE)
  if (is.null(inner_radius)) {
    inner_radius <- -max_count * 20
  }

  # 5. 环形树状图（内圈）----
  g <- ggraph(gr, layout = "dendrogram", circular = TRUE) +
    geom_edge_diagonal(aes(colour = colour), show.legend = FALSE) +
    geom_node_point(
      aes(filter = level == 1, colour = name),
      size = node_size, show.legend = FALSE
    ) +
    geom_node_text(
      aes(filter = level == 1, label = label_wrapped, colour = name),
      lineheight = 0.8, repel = TRUE,
      size = node_label_size, show.legend = FALSE
    ) +
    coord_fixed() +
    scale_edge_colour_manual(values = cols_light) +
    scale_colour_manual(values = cols) +
    theme_void()

  if (requireNamespace("ggnewscale", quietly = TRUE)) {
    g <- g + ggnewscale::new_scale_colour()
  }

  # 6. 极坐标条形图（外圈）----
  bar <- ggplot(leaf, aes(x = name, y = Count)) +
    geom_col(
      aes(fill = from, alpha = -log10(pvalue)),
      show.legend = FALSE
    ) +
    geom_text(
      aes(
        y     = Count + max_count * 0.06,
        label = name,
        angle = angle,
        hjust = hjust
      ),
      size = bar_label_size
    ) +
    coord_polar() +
    ylim(inner_radius, NA) +
    scale_fill_manual(values = cols) +
    theme_void() +
    theme(
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )

  # 7. 叠合 ----
  result <- bar + inset_element(
    g,
    left     = inset_left,
    bottom   = inset_bottom,
    right    = inset_right,
    top      = inset_top,
    align_to = "full"
  )

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


# 注册
if (exists("register_enrich_plot", mode = "function")) {
  register_enrich_plot(
    type        = "radial",
    fn          = plot_enrich_radial,
    description = "环形树状 + 极坐标条形图（ggraph + patchwork inset）",
    overwrite   = TRUE
  )
}
