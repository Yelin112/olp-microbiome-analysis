# ==============================================================================
# plot_enrich_bubble(): 富集分析气泡图
# 创建: 2026-04-24
# ==============================================================================

# 依赖
# - core/std_format.R 中的 as_enrich_std(), prep_enrich_data()
# - ggplot2

# ------------------------------------------------------------------------------
# 默认配色与主题（与模板一致，后续统一进 registry）
# ------------------------------------------------------------------------------

library(ggplot2)

.bubble_default_palette <- c('#ffdbdb', '#f29c98')

.enrich_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12),
      strip.placement = "outside",
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 16)
    )
}


# ------------------------------------------------------------------------------
# plot_enrich_bubble()
# ------------------------------------------------------------------------------

#' 富集分析气泡图
#'
#' 以 GeneRatio（或 Count/RichFactor）为 x 轴，通路名称为 y 轴，
#' 气泡大小代表富集基因数量，颜色深浅代表显著性（-log10 p.adjust）。
#'
#' @param data 富集分析结果，接受：
#'   - EnrichStd 对象（由 as_enrich_std() 转换）
#'   - enrichResult 对象（clusterProfiler 直接输出，会自动转换）
#'   - data.frame（需包含标准列或通过 col_map 指定）
#' @param col_map 列名映射，仅当 data 为非标准 data.frame 时使用，
#'   见 as_enrich_std() 的 col_map 参数说明
#' @param top_n 显示前 N 个 term（按 p_adjust 升序），默认 15
#' @param filter_padj p_adjust 阈值过滤，默认 NULL（不过滤）
#' @param terms 手动指定要显示的 term_name 向量，会覆盖 top_n
#' @param x_var x 轴变量，可选 "gene_ratio"（默认）、"rich_factor"、"gene_count"
#' @param palette 渐变色板，接受任意长度的颜色向量，颜色顺序对应从低值（不显著）到高值（显著）。
#'   默认 c('#ffdbdb', '#f29c98')（两端粉色渐变）。
#'   可传入 RColorBrewer 调色板：RColorBrewer::brewer.pal(9, "Reds")，
#'   或截取其中几个色：brewer.pal(11, "Spectral")[c(1,4,8,11)]
#' @param size_range 气泡大小范围，默认 c(2, 8)
#' @param order_by y 轴排序依据，默认 "p_adjust"（升序 = 最显著在顶）
#'   可选："gene_count", "rich_factor", "gene_ratio"
#' @param x_label x 轴标签，NULL 时根据 x_var 自动生成
#' @param title 图标题，默认 NULL
#' @param legend_title_color 颜色图例标题，默认 "-log10(p.adjust)"
#' @param legend_title_size  大小图例标题，默认 "Gene Count"
#' @param base_size 基础字号，影响整体字体大小，默认 11
#' @param ... 传递给 ggplot2::theme() 的额外主题参数
#'
#' @return ggplot 对象
#'
#' @examples
#' \dontrun{
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#' source("core/std_format.R")
#' source("styles/bubble.R")
#'
#' kegg <- enrichKEGG(gene = my_entrez_ids, organism = "hsa",
#'                    pvalueCutoff = 0.05)
#'
#' # 最简用法（直接传 enrichResult）
#' plot_enrich_bubble(kegg)
#'
#' # 先转换为 EnrichStd，再绘图
#' std <- as_enrich_std(kegg)
#' plot_enrich_bubble(std, top_n = 20, x_var = "rich_factor")
#'
#' # 自定义两色渐变
#' plot_enrich_bubble(kegg,
#'   palette    = c("#d6eaf8", "#1a5276"),
#'   size_range = c(3, 12)
#' )
#'
#' # 使用 RColorBrewer 多色板
#' plot_enrich_bubble(kegg,
#'   palette = RColorBrewer::brewer.pal(9, "YlOrRd")
#' )
#'
#' # 截取 Spectral 中的几个颜色
#' plot_enrich_bubble(kegg,
#'   palette = RColorBrewer::brewer.pal(11, "Spectral")[c(2, 6, 10)]
#' )
#' }
#'
#' @export
plot_enrich_bubble <- function(
  data,
  col_map = NULL,
  top_n = 15,
  filter_padj = NULL,
  terms = NULL,
  x_var = c("gene_ratio", "rich_factor", "gene_count"),
  palette = .bubble_default_palette,
  size_range = c(2, 8),
  order_by = "p_adjust",
  x_label = NULL,
  title = NULL,
  legend_title_color = "-log10(p.adjust)",
  legend_title_size = "Gene Count",
  base_size = 11
) {
  x_var <- match.arg(x_var)

  # 1. 数据标准化 ----
  std <- as_enrich_std(data, col_map = col_map)

  # 2. 筛选 / 排序 ----
  std <- prep_enrich_data(
    std,
    top_n = top_n,
    filter_padj = filter_padj,
    terms = terms,
    order_by = order_by,
    decreasing = FALSE
  )

  # 3. 检查 x_var 可用性 ----
  if (!x_var %in% colnames(std) || all(is.na(std[[x_var]]))) {
    fallback <- c("gene_ratio", "rich_factor", "gene_count")
    fallback <- fallback[fallback != x_var]
    ok <- Filter(
      function(v) v %in% colnames(std) && any(!is.na(std[[v]])),
      fallback
    )
    if (length(ok) == 0) {
      stop("x_var '", x_var, "' 的值全为 NA，且无可用替代列")
    }
    warning("x_var '", x_var, "' 全为 NA，自动切换为 '", ok[1], "'")
    x_var <- ok[1]
  }

  # 4. 固定 term_name 的显示顺序（按 x_var 降序：最大值在顶，最小值在底）
  std <- std[order(std[[x_var]], decreasing = TRUE), ]
  std$term_name <- factor(std$term_name, levels = rev(std$term_name))

  # 5. x 轴标签 ----
  if (is.null(x_label)) {
    x_label <- switch(
      x_var,
      gene_ratio = "Gene Ratio",
      rich_factor = "Rich Factor",
      gene_count = "Gene Count"
    )
  }

  # 6. 绘图 ----
  p <- ggplot(
    std,
    aes(
      x = .data[[x_var]],
      y = term_name,
      fill = -log10(p_adjust),
      size = gene_count
    )
  ) +
    geom_point(shape = 21, color = "white", stroke = 0.3) +
    scale_fill_gradientn(
      colours = colorRampPalette(palette)(50),
      name = legend_title_color
    ) +
    scale_size_continuous(
      range = size_range,
      name = legend_title_size
    ) +
    labs(
      x = x_label,
      y = NULL,
      title = title
    ) +
    .enrich_theme() +
    theme(
      axis.text.y  = element_text(size = base_size),
      axis.text.x  = element_text(size = base_size - 1),
      axis.title.x = element_text(size = base_size + 2),
      legend.text  = element_text(size = base_size - 1),
      plot.title   = element_text(size = base_size + 4, face = "bold")
    )

  p
}


# 文件被 source 时自动注册到主系统
if (exists("register_enrich_plot", mode = "function")) {
  register_enrich_plot(
    type        = "bubble",
    fn          = plot_enrich_bubble,
    description = "气泡图（x=GeneRatio/RichFactor, y=通路, 颜色=p.adjust, 大小=Count）",
    overwrite   = TRUE
  )
}
