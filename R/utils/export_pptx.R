# ============================================================================
# export_pptx.R — ggplot2 → 可编辑 PPTX 导出工具
# 版本: 1.0
# ============================================================================
#
# 用法:
#   source("R/utils/helpers.R")
#   source("R/utils/panel_fix.R")
#   source("R/utils/export_pptx.R")
#
# 依赖: officer, rvg, ggplot2, grid, panel_fix.R (fix_panel_size)
#
# 设计原则:
#   - 用 officer + rvg::dml() 把 ggplot 渲染为 Office DrawingML 矢量图形，
#     插入 PPTX 后图形元素（线条/填充/文字）在 PowerPoint 中完全可编辑，
#     而不是像插入 PNG 那样是不可编辑的位图。
#   - 复用 panel_fix.R 的 fix_panel_size()，保证单图 PPTX 版本与同一图的
#     PDF/PNG 版本 panel 尺寸/比例一致。patchwork 等拼图对象 panel 结构
#     不适合强制统一尺寸，可设 fix_panel = FALSE 按整体宽高直接导出。

#' 保存单个 ggplot 为可编辑 PPTX（单页幻灯片）
#'
#' @description
#' 用 rvg::dml() 将 ggplot 渲染为矢量 DrawingML，插入单页 PPTX。
#'
#' @param plot ggplot 对象（也可以是 patchwork 拼图对象，此时建议
#'   fix_panel = FALSE）。
#' @param file 输出路径，需以 .pptx 结尾。
#' @param panel_width,panel_height 尺寸（cm）。fix_panel = TRUE 时为传给
#'   fix_panel_size() 的 panel 尺寸；fix_panel = FALSE 时为整图尺寸。默认 8x6。
#' @param title 幻灯片标题（可选）。NULL 则不显示标题占位符。
#' @param fix_panel 是否用 fix_panel_size() 固定 panel 尺寸（与 PDF/PNG 版本
#'   保持一致）。对单个 compare_plot() 结果建议 TRUE；对 patchwork 拼图建议
#'   FALSE。默认 TRUE。
#'
#' @return 不可见地返回 file 路径；缺 officer/rvg 时警告并返回 NULL（不报错）。
#' @export
#'
#' @examples
#' p <- ggplot(iris, aes(Species, Sepal.Length, fill = Species)) + geom_boxplot()
#' save_plot_pptx(p, "output.pptx", panel_width = 8, panel_height = 6, title = "Sepal Length")
save_plot_pptx <- function(
  plot,
  file,
  panel_width = 8,
  panel_height = 6,
  title = NULL,
  fix_panel = TRUE
) {
  if (!.check_pptx_deps()) {
    return(invisible(NULL))
  }

  slide_content <- .prepare_slide_content(
    plot,
    panel_width,
    panel_height,
    fix_panel
  )

  doc <- officer::read_pptx()
  doc <- .add_plot_slide(
    doc,
    slide_content$value,
    slide_content$w_in,
    slide_content$h_in,
    title
  )

  print(doc, target = file)
  invisible(file)
}

#' 保存多个 ggplot 为一个多页 PPTX（每图一页）
#'
#' @description
#' 适合把一组指标图（如各细胞类型占比图）打包成一份可编辑的 PPTX，
#' 方便后续在 PowerPoint 中调整文字、颜色、图例位置等用于汇报/投稿排版。
#'
#' @param plots 命名 list（ggplot 对象），名称将作为对应幻灯片标题；
#'   未命名时不显示标题。
#' @param file 输出路径，需以 .pptx 结尾。
#' @param panel_width,panel_height 尺寸（cm），对所有图统一应用。默认 8x6。
#' @param fix_panel 是否用 fix_panel_size() 固定 panel 尺寸，见 save_plot_pptx()。默认 TRUE。
#'
#' @return 不可见地返回 file 路径。
#' @export
#'
#' @examples
#' plots <- list("Sepal Length" = p1, "Petal Length" = p2)
#' save_plots_pptx(plots, "output/all_plots.pptx", panel_width = 9, panel_height = 7)
save_plots_pptx <- function(
  plots,
  file,
  panel_width = 8,
  panel_height = 6,
  fix_panel = TRUE
) {
  if (!.check_pptx_deps()) {
    return(invisible(NULL))
  }
  stopifnot(is.list(plots), length(plots) > 0)

  titles <- names(plots) %||% rep("", length(plots))

  doc <- officer::read_pptx()
  for (i in seq_along(plots)) {
    slide_content <- .prepare_slide_content(
      plots[[i]],
      panel_width,
      panel_height,
      fix_panel
    )
    ttl <- if (nzchar(titles[i])) titles[i] else NULL
    doc <- .add_plot_slide(
      doc,
      slide_content$value,
      slide_content$w_in,
      slide_content$h_in,
      ttl
    )
  }

  print(doc, target = file)
  invisible(file)
}

# ==============================================================================
# 内部辅助函数
# ==============================================================================

#' 检查 officer/rvg 是否可用（内部）
#'
#' 缺包时只给出警告并返回 FALSE，调用方跳过 PPTX 导出，
#' 这样在未安装这两个包的服务器上，PDF/PNG 导出不会被中断。
#' @keywords internal
.check_pptx_deps <- function() {
  ok <- requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("rvg", quietly = TRUE)
  if (!ok) {
    warning(
      "缺少 officer/rvg 包，已跳过 PPTX 导出（PDF/PNG 不受影响）。",
      "如需导出：install.packages(c('officer', 'rvg'))",
      call. = FALSE
    )
  }
  ok
}

#' 准备单页幻灯片的矢量内容与尺寸（内部）
#' @keywords internal
.prepare_slide_content <- function(plot, width, height, fix_panel) {
  if (fix_panel) {
    g <- fix_panel_size(plot, width = width, height = height)
    list(
      value = rvg::dml(code = grid::grid.draw(g)),
      w_in = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
      h_in = grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)
    )
  } else {
    list(
      value = rvg::dml(ggobj = plot),
      w_in = width / 2.54,
      h_in = height / 2.54
    )
  }
}

#' 添加一页含矢量图的幻灯片（内部）
#' @keywords internal
.add_plot_slide <- function(doc, value, w_in, h_in, title = NULL) {
  doc <- officer::add_slide(
    doc,
    layout = "Title and Content",
    master = "Office Theme"
  )

  top <- 0.5
  if (!is.null(title)) {
    doc <- officer::ph_with(
      doc,
      value = title,
      location = officer::ph_location_type(type = "title")
    )
    top <- 1.3
  }

  doc <- officer::ph_with(
    doc,
    value = value,
    location = officer::ph_location(
      width = w_in,
      height = h_in,
      left = 0.5,
      top = top
    )
  )

  doc
}
