# ============================================================================
# panel_fix.R — 固定 ggplot2 面板尺寸工具
# 版本: 2.0
# ============================================================================
#
# 用法:
#   source("R/utils/panel_fix.R")
#
# 依赖: ggplot2, grid
#
# 解决的问题:
#   ggplot2 的 panel 区域大小由数据和图例共同决定，无法固定。
#   当多个图需要拼接时（patchwork / cowplot），子图 panel 大小不一致导致对齐错误。
#   本工具通过对 gtable 的 panel 尺寸做强制覆盖来解决这个问题。
#
# 两种使用方式:
#
#   1. 交互式 (直接在 RStudio 看图):
#      p %>% render_fixed(width = 4, height = 4)
#
#   2. 保存文件 (配合 ggsave):
#      g <- fix_panel_size(p, width = 10, height = 7)
#      dw <- convertWidth(sum(g$widths), "cm", valueOnly = TRUE)
#      dh <- convertHeight(sum(g$heights), "cm", valueOnly = TRUE)
#      ggsave("output.pdf", g, width = dw, height = dh)

library(ggplot2)
library(grid)

#' 固定面板尺寸（仅修改 gtable，不绘制）
#'
#' @description
#' 将 ggplot 对象转换为 gtable，强制设置 panel 区域的物理尺寸（cm），
#' 返回修改后的 gtable。不绘制到设备——适合配合 ggsave 使用。
#'
#' @param plot ggplot 对象。
#' @param width Panel 宽度（cm），默认 4。
#' @param height Panel 高度（cm），默认 4。
#'
#' @return 修改后的 gtable 对象。
#' @export
#'
#' @examples
#' p <- ggplot(iris, aes(Species, Sepal.Length, fill = Species)) + geom_boxplot()
#'
#' # 固定 panel 后保存
#' g <- fix_panel_size(p, width = 8, height = 6)
#' ggsave("output.pdf", g,
#'   width  = convertWidth(sum(g$widths), "cm", valueOnly = TRUE),
#'   height = convertHeight(sum(g$heights), "cm", valueOnly = TRUE)
#' )
fix_panel_size <- function(plot, width = 4, height = 4) {
  w <- unit(width, "cm")
  h <- unit(height, "cm")

  g <- ggplotGrob(plot)

  panels <- grep("panel", g$layout$name)
  if (length(panels) == 0) {
    warning("未找到 panel，无法固定尺寸。返回原始 gtable。")
    return(g)
  }

  panel_cols <- unique(g$layout$l[panels])
  panel_rows <- unique(g$layout$t[panels])

  # 覆盖 panel 的宽度和高度
  if (length(panel_cols) == 1) {
    g$widths[panel_cols] <- w
  } else {
    g$widths[panel_cols] <- rep(w, length(panel_cols))
  }
  if (length(panel_rows) == 1) {
    g$heights[panel_rows] <- h
  } else {
    g$heights[panel_rows] <- rep(h, length(panel_rows))
  }

  return(g)
}


#' 固定面板尺寸并渲染到当前设备
#'
#' @description
#' 调用 fix_panel_size() 修改 gtable，然后绘制到当前设备。
#' 适合 RStudio 中交互式查看。
#'
#' @inheritParams fix_panel_size
#' @return 不可见地返回 gtable 对象。
#' @export
render_fixed <- function(plot, width = 4, height = 4) {
  g <- fix_panel_size(plot, width, height)
  grid.newpage()
  grid.draw(g)
  invisible(g)
}
