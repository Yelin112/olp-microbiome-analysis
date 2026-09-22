#' 柱状图+散点图层构建器（适用于 threshold_small ≤ 最大组样本量 < threshold_medium）
#'
#' @description
#' 返回均值柱状图 + 标准误 errorbar + 全部散点的几何图层列表。
#' 散点使用黑色填充（shape=21）、白色极细边框，覆盖在柱状图上方形成
#' "柱状图骨架 + 个体数据叠加"的视觉层次。
#'
#' @param bar_width 柱宽，默认 0.7
#' @param bar_alpha 柱子透明度，默认 0.85
#' @param point_size 散点大小，默认 2
#' @param point_alpha 散点透明度，默认 0.8
#' @param point_jitter 散点抖动幅度，默认 0.1
#'
#' @return ggplot2 layer 对象列表
#' @keywords internal
build_layers_bar_points <- function(bar_width    = 0.7,
                                     bar_alpha    = 0.85,
                                     point_size   = 2,
                                     point_alpha  = 0.8,
                                     point_jitter = 0.1) {
  list(
    # 均值柱状图
    ggplot2::stat_summary(
      fun       = mean,
      geom      = "bar",
      width     = bar_width,
      alpha     = bar_alpha,
      color     = "black",
      linewidth = 0.5
    ),
    # 标准误 errorbar
    ggplot2::stat_summary(
      fun.data  = mean_se,
      geom      = "errorbar",
      width     = bar_width * 0.35,
      linewidth = 0.7,
      color     = "black"
    ),
    # 散点层：黑色填充 + 白色边框，覆盖在柱状图上方
    ggplot2::geom_jitter(
      shape  = 21,
      size   = point_size,
      stroke = 0.4,
      fill   = "black",
      color  = "white",
      alpha  = point_alpha,
      width  = point_jitter,
      height = 0
    )
  )
}
