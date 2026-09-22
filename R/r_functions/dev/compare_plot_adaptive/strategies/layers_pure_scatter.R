#' 纯散点图层构建器（适用于最大组样本量 < threshold_small）
#'
#' @description
#' 返回在统一 aes 映射基础上叠加的几何图层列表。
#' 包含：全部散点 + 均值 crossbar + 标准误 errorbar。
#' 适用于极小样本（n < 5），展示每个数据点同时标注均值和误差。
#'
#' @param point_size 散点大小，默认 3
#' @param point_stroke 散点边框粗细，默认 0.5
#' @param errorbar_size 误差线粗细，默认 1.2
#' @param errorbar_width 误差线帽宽，默认 0.3
#'
#' @return ggplot2 layer 对象列表
#' @keywords internal
#'
#' @note
#' 该函数假设 ggplot 已通过 `aes(x = group, y = value, fill = group)` 初始化，
#' 各图层自动继承此映射。不在此函数内设置 scale_fill_manual 或 theme，
#' 这些由管线统一处理。
build_layers_pure_scatter <- function(point_size    = 3,
                                       point_stroke  = 0.5,
                                       errorbar_size = 1.2,
                                       errorbar_width = 0.3) {
  list(
    # 散点层：全部数据点，轻微水平抖动
    ggplot2::geom_jitter(
      shape  = 21,
      size   = point_size,
      stroke = point_stroke,
      color  = "black",
      width  = 0.15,
      height = 0
    ),
    # 均值线（crossbar 风格，fatten=1 即粗线段）
    ggplot2::stat_summary(
      fun       = mean,
      geom      = "crossbar",
      width     = errorbar_width,
      linewidth = errorbar_size,
      color     = "black",
      fatten    = 1
    ),
    # 标准误 errorbar
    ggplot2::stat_summary(
      fun.data  = mean_se,
      geom      = "errorbar",
      width     = errorbar_width * 0.6,
      linewidth = errorbar_size * 0.7,
      color     = "black"
    )
  )
}
