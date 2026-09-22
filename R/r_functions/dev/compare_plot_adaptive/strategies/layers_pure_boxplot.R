#' 纯箱线图层构建器（适用于最大组样本量 ≥ threshold_large）
#'
#' @description
#' 返回箱线图（或小提琴图+内置箱线图）+ 均值点标记的几何图层列表。
#' 大样本量时默认不叠加散点以避免过度绘制（可通过 show_points = TRUE 开启）。
#'
#' @param use_violin 是否使用小提琴图，默认 FALSE
#' @param show_points 是否叠加散点，默认 FALSE
#' @param box_width 箱/琴宽，默认 0.7
#' @param box_alpha 箱/琴体透明度，默认 0.8
#' @param point_size 散点大小（show_points = TRUE 时生效），默认 1
#' @param point_alpha 散点透明度（show_points = TRUE 时生效），默认 0.4
#'
#' @return ggplot2 layer 对象列表
#' @keywords internal
build_layers_pure_boxplot <- function(use_violin  = FALSE,
                                       show_points = FALSE,
                                       box_width   = 0.7,
                                       box_alpha   = 0.8,
                                       point_size  = 1,
                                       point_alpha = 0.4) {
  layers <- list()

  if (use_violin) {
    # 小提琴图 + 内置极细箱线图
    layers <- c(layers, list(
      ggplot2::geom_violin(
        scale     = "width",
        trim      = FALSE,
        width     = box_width,
        alpha     = box_alpha,
        color     = "black",
        linewidth = 0.5
      ),
      ggplot2::geom_boxplot(
        width         = box_width * 0.2,
        fill          = "white",
        outlier.shape = NA,
        color         = "black",
        linewidth     = 0.4,
        alpha         = 0.9
      )
    ))
  } else {
    # 箱线图：根据 show_points 决定是否显示离群点
    layers <- c(layers, list(
      ggplot2::geom_boxplot(
        width         = box_width,
        alpha         = box_alpha,
        outlier.shape = if (show_points) NA else 19,
        outlier.size  = 1.5,
        outlier.alpha = 0.5,
        color         = "black",
        linewidth     = 0.5
      )
    ))
  }

  # 均值点标记（白底黑边菱形）
  layers <- c(layers, list(
    ggplot2::stat_summary(
      fun         = mean,
      geom        = "point",
      shape       = 21,
      size        = 2.5,
      fill        = "white",
      color       = "black",
      stroke      = 0.8,
      show.legend = FALSE
    )
  ))

  # 可选散点（大数据量时建议关闭，使用实心圆 shape=16 减少视觉重量）
  if (show_points) {
    layers <- c(layers, list(
      ggplot2::geom_jitter(
        shape  = 16,
        size   = point_size,
        alpha  = point_alpha,
        color  = "black",
        width  = 0.2,
        height = 0
      )
    ))
  }

  return(layers)
}
