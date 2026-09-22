#' 箱线图+散点图层构建器（适用于 threshold_medium ≤ 最大组样本量 < threshold_large）
#'
#' @description
#' 返回箱线图（或小提琴图+内置箱线图）+ 散点的几何图层列表。
#' 箱线图隐藏离群点（由散点层统一展示，避免重复标记）。
#'
#' @param use_violin 是否使用小提琴图替代箱线图，默认 FALSE
#' @param box_width 箱/琴宽，默认 0.6
#' @param box_alpha 箱/琴体透明度，默认 0.7
#' @param point_size 散点大小，默认 1.5
#' @param point_alpha 散点透明度，默认 0.6
#' @param point_jitter 散点抖动幅度，默认 0.15
#'
#' @return ggplot2 layer 对象列表
#' @keywords internal
build_layers_boxplot_points <- function(use_violin   = FALSE,
                                         box_width    = 0.6,
                                         box_alpha    = 0.7,
                                         point_size   = 1.5,
                                         point_alpha  = 0.6,
                                         point_jitter = 0.15) {
  layers <- list()

  if (use_violin) {
    # 小提琴图 + 内置细箱线图
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
        width         = box_width * 0.25,
        fill          = "white",
        outlier.shape = NA,
        color         = "black",
        linewidth     = 0.4,
        alpha         = 0.9
      )
    ))
  } else {
    # 箱线图：隐藏离群点（散点层统一展示）
    layers <- c(layers, list(
      ggplot2::geom_boxplot(
        width         = box_width,
        alpha         = box_alpha,
        outlier.shape = NA,
        color         = "black",
        linewidth     = 0.5
      )
    ))
  }

  # 散点层：白底黑边，叠加在箱体/琴体上方
  layers <- c(layers, list(
    ggplot2::geom_jitter(
      shape  = 21,
      size   = point_size,
      stroke = 0.3,
      fill   = "white",
      color  = "black",
      alpha  = point_alpha,
      width  = point_jitter,
      height = 0
    )
  ))

  return(layers)
}
