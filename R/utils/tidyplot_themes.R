# =============================================================================
# Tidyplots 主题管理系统
# 作者: Assistant
# 版本: 1.0
# 说明: 提供多种预设科研主题，支持灵活自定义
# =============================================================================

library(tidyplots)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1. 配色方案库
# -----------------------------------------------------------------------------

COLOR_PALETTES <- list(
  # Nature 系列配色
  nature = c(
    "#E64B35",
    "#4DBBD5",
    "#00A087",
    "#3C5488",
    "#F39B7F",
    "#8491B4",
    "#91D1C2",
    "#DC0000",
    "#7E6148"
  ),

  # JAMA 医学期刊配色
  jama = c(
    "#374E55",
    "#DF8F44",
    "#00A1D5",
    "#B24745",
    "#79AF97",
    "#6A6599",
    "#80796B"
  ),

  # Lancet 配色
  lancet = c(
    "#00468B",
    "#ED0000",
    "#42B540",
    "#0099B4",
    "#925E9F",
    "#FDAF91",
    "#AD002A"
  ),

  # NEJM 配色
  nejm = c(
    "#BC3C29",
    "#0072B5",
    "#E18727",
    "#20854E",
    "#7876B1",
    "#6F99AD",
    "#FFDC91",
    "#EE4C97"
  ),

  # Science 配色
  science = c(
    "#3B4992",
    "#EE0000",
    "#008B45",
    "#631879",
    "#008280",
    "#BB0021",
    "#5F559B"
  ),

  # NPG (Nature Publishing Group)
  npg = c(
    "#E64B35",
    "#4DBBD5",
    "#00A087",
    "#3C5488",
    "#F39B7F",
    "#8491B4",
    "#91D1C2",
    "#DC0000"
  ),

  # Cell Press 配色
  cell = c(
    "#E41A1C",
    "#377EB8",
    "#4DAF4A",
    "#984EA3",
    "#FF7F00",
    "#FFFF33",
    "#A65628",
    "#F781BF"
  ),

  # 简约灰度
  grayscale = c(
    "#000000",
    "#404040",
    "#808080",
    "#BFBFBF",
    "#E0E0E0",
    "#F5F5F5"
  ),

  # 彩虹色（数据多时使用）
  rainbow = c(
    "#FF6B6B",
    "#4ECDC4",
    "#45B7D1",
    "#FFA07A",
    "#98D8C8",
    "#F7DC6F",
    "#BB8FCE",
    "#85C1E2"
  ),

  # 色盲友好配色
  colorblind = c(
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
  ),

  # 暗色主题配色
  dark = c(
    "#FF6B9D",
    "#C44569",
    "#F8B500",
    "#54A0FF",
    "#00D2D3",
    "#5F27CD",
    "#01A3A4"
  )
)

# -----------------------------------------------------------------------------
# 2. 主题预设配置
# -----------------------------------------------------------------------------

THEME_PRESETS <- list(
  # 【Nature 期刊风格】- 最常用
  nature = list(
    name = "Nature Style",
    palette = COLOR_PALETTES$nature,
    font_family = "Arial",
    font_face = "plain",
    font_size = 8,
    border_width = 0.5,
    border_color = "black",
    grid_major = "grey92",
    grid_minor = FALSE,
    axis_line = FALSE, # 使用边框代替轴线
    legend_position = "right",
    background = "white"
  ),

  # 【JAMA 医学风格】- 医学研究
  jama = list(
    name = "JAMA Medical Style",
    palette = COLOR_PALETTES$jama,
    font_family = "Times",
    font_face = "plain",
    font_size = 9,
    border_width = 0.8,
    border_color = "black",
    grid_major = "grey85",
    grid_minor = FALSE,
    axis_line = FALSE,
    legend_position = "top",
    background = "white"
  ),

  # 【Science 期刊风格】- 简洁现代
  science = list(
    name = "Science Magazine Style",
    palette = COLOR_PALETTES$science,
    font_family = "Helvetica",
    font_face = "bold",
    font_size = 10,
    border_width = 1,
    border_color = "black",
    grid_major = FALSE, # Science 通常不要网格线
    grid_minor = FALSE,
    axis_line = TRUE,
    legend_position = "right",
    background = "white"
  ),

  # 【简约现代】- 演示/报告
  minimal = list(
    name = "Minimal Modern",
    palette = COLOR_PALETTES$colorblind,
    font_family = "sans",
    font_face = "plain",
    font_size = 11,
    border_width = 0,
    border_color = NA,
    grid_major = "grey90",
    grid_minor = FALSE,
    axis_line = TRUE,
    legend_position = "bottom",
    background = "white"
  ),

  # 【经典黑白】- 打印友好
  classic = list(
    name = "Classic Black & White",
    palette = COLOR_PALETTES$grayscale,
    font_family = "serif",
    font_face = "plain",
    font_size = 10,
    border_width = 1,
    border_color = "black",
    grid_major = "grey80",
    grid_minor = "grey95",
    axis_line = FALSE,
    legend_position = "right",
    background = "white"
  ),

  # 【暗色主题】- 演示/屏幕展示
  dark = list(
    name = "Dark Theme",
    palette = COLOR_PALETTES$dark,
    font_family = "sans",
    font_face = "bold",
    font_size = 12,
    border_width = 0,
    border_color = NA,
    grid_major = "grey30",
    grid_minor = FALSE,
    axis_line = FALSE,
    legend_position = "right",
    background = "grey20",
    text_color = "white"
  ),

  # 【Cell Press 风格】- 生物医学
  cell = list(
    name = "Cell Press Style",
    palette = COLOR_PALETTES$cell,
    font_family = "Arial",
    font_face = "plain",
    font_size = 8,
    border_width = 0.75,
    border_color = "black",
    grid_major = "grey90",
    grid_minor = FALSE,
    axis_line = FALSE,
    legend_position = "right",
    background = "white"
  ),

  # 【演示专用】- PPT/Keynote
  presentation = list(
    name = "Presentation Style",
    palette = COLOR_PALETTES$rainbow,
    font_family = "sans",
    font_face = "bold",
    font_size = 14,
    border_width = 1.5,
    border_color = "black",
    grid_major = "grey85",
    grid_minor = FALSE,
    axis_line = FALSE,
    legend_position = "bottom",
    background = "white"
  )
)

# -----------------------------------------------------------------------------
# 3. 核心主题应用函数
# -----------------------------------------------------------------------------

#' 应用 Tidyplots 主题
#'
#' @param plot tidyplot 对象
#' @param theme 主题名称或自定义配置列表
#' @param palette 可选，覆盖主题的配色方案
#' @param font_size 可选，覆盖字体大小
#' @param ... 其他自定义参数
#'
#' @return 应用主题后的 tidyplot 对象
#'
#' @examples
#' # 使用预设主题
#' plot %>% apply_tidyplot_theme("nature")
#'
#' # 自定义配色
#' plot %>% apply_tidyplot_theme("nature", palette = c("#FF0000", "#00FF00"))
#'
#' # 调整字体大小
#' plot %>% apply_tidyplot_theme("jama", font_size = 12)

apply_tidyplot_theme <- function(
  plot,
  theme = "nature",
  palette = NULL,
  font_size = NULL,
  legend_position = NULL,
  border_width = NULL,
  ...
) {
  # 获取主题配置
  if (is.character(theme)) {
    if (!theme %in% names(THEME_PRESETS)) {
      stop(paste(
        "Unknown theme:",
        theme,
        "\nAvailable themes:",
        paste(names(THEME_PRESETS), collapse = ", ")
      ))
    }
    config <- THEME_PRESETS[[theme]]
  } else if (is.list(theme)) {
    config <- theme
  } else {
    stop("theme must be a character string or list")
  }

  # 参数覆盖
  if (!is.null(palette)) {
    config$palette <- palette
  }
  if (!is.null(font_size)) {
    config$font_size <- font_size
  }
  if (!is.null(legend_position)) {
    config$legend_position <- legend_position
  }
  if (!is.null(border_width)) {
    config$border_width <- border_width
  }

  # 处理额外参数
  extra_params <- list(...)
  for (param_name in names(extra_params)) {
    config[[param_name]] <- extra_params[[param_name]]
  }

  # 应用配色
  color_scheme <- new_color_scheme(config$palette, name = config$name)
  plot <- plot %>% adjust_colors(color_scheme)

  # 应用字体
  plot <- plot %>%
    adjust_font(
      family = config$font_family,
      face = config$font_face
    )

  # 构建主题细节参数
  theme_params <- list()

  # 文字颜色（暗色主题特殊处理）
  text_color <- if (!is.null(config$text_color)) config$text_color else "black"

  # 面板设置
  theme_params$panel.background <- element_rect(
    fill = config$background,
    colour = NA
  )

  # 边框
  if (!is.na(config$border_color) && config$border_width > 0) {
    theme_params$panel.border <- element_rect(
      colour = config$border_color,
      fill = NA,
      linewidth = config$border_width
    )
  } else {
    theme_params$panel.border <- element_blank()
  }

  # 网格线
  if (is.character(config$grid_major)) {
    theme_params$panel.grid.major <- element_line(
      colour = config$grid_major,
      linewidth = 0.3
    )
  } else if (config$grid_major == FALSE) {
    theme_params$panel.grid.major <- element_blank()
  }

  if (is.character(config$grid_minor)) {
    theme_params$panel.grid.minor <- element_line(
      colour = config$grid_minor,
      linewidth = 0.2
    )
  } else {
    theme_params$panel.grid.minor <- element_blank()
  }

  # 坐标轴线
  if (config$axis_line) {
    theme_params$axis.line <- element_line(
      colour = text_color,
      linewidth = 0.5
    )
  } else {
    theme_params$axis.line <- element_blank()
  }

  # 坐标轴刻度
  theme_params$axis.ticks <- element_line(
    colour = text_color,
    linewidth = 0.5
  )
  theme_params$axis.ticks.length <- unit(0.15, "cm")

  # 坐标轴文字
  theme_params$axis.text.x <- element_text(
    colour = text_color,
    size = config$font_size,
    margin = margin(t = 5)
  )
  theme_params$axis.text.y <- element_text(
    colour = text_color,
    size = config$font_size,
    margin = margin(r = 5)
  )
  theme_params$axis.title <- element_text(
    colour = text_color,
    size = config$font_size + 1,
    face = "bold"
  )

  # 图例
  theme_params$legend.position <- config$legend_position
  theme_params$legend.background <- element_blank()
  theme_params$legend.key <- element_blank()
  theme_params$legend.title <- element_text(
    colour = text_color,
    size = config$font_size,
    face = "bold"
  )
  theme_params$legend.text <- element_text(
    colour = text_color,
    size = config$font_size - 1
  )

  # 标题
  theme_params$plot.title <- element_text(
    colour = text_color,
    size = config$font_size + 2,
    face = "bold",
    hjust = 0.5,
    margin = margin(b = 10)
  )

  # 背景
  theme_params$plot.background <- element_rect(
    fill = config$background,
    colour = NA
  )

  # 应用主题
  plot <- do.call(adjust_theme_details, c(list(plot), theme_params))

  return(plot)
}

# -----------------------------------------------------------------------------
# 4. 便捷别名函数
# -----------------------------------------------------------------------------

# 为常用主题创建快捷函数
theme_nature <- function(plot, ...) apply_tidyplot_theme(plot, "nature", ...)
theme_jama <- function(plot, ...) apply_tidyplot_theme(plot, "jama", ...)
theme_science <- function(plot, ...) apply_tidyplot_theme(plot, "science", ...)
theme_minimal <- function(plot, ...) apply_tidyplot_theme(plot, "minimal", ...)
theme_classic <- function(plot, ...) apply_tidyplot_theme(plot, "classic", ...)
theme_dark <- function(plot, ...) apply_tidyplot_theme(plot, "dark", ...)
theme_cell <- function(plot, ...) apply_tidyplot_theme(plot, "cell", ...)
theme_presentation <- function(plot, ...) {
  apply_tidyplot_theme(plot, "presentation", ...)
}

# -----------------------------------------------------------------------------
# 5. 主题预览函数
# -----------------------------------------------------------------------------

#' 预览所有可用主题
#'
#' @param data 示例数据框
#' @param save_path 可选，保存路径
#'
#' @examples
#' preview_all_themes(mtcars)

preview_all_themes <- function(data = NULL, save_path = NULL) {
  # 使用默认数据
  if (is.null(data)) {
    data <- data.frame(
      Group = rep(c("Control", "Treatment A", "Treatment B"), each = 10),
      Value = c(rnorm(10, 5, 1), rnorm(10, 7, 1.2), rnorm(10, 6.5, 0.8))
    )
  }

  cat("📊 Available Tidyplot Themes:\n")
  cat(rep("=", 60), "\n\n", sep = "")

  for (theme_name in names(THEME_PRESETS)) {
    config <- THEME_PRESETS[[theme_name]]
    cat(sprintf("✓ %s (%s)\n", theme_name, config$name))
    cat(sprintf(
      "  Palette: %d colors | Font: %s | Size: %dpt\n",
      length(config$palette),
      config$font_family,
      config$font_size
    ))

    # 显示配色
    cat("  Colors: ")
    cat(paste(
      config$palette[1:min(5, length(config$palette))],
      collapse = ", "
    ))
    if (length(config$palette) > 5) {
      cat(", ...")
    }
    cat("\n\n")
  }

  cat(rep("=", 60), "\n", sep = "")
  cat("\nUsage examples:\n")
  cat("  plot %>% apply_tidyplot_theme('nature')\n")
  cat("  plot %>% theme_jama()\n")
  cat("  plot %>% theme_science(font_size = 12)\n")
}

# -----------------------------------------------------------------------------
# 6. 自定义主题创建函数
# -----------------------------------------------------------------------------

#' 创建自定义主题
#'
#' @param name 主题名称
#' @param palette 配色方案向量
#' @param ... 其他主题参数
#'
#' @return 主题配置列表
#'
#' @examples
#' my_theme <- create_custom_theme(
#'   name = "My Custom Theme",
#'   palette = c("#FF5733", "#33FF57", "#3357FF"),
#'   font_family = "Arial",
#'   font_size = 10
#' )
#' plot %>% apply_tidyplot_theme(my_theme)

create_custom_theme <- function(
  name = "Custom Theme",
  palette = COLOR_PALETTES$nature,
  font_family = "Arial",
  font_face = "plain",
  font_size = 10,
  border_width = 0.5,
  border_color = "black",
  grid_major = "grey92",
  grid_minor = FALSE,
  axis_line = FALSE,
  legend_position = "right",
  background = "white",
  text_color = "black",
  ...
) {
  config <- list(
    name = name,
    palette = palette,
    font_family = font_family,
    font_face = font_face,
    font_size = font_size,
    border_width = border_width,
    border_color = border_color,
    grid_major = grid_major,
    grid_minor = grid_minor,
    axis_line = axis_line,
    legend_position = legend_position,
    background = background,
    text_color = text_color
  )

  # 添加额外参数
  extra_params <- list(...)
  config <- c(config, extra_params)

  class(config) <- c("tidyplot_theme", "list")
  return(config)
}

# -----------------------------------------------------------------------------
# 7. 配色方案管理
# -----------------------------------------------------------------------------

#' 列出所有可用配色
list_color_palettes <- function() {
  cat("🎨 Available Color Palettes:\n")
  cat(rep("=", 60), "\n\n", sep = "")

  for (pal_name in names(COLOR_PALETTES)) {
    colors <- COLOR_PALETTES[[pal_name]]
    cat(sprintf("%-15s: %d colors\n", pal_name, length(colors)))
    cat(sprintf(
      "                %s\n\n",
      paste(colors[1:min(5, length(colors))], collapse = ", ")
    ))
  }
}

#' 添加自定义配色方案
#'
#' @param name 配色方案名称
#' @param colors 颜色向量

add_color_palette <- function(name, colors) {
  COLOR_PALETTES[[name]] <<- colors
  message(sprintf(
    "✓ Added color palette '%s' with %d colors",
    name,
    length(colors)
  ))
}

# -----------------------------------------------------------------------------
# 8. 使用说明
# -----------------------------------------------------------------------------

print_usage_guide <- function() {
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════╗\n")
  cat("║         Tidyplots 主题管理系统 - 使用指南                 ║\n")
  cat("╚════════════════════════════════════════════════════════════╝\n\n")

  cat("【基础用法】\n")
  cat("  # 使用预设主题\n")
  cat("  plot %>% apply_tidyplot_theme('nature')\n")
  cat("  plot %>% theme_jama()  # 快捷函数\n\n")

  cat("【自定义参数】\n")
  cat("  # 修改字体大小\n")
  cat("  plot %>% theme_nature(font_size = 12)\n\n")
  cat("  # 修改配色\n")
  cat("  plot %>% theme_science(palette = c('#FF0000', '#00FF00'))\n\n")
  cat("  # 修改图例位置\n")
  cat("  plot %>% theme_cell(legend_position = 'bottom')\n\n")

  cat("【创建自定义主题】\n")
  cat("  my_theme <- create_custom_theme(\n")
  cat("    name = 'My Theme',\n")
  cat("    palette = c('#E64B35', '#4DBBD5'),\n")
  cat("    font_size = 11,\n")
  cat("    border_width = 1\n")
  cat("  )\n")
  cat("  plot %>% apply_tidyplot_theme(my_theme)\n\n")

  cat("【查看可用资源】\n")
  cat("  preview_all_themes()      # 查看所有主题\n")
  cat("  list_color_palettes()     # 查看所有配色\n\n")

  cat("【可用主题】\n")
  cat("  ", paste(names(THEME_PRESETS), collapse = ", "), "\n\n")
}

# 加载时显示使用指南
print_usage_guide()
