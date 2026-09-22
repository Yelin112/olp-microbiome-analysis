# ============================================================================
# palette_system.R — 统一颜色管理系统
# 版本: 1.0
# ============================================================================
#
# 用法:
#   source("R/utils/helpers.R")
#   source("R/utils/palette_system.R")
#
# 依赖: ggplot2, scales, grDevices, helpers.R
#
# 设计原则:
#   - 本文件是整个项目中"颜色的唯一真相来源"
#   - 所有绘图函数通过 scale_fill_pub_d() 或 get_colors() 获取颜色
#   - 不在此文件之外的任何地方定义调色板列表
#   - 运行时注册自定义配色（set_palette），支持跨脚本共享

library(ggplot2)

# ============================================================================
# 1. 内置调色板库（可随时扩展）
# ============================================================================

.builtin_palettes <- list(
  # ── 综合期刊配色 ──
  "NPG" = c(
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
  "AAAS" = c(
    "#3B4992",
    "#EE0000",
    "#008B45",
    "#631879",
    "#008280",
    "#BB0021",
    "#5F559B",
    "#A20056",
    "#808180"
  ),
  "NEJM" = c(
    "#BC3C29",
    "#0072B5",
    "#E18727",
    "#20854E",
    "#7876B1",
    "#6F99AD",
    "#FFDC91",
    "#EE4C97"
  ),
  "Lancet" = c(
    "#00468B",
    "#ED0000",
    "#42B540",
    "#0099B4",
    "#925E9F",
    "#FDAF91",
    "#AD002A",
    "#ADB6B6"
  ),
  "JAMA" = c(
    "#374E55",
    "#DF8F44",
    "#00A1D5",
    "#B24745",
    "#79AF97",
    "#6A6599",
    "#80796B"
  ),
  "JCO" = c(
    "#0073C2",
    "#EFC000",
    "#868686",
    "#CD534C",
    "#7AA6DC",
    "#003C67",
    "#8F7700",
    "#3B3B3B"
  ),
  # ── 数据可视化通用 ──
  "D3" = c(
    "#1F77B4",
    "#FF7F0E",
    "#2CA02C",
    "#D62728",
    "#9467BD",
    "#8C564B",
    "#E377C2",
    "#7F7F7F",
    "#BCBD22",
    "#17BECF"
  ),
  # ── 色盲友好 ──
  "ColorBlind" = c(
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
  ),
  # ── 连续型基准 ──
  "RWB" = c("#313695", "#FFFFFF", "#A50026"),
  "Viridis" = c("#440154", "#21908C", "#FDE725")
)

# ============================================================================
# 2. 单例调色板管理器
# ============================================================================

.palette_registry <- new.env(parent = emptyenv())
.palette_registry$custom <- list() # 用户注册的色板
.palette_registry$default <- "NPG" # 全局默认

# ============================================================================
# 3. 颜色向量校验
# ============================================================================

#' 验证颜色格式（内部）
#' @keywords internal
.is_valid_color <- function(colors) {
  is_hex <- grepl("^#([0-9A-Fa-f]{3}){1,2}([0-9A-Fa-f]{2})?$", colors)
  is_named <- colors %in% grDevices::colors()
  is_hex | is_named
}

# ============================================================================
# 4. 调色板 CRUD
# ============================================================================

#' 注册自定义调色板
#'
#' @param name_or_list 字符串（单个色板名）或命名 list（批量注册）。
#' @param colors 颜色向量。批量模式时此参数省略。
#' @param overwrite 是否允许覆盖已存在的同名色板。
#'
#' @examples
#' # 单个注册
#' set_palette("my_colors", c("#FF0000", "#00FF00", "#0000FF"))
#'
#' # 批量注册
#' set_palette(list(
#'   "phyla"    = c("#E64B35", "#4DBBD5", "#00A087"),
#'   "subtypes" = c("#BC3C29", "#0072B5", "#E18727")
#' ))
#' @export
set_palette <- function(name_or_list, colors = NULL, overwrite = FALSE) {
  # ── 批量模式 ──
  if (is.list(name_or_list)) {
    pal_list <- name_or_list
    if (is.null(names(pal_list)) || any(names(pal_list) == "")) {
      stop("批量导入时，list 必须为每个元素命名")
    }
    for (nm in names(pal_list)) {
      set_palette(nm, pal_list[[nm]], overwrite = overwrite)
    }
    message("✓ 批量注册了 ", length(pal_list), " 个色板")
    return(invisible(pal_list))
  }

  # ── 单个模式 ──
  name <- name_or_list
  if (!is.character(name) || nchar(name) == 0) {
    stop("色板名必须是非空字符串")
  }
  if (!is.character(colors) || length(colors) == 0) {
    stop("颜色必须是非空字符串向量")
  }

  invalid <- !.is_valid_color(colors)
  if (any(invalid)) {
    warning("可能存在无效颜色: ", paste(colors[invalid], collapse = ", "))
  }

  if (name %in% names(.palette_registry$custom) && !overwrite) {
    stop("色板 '", name, "' 已存在。使用 overwrite = TRUE 覆盖。")
  }

  .palette_registry$custom[[name]] <- colors
  message("✓ 色板 '", name, "' 已注册 (", length(colors), " 色)")
  invisible(colors)
}

#' 列出可用调色板
#'
#' @param type "all"（全部）、"builtin"（内置）、"custom"（自定义）。
#' @export
list_palettes <- function(type = c("all", "builtin", "custom")) {
  type <- match.arg(type)
  switch(
    type,
    all = c(names(.builtin_palettes), names(.palette_registry$custom)),
    builtin = names(.builtin_palettes),
    custom = names(.palette_registry$custom)
  )
}

#' 删除自定义调色板
#' @export
remove_palette <- function(name) {
  if (!name %in% names(.palette_registry$custom)) {
    warning("色板 '", name, "' 不存在")
    return(invisible(FALSE))
  }
  .palette_registry$custom[[name]] <- NULL
  message("✓ 色板 '", name, "' 已移除")
  invisible(TRUE)
}

#' 设置全局默认调色板
#' @export
set_default_palette <- function(palette) {
  if (!palette %in% list_palettes("all")) {
    stop("色板 '", palette, "' 不存在")
  }
  .palette_registry$default <- palette
  message("✓ 默认色板设为 '", palette, "'")
}

# ============================================================================
# 5. 核心颜色解析引擎
# ============================================================================

#' 获取颜色向量
#'
#' @description
#' 四级优先级解析: 自定义注册 → 内置库 → RColorBrewer → paletteer → 回退默认。
#' 这是整个颜色系统唯一的解析入口。
#'
#' @param palette 色板名（字符串）、颜色向量、或 NULL（使用全局默认）。
#' @param n 需要的颜色数量。NULL = 返回全部。
#' @param type "discrete"（离散，默认）或 "continuous"（连续）。
#'
#' @return 颜色字符向量
#' @export
get_colors <- function(
  palette = NULL,
  n = NULL,
  type = c("discrete", "continuous")
) {
  type <- match.arg(type)

  # 使用全局默认
  if (is.null(palette)) {
    palette <- .palette_registry$default
  }

  # ── 情况 1: 用户直接传入颜色向量 ──
  if (length(palette) > 1) {
    cols <- as.character(palette)

    # ── 情况 2: 字符串名，按优先级解析 ──
  } else {
    # 优先级 1: 自定义注册
    cols <- .palette_registry$custom[[palette]]

    # 优先级 2: 内置库
    if (is.null(cols)) {
      cols <- .builtin_palettes[[palette]]
    }

    # 优先级 3: RColorBrewer
    if (
      is.null(cols) &&
        requireNamespace("RColorBrewer", quietly = TRUE) &&
        palette %in% rownames(RColorBrewer::brewer.pal.info)
    ) {
      max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
      cols <- RColorBrewer::brewer.pal(
        max(3, min(n %||% 3, max_n)),
        palette
      )
    }

    # 优先级 4: paletteer（离散 + 连续双重回退）
    if (
      is.null(cols) &&
        requireNamespace("paletteer", quietly = TRUE)
    ) {
      cols <- tryCatch(
        # 尝试离散型
        as.character(paletteer::paletteer_d(palette)),
        error = function(e1) {
          tryCatch(
            # 尝试连续型（提取 100 色作为基准池）
            as.character(paletteer::paletteer_c(palette, n = 100)),
            error = function(e2) NULL
          )
        }
      )
    }

    # 终极回退: 全局默认
    if (is.null(cols)) {
      warning(
        "色板 '",
        palette,
        "' 未找到，回退至默认 '",
        .palette_registry$default,
        "'"
      )
      cols <- .builtin_palettes[[.palette_registry$default]]
      if (is.null(cols)) {
        cols <- .palette_registry$custom[[.palette_registry$default]]
      }
    }
  }

  # ── 命名向量严格映射: 不截取不插值 ──
  if (!is.null(names(cols)) && type == "discrete") {
    return(cols)
  }

  # ── 清除 names ──
  cols <- unname(cols)

  # ── 不需要特定数量，直接返回 ──
  if (is.null(n)) {
    return(cols)
  }

  # ── 颜色不足时插值 ──
  if (n > length(cols) || type == "continuous") {
    return(grDevices::colorRampPalette(cols)(n))
  }

  return(cols[seq_len(n)])
}

# ============================================================================
# 6. 可视化预览
# ============================================================================

#' 预览调色板
#'
#' @param palette 色板名或颜色向量。NULL = 全局默认。
#' @param n 显示的颜色数
#' @export
preview_palette <- function(palette = NULL, n = NULL) {
  name <- if (is.null(palette)) {
    .palette_registry$default
  } else if (length(palette) == 1 && is.character(palette)) {
    palette
  } else {
    "Custom"
  }

  cols <- get_colors(palette, n = n, type = "discrete")
  labels <- names(cols) %||% seq_along(cols)

  df <- data.frame(
    x = seq_along(cols),
    y = 1,
    label = labels,
    color = cols
  )

  ggplot(df, aes(x = x, y = y, fill = color)) +
    geom_tile(color = "white", linewidth = 2) +
    geom_text(aes(label = label), size = 3.5, fontface = "bold") +
    scale_fill_identity() +
    labs(
      title = paste0("Palette: ", name, " (", length(cols), " colors)"),
      x = NULL,
      y = NULL
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
}

# ============================================================================
# 7. ggplot2 Scale 函数（上层绘图函数的对接点）
# ============================================================================

#' 离散型填充色标
#'
#' @param palette 色板名、颜色向量、或 NULL
#' @param alpha 透明度 (0–1)
#' @param ... 传递给 scale_fill_manual 或 discrete_scale
#' @export
scale_fill_pub_d <- function(palette = NULL, alpha = 1, ...) {
  cols <- get_colors(palette, type = "discrete")

  if (!is.null(names(cols))) {
    # 命名向量 → 严格映射
    if (alpha < 1) {
      cols <- setNames(scales::alpha(cols, alpha), names(cols))
    }
    return(scale_fill_manual(values = cols, ...))
  }

  # 无命名 → 灵活调色
  pal_func <- function(n) {
    scales::alpha(get_colors(palette, n = n, type = "discrete"), alpha)
  }
  discrete_scale("fill", "pub_d", palette = pal_func, ...)
}

#' 离散型颜色标（线条/点）
#' @rdname scale_fill_pub_d
#' @export
scale_color_pub_d <- function(palette = NULL, alpha = 1, ...) {
  cols <- get_colors(palette, type = "discrete")

  if (!is.null(names(cols))) {
    if (alpha < 1) {
      cols <- setNames(scales::alpha(cols, alpha), names(cols))
    }
    return(scale_color_manual(values = cols, ...))
  }

  pal_func <- function(n) {
    scales::alpha(get_colors(palette, n = n, type = "discrete"), alpha)
  }
  discrete_scale("colour", "pub_d", palette = pal_func, ...)
}

#' 连续型填充色标
#'
#' @param palette 色板名（默认 "RWB"）
#' @param alpha 透明度
#' @param reverse 是否反转色序
#' @param ... 传递给 scale_fill_gradientn
#' @export
scale_fill_pub_c <- function(palette = NULL, alpha = 1, reverse = FALSE, ...) {
  palette <- palette %||% "RWB"
  cols <- get_colors(palette, type = "continuous")
  if (reverse) {
    cols <- rev(cols)
  }
  if (alpha < 1) {
    cols <- scales::alpha(cols, alpha)
  }
  scale_fill_gradientn(colours = cols, ...)
}

#' 连续型颜色标
#' @rdname scale_fill_pub_c
#' @export
scale_color_pub_c <- function(palette = NULL, alpha = 1, reverse = FALSE, ...) {
  palette <- palette %||% "RWB"
  cols <- get_colors(palette, type = "continuous")
  if (reverse) {
    cols <- rev(cols)
  }
  if (alpha < 1) {
    cols <- scales::alpha(cols, alpha)
  }
  scale_color_gradientn(colours = cols, ...)
}
