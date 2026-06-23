library(circlize)
library(RColorBrewer)
library(dplyr)
library(scales)

# ==============================================================================
# 问题修复说明
# ==============================================================================
# 原脚本存在以下未定义的辅助函数（运行时会直接报错）：
#   1. get_palette_colors()       —— 调色板名称转颜色向量
#   2. .extract_row_order()       —— 提取热图行顺序（支持聚类）
#   3. .process_track_data()      —— 处理/转换单条轨道数据
#   4. .draw_track()              —— 根据 type 分发绘图逻辑
#   5. .calculate_feature_positions() —— 计算每个 feature 在环上的 x 坐标（用于连线）
#   6. .add_heatmap_legend()      —— 在图外绘制热图颜色图例
#
# 额外结构缺陷：
#   a. %||% 运算符未定义（仅 rlang 包内置，此处需自行声明）
#   b. .prepare_heatmap_colors() 情况2中，log message 误用了已被覆盖的 heatmap_col
#      （此时 heatmap_col 已是函数，sprintf 会打印 "function"）
#   c. circos.clear() 在 tryCatch finally 中执行，会清除刚绘制完的图形
#      正确做法：仅在 error 分支清除，finally 不应无条件 clear
# ==============================================================================

# ==============================================================================
# 0. 基础运算符 & 工具
# ==============================================================================

#' Null coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b


# ==============================================================================
# 1. 补全：get_palette_colors
#    支持 RColorBrewer、ggplot2 内置名、以及 scales::hue_pal 作为兜底
# ==============================================================================
get_palette_colors <- function(palette_name, n) {
  # 内置映射表：常见调色板名 -> brewer 或 hue
  brewer_palettes <- c(
    "Set1",
    "Set2",
    "Set3",
    "Pastel1",
    "Pastel2",
    "Paired",
    "Dark2",
    "Accent",
    "Blues",
    "Reds",
    "Greens",
    "Purples",
    "Oranges",
    "RdBu",
    "RdYlBu",
    "Spectral",
    "PiYG",
    "PRGn",
    "BrBG",
    "PuOr",
    "RdGy",
    "RdYlGn"
  )

  if (palette_name %in% brewer_palettes) {
    max_n <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
    base_colors <- RColorBrewer::brewer.pal(min(n, max_n), palette_name)
    if (n > max_n) {
      return(colorRampPalette(base_colors)(n))
    }
    return(base_colors[seq_len(n)])
  }

  # 常用期刊/风格调色板（手动内置）
  custom_palettes <- list(
    npg = c(
      "#E64B35",
      "#4DBBD5",
      "#00A087",
      "#3C5488",
      "#F39B7F",
      "#8491B4",
      "#91D1C2",
      "#DC0000",
      "#7E6148",
      "#B09C85"
    ),
    lancet = c(
      "#00468B",
      "#ED0000",
      "#42B540",
      "#0099B4",
      "#925E9F",
      "#FDAF91",
      "#AD002A",
      "#ADB6B6",
      "#1B1919"
    ),
    jco = c(
      "#0073C2",
      "#EFC000",
      "#868686",
      "#CD534C",
      "#7AA6DC",
      "#003C67",
      "#8F7700",
      "#3B3B3B",
      "#A73030",
      "#4A6990"
    ),
    d3 = c(
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
    tableau10 = c(
      "#4E79A7",
      "#F28E2B",
      "#E15759",
      "#76B7B2",
      "#59A14F",
      "#EDC948",
      "#B07AA1",
      "#FF9DA7",
      "#9C755F",
      "#BAB0AC"
    )
  )

  if (palette_name %in% names(custom_palettes)) {
    base_colors <- custom_palettes[[palette_name]]
    if (n <= length(base_colors)) {
      return(base_colors[seq_len(n)])
    }
    return(colorRampPalette(base_colors)(n))
  }

  # 兜底：使用 scales::hue_pal
  warning(sprintf(
    "Unknown palette '%s', falling back to hue_pal.",
    palette_name
  ))
  scales::hue_pal()(n)
}


# ==============================================================================
# 2. 补全：.extract_row_order
#    返回一个 list，key = group名称，value = 该组内 feature 的有序索引
# ==============================================================================
.extract_row_order <- function(data, group_col, clustered = FALSE) {
  groups <- unique(data[[group_col]])
  row_order <- lapply(setNames(groups, groups), function(g) {
    idx <- which(data[[group_col]] == g)
    if (clustered && length(idx) > 1) {
      # 对该组的 numeric 列进行层次聚类排序
      sub_mat <- data[idx, sapply(data, is.numeric), drop = FALSE]
      if (ncol(sub_mat) > 0 && nrow(sub_mat) > 1) {
        hc <- hclust(dist(sub_mat), method = "ward.D2")
        idx <- idx[hc$order]
      }
    }
    idx
  })
  row_order
}


# ==============================================================================
# 3. 补全：.process_track_data
#    对轨道列进行类型特定的预处理（-log10、rescale、zscore 等）
# ==============================================================================
.process_track_data <- function(data, config) {
  col <- config$column
  type <- config$type
  vals <- data[[col]]

  transform <- config$transform %||% "none"

  if (transform == "log10" || type == "pvalue") {
    vals <- suppressWarnings(-log10(as.numeric(vals) + 1e-300))
  } else if (transform == "log2") {
    vals <- suppressWarnings(log2(as.numeric(vals) + 1))
  } else if (transform == "zscore") {
    mu <- mean(vals, na.rm = TRUE)
    sd <- sd(vals, na.rm = TRUE)
    if (sd > 0) vals <- (vals - mu) / sd
  }

  # 用 0 填充 NA（避免绘图崩溃）
  vals[is.na(vals)] <- 0
  as.numeric(vals)
}


# ==============================================================================
# 4. 补全：.draw_track
#    根据 config$type 分发到不同的 circlize 绘图函数
# ==============================================================================
.draw_track <- function(
  sector_data,
  config,
  group_col,
  n_groups,
  row_order_info,
  plot_params
) {
  col_name <- config$column
  track_type <- config$type %||% "points"
  color <- config$color %||% "#333333"
  group <- CELL_META$sector.index

  # 按照 row_order_info 排序
  if (!is.null(row_order_info[[group]])) {
    ordered_idx <- row_order_info[[group]]
    # 取与当前 sector_data 行对应的排序
    local_order <- match(ordered_idx, as.integer(rownames(sector_data)))
    local_order <- local_order[!is.na(local_order)]
    if (length(local_order) > 0) sector_data <- sector_data[local_order, ]
  }

  n_feat <- nrow(sector_data)
  if (n_feat == 0) {
    return(invisible(NULL))
  }

  x_coords <- seq_len(n_feat) - 0.5 # 每个 feature 的 x 中心
  y_vals <- sector_data[[col_name]]

  xlim <- CELL_META$cell.xlim
  ylim <- CELL_META$cell.ylim

  switch(
    track_type,

    # ---- 散点 ----------------------------------------------------------------
    points = {
      col_vec <- if (length(color) > 1) color[seq_len(n_feat)] else color
      circos.points(
        x_coords,
        y_vals,
        col = col_vec,
        pch = config$pch %||% 16,
        cex = plot_params$track_point_cex %||% 0.8
      )
    },

    # ---- 折线 ----------------------------------------------------------------
    line = {
      circos.lines(
        x_coords,
        y_vals,
        col = color,
        lwd = plot_params$track_line_lwd %||% 1.5
      )
    },

    # ---- 柱状 ----------------------------------------------------------------
    bars = {
      col_vec <- if (length(color) > 1) color[seq_len(n_feat)] else color
      for (k in seq_len(n_feat)) {
        circos.rect(
          x_coords[k] - 0.4,
          0,
          x_coords[k] + 0.4,
          y_vals[k],
          col = col_vec[k],
          border = NA
        )
      }
    },

    # ---- 棒棒糖 --------------------------------------------------------------
    lollipop = {
      col_vec <- if (length(color) > 1) color[seq_len(n_feat)] else color
      for (k in seq_len(n_feat)) {
        circos.lines(
          c(x_coords[k], x_coords[k]),
          c(0, y_vals[k]),
          col = col_vec[k],
          lwd = 1
        )
        circos.points(
          x_coords[k],
          y_vals[k],
          col = col_vec[k],
          pch = 19,
          cex = 0.9
        )
      }
    },

    # ---- 面积 ----------------------------------------------------------------
    area = {
      circos.lines(
        x_coords,
        y_vals,
        col = adjustcolor(color, alpha.f = 0.6),
        lwd = 1,
        area = TRUE,
        baseline = 0
      )
    },

    # ---- p 值（-log10 dots，分级着色）---------------------------------------
    pvalue = {
      thresholds <- c(config$p1 %||% 2, config$p2 %||% 3) # -log10 scale
      col_low <- config$col_low %||% "grey60"
      col_mid <- config$col_mid %||% "#FFA500"
      col_high <- config$col_high %||% "#E63946"
      col_vec <- ifelse(
        y_vals >= thresholds[2],
        col_high,
        ifelse(y_vals >= thresholds[1], col_mid, col_low)
      )
      circos.points(
        x_coords,
        y_vals,
        col = col_vec,
        pch = 16,
        cex = plot_params$track_point_cex %||% 0.8
      )
      # 参考线
      for (thr in thresholds) {
        circos.lines(xlim, c(thr, thr), col = "grey50", lty = 2, lwd = 0.8)
      }
    },

    # ---- 默认回退 ------------------------------------------------------------
    {
      warning(sprintf("Unknown track type '%s', using points.", track_type))
      circos.points(x_coords, y_vals, col = color, pch = 16)
    }
  )

  invisible(NULL)
}


# ==============================================================================
# 5. 补全：.calculate_feature_positions
#    为每个 feature 计算其在对应 sector 中的 x 坐标（用于 circos.link）
# ==============================================================================
.calculate_feature_positions <- function(
  data,
  group_col,
  feature_col,
  row_order_info
) {
  result <- lapply(unique(data[[group_col]]), function(g) {
    idx_in_data <- which(data[[group_col]] == g)
    ordered_idx <- row_order_info[[g]] %||% idx_in_data
    # 对齐
    ordered_idx <- ordered_idx[ordered_idx %in% idx_in_data]
    features <- data[[feature_col]][ordered_idx]
    n <- length(features)
    data.frame(
      feature = features,
      sector = g,
      x_pos = seq_len(n) - 0.5,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}


# ==============================================================================
# 6. 补全：.add_heatmap_legend
#    在图形右侧绘制渐变色条图例（使用 base graphics）
# ==============================================================================
.add_heatmap_legend <- function(
  col_fun,
  value_range,
  title = "Value",
  x = 0.85,
  y = 0.85,
  width = 0.08,
  height = 0.5
) {
  n_steps <- 100
  breaks <- seq(value_range[1], value_range[2], length.out = n_steps)
  colors <- col_fun(breaks)

  # 转换到 npc 坐标（图形设备左下角为原点）
  legend_x <- grconvertX(x, from = "npc", to = "user")
  legend_y <- grconvertY(y, from = "npc", to = "user")
  bar_w <- grconvertX(x + width, from = "npc", to = "user") - legend_x
  bar_h <- grconvertY(y - height, from = "npc", to = "user") - legend_y

  # 绘制色条（从下到上）
  step_h <- bar_h / n_steps
  for (i in seq_len(n_steps)) {
    rect(
      xleft = legend_x,
      ybottom = legend_y + (i - 1) * step_h,
      xright = legend_x + bar_w,
      ytop = legend_y + i * step_h,
      col = colors[i],
      border = NA
    )
  }

  # 外框
  rect(
    legend_x,
    legend_y,
    legend_x + bar_w,
    legend_y + bar_h,
    col = NA,
    border = "black",
    lwd = 0.8
  )

  # 刻度标签
  n_labels <- 5
  label_at <- seq(value_range[1], value_range[2], length.out = n_labels)
  for (i in seq_len(n_labels)) {
    frac <- (i - 1) / (n_labels - 1)
    yy <- legend_y + frac * bar_h
    text(
      legend_x + bar_w * 1.15,
      yy,
      labels = sprintf("%.2f", label_at[i]),
      cex = 0.7,
      adj = c(0, 0.5)
    )
  }

  # 标题
  text(
    legend_x + bar_w / 2,
    legend_y + bar_h + abs(bar_h) * 0.06,
    labels = title,
    cex = 0.75,
    font = 2,
    adj = c(0.5, 0)
  )

  invisible(NULL)
}


# ==============================================================================
# 7. 验证辅助函数（保持原有）
# ==============================================================================
.validate_track_heights <- function(
  track_configs,
  heatmap_height = 0.2,
  show_colnames = TRUE,
  show_group_labels = TRUE
) {
  fixed_space <- heatmap_height
  if (show_colnames) {
    fixed_space <- fixed_space + 0.05
  }
  if (show_group_labels) {
    fixed_space <- fixed_space + 0.08
  }

  track_heights <- sapply(track_configs, function(x) x$height %||% 0.1)
  total_track_height <- sum(track_heights)
  available_space <- 0.8 - fixed_space

  if (total_track_height > available_space) {
    warning(sprintf(
      "Track heights (%.2f) exceed available space (%.2f). Auto-adjusting...",
      total_track_height,
      available_space
    ))
    scale_factor <- available_space / total_track_height * 0.95
    for (i in seq_along(track_configs)) {
      original_height <- track_configs[[i]]$height %||% 0.1
      track_configs[[i]]$height <- original_height * scale_factor
      message(sprintf(
        "  Track %d (%s): %.3f -> %.3f",
        i,
        track_configs[[i]]$label %||% paste0("track", i),
        original_height,
        track_configs[[i]]$height
      ))
    }
  }
  return(track_configs)
}

.suggest_track_heights <- function(n_tracks, total_available = 0.5) {
  if (n_tracks == 0) {
    return(numeric(0))
  }
  if (n_tracks <= 3) {
    base_height <- 0.12
  } else if (n_tracks <= 5) {
    base_height <- total_available / n_tracks * 0.9
  } else {
    base_height <- total_available / n_tracks * 0.85
  }
  base_height <- max(base_height, 0.04)
  rep(base_height, n_tracks)
}

.validate_circos_input <- function(data, group_col, feature_col, value_cols) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }
  if (nrow(data) == 0) {
    stop("'data' is empty")
  }
  missing_cols <- setdiff(c(group_col, feature_col, value_cols), colnames(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns: %s", paste(missing_cols, collapse = ", ")))
  }
  if (anyDuplicated(data[[feature_col]])) {
    stop(sprintf("Duplicate values in '%s'", feature_col))
  }
  non_numeric <- value_cols[!sapply(data[value_cols], is.numeric)]
  if (length(non_numeric) > 0) {
    stop(sprintf(
      "Non-numeric columns: %s",
      paste(non_numeric, collapse = ", ")
    ))
  }
  invisible(TRUE)
}

.prepare_group_colors <- function(data, group_col, group_palette = NULL) {
  data[[group_col]] <- as.factor(data[[group_col]])
  group_levels <- levels(data[[group_col]])
  n_groups <- length(group_levels)

  # 用户提供命名向量
  if (
    !is.null(group_palette) &&
      is.character(group_palette) &&
      !is.null(names(group_palette))
  ) {
    data[[group_col]] <- factor(
      data[[group_col]],
      levels = names(group_palette)
    )
    missing_levels <- setdiff(group_levels, names(group_palette))
    if (length(missing_levels) > 0) {
      extra_colors <- get_palette_colors("npg", length(missing_levels))
      names(extra_colors) <- missing_levels
      group_palette <- c(group_palette, extra_colors)
    }
    message("Using user-defined color palette (named vector)")
    return(list(data = data, colors = group_palette))
  }

  # 用户提供调色板名称（单字符串）
  if (
    !is.null(group_palette) &&
      is.character(group_palette) &&
      length(group_palette) == 1
  ) {
    colors <- get_palette_colors(group_palette, n_groups)
    group_colors <- setNames(colors, group_levels)
    message(sprintf("Using built-in palette: '%s'", group_palette))
    return(list(data = data, colors = group_colors))
  }

  # 用户提供无名称颜色向量
  if (
    !is.null(group_palette) &&
      is.character(group_palette) &&
      length(group_palette) > 1
  ) {
    if (length(group_palette) < n_groups) {
      warning(sprintf(
        "Provided %d colors but have %d groups. Interpolating...",
        length(group_palette),
        n_groups
      ))
      colors <- colorRampPalette(group_palette)(n_groups)
    } else {
      colors <- group_palette[seq_len(n_groups)]
    }
    group_colors <- setNames(colors, group_levels)
    message("Using user-defined color vector")
    return(list(data = data, colors = group_colors))
  }

  # 默认
  if (n_groups <= 12) {
    colors <- RColorBrewer::brewer.pal(max(3, n_groups), "Set3")[seq_len(
      n_groups
    )]
  } else {
    colors <- scales::hue_pal()(n_groups)
  }
  message("Using default color palette")
  return(list(data = data, colors = setNames(colors, group_levels)))
}

.prepare_heatmap_colors <- function(
  heatmap_col = NULL,
  value_range = c(0, 0.5)
) {
  if (is.function(heatmap_col)) {
    message("Using user-defined color mapping function")
    return(heatmap_col)
  }
  if (is.character(heatmap_col) && length(heatmap_col) == 1) {
    colors <- get_palette_colors(heatmap_col, 11)
    # 修复：先构建函数，再打印 message（原代码此处 heatmap_col 已被覆盖为函数）
    palette_name <- heatmap_col
    heatmap_col <- circlize::colorRamp2(
      seq(value_range[1], value_range[2], length.out = 11),
      colors
    )
    message(sprintf("Using built-in heatmap palette: '%s'", palette_name))
    return(heatmap_col)
  }
  if (is.character(heatmap_col) && length(heatmap_col) > 1) {
    n_colors <- length(heatmap_col)
    heatmap_col <- circlize::colorRamp2(
      seq(value_range[1], value_range[2], length.out = n_colors),
      heatmap_col
    )
    message(sprintf("Using user-defined color vector (%d colors)", n_colors))
    return(heatmap_col)
  }
  # 默认
  colors <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
  heatmap_col <- circlize::colorRamp2(
    seq(value_range[1], value_range[2], length.out = length(colors)),
    colors
  )
  message("Using default heatmap colors (RdBu)")
  return(heatmap_col)
}


# ==============================================================================
# 8. 主函数（修复 finally 无条件 circos.clear() 的问题）
# ==============================================================================
plot_advanced_circos <- function(
  data,
  group_col = "Group",
  feature_col = "Gene",
  value_cols,
  group_palette = NULL,
  heatmap_col = NULL,
  track_configs = list(),
  links_data = NULL,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gap_degree = 90,
  gap_between = 2,
  start_degree = 135,
  value_range = c(0, 0.5),
  plot_params = list(),
  group_label_params = list(),
  show_legend = TRUE,
  legend_params = list(),
  auto_adjust_heights = TRUE,
  ...
) {
  tryCatch(circos.clear(), error = function(e) NULL)

  result <- tryCatch(
    {
      # ---- 1. 输入验证 -------------------------------------------------------
      .validate_circos_input(data, group_col, feature_col, value_cols)

      # ---- 参数合并 -----------------------------------------------------------
      default_plot_params <- list(
        heatmap_rownames_cex = 0.7,
        heatmap_cell_border = "white",
        heatmap_cell_lwd = 0.5,
        heatmap_bg_border = "black",
        heatmap_bg_lwd = 1,
        colnames_cex = 0.8,
        track_point_cex = 1.0,
        track_line_lwd = 1.5,
        track_label_cex = 0.8,
        group_label_cex = 1.2,
        group_bg_alpha = 0.5,
        heatmap_track_height = 0.2
      )
      plot_params <- utils::modifyList(default_plot_params, plot_params)

      default_group_label_params <- list(
        cex = NULL,
        font = 2,
        col = "white",
        bg_alpha = 0.5,
        facing = "bending.inside",
        niceFacing = TRUE,
        track_height = 0.08,
        show = TRUE,
        auto_adjust = TRUE
      )
      group_label_params <- utils::modifyList(
        default_group_label_params,
        group_label_params
      )

      default_legend_params <- list(
        x = 0.85,
        y = 0.85,
        width = 0.08,
        height = 0.5,
        title = "Value"
      )
      legend_params <- utils::modifyList(default_legend_params, legend_params)

      # ---- 1.1 聚类检查 -------------------------------------------------------
      user_args <- list(...)
      will_cluster <- isTRUE(user_args$cluster)
      if (will_cluster) {
        if (plot_params$heatmap_track_height > 0.18) {
          original_height <- plot_params$heatmap_track_height
          plot_params$heatmap_track_height <- 0.15
          message(sprintf(
            "Clustering: reducing heatmap height %.2f -> %.2f",
            original_height,
            plot_params$heatmap_track_height
          ))
        }
        if (show_legend) show_legend <- FALSE
      }

      # ---- 1.2 自适应字体 ----------------------------------------------------
      n_tracks <- length(track_configs)
      if (group_label_params$auto_adjust && is.null(group_label_params$cex)) {
        group_label_params$cex <- dplyr::case_when(
          n_tracks == 0 ~ 1.2,
          n_tracks <= 2 ~ 1.1,
          n_tracks <= 4 ~ 1.0,
          n_tracks <= 6 ~ 0.9,
          TRUE ~ 0.8
        )
        message(sprintf(
          "Auto-adjusted group label size to %.1f (%d tracks)",
          group_label_params$cex,
          n_tracks
        ))
      }

      # ---- 1.3 轨道高度验证 --------------------------------------------------
      if (auto_adjust_heights && length(track_configs) > 0) {
        extra_space <- if (will_cluster) {
          (user_args$dend.track.height %||% 0.1) + 0.02
        } else {
          0
        }
        track_configs <- .validate_track_heights(
          track_configs,
          heatmap_height = plot_params$heatmap_track_height + extra_space,
          show_colnames = show_colnames,
          show_group_labels = group_label_params$show
        )
      }

      # ---- 2. 数据准备 -------------------------------------------------------
      color_result <- .prepare_group_colors(data, group_col, group_palette)
      data <- color_result$data
      group_colors <- color_result$colors
      heatmap_col <- .prepare_heatmap_colors(heatmap_col, value_range)

      data_normalized <- data %>%
        dplyr::group_by(!!sym(group_col)) %>%
        dplyr::mutate(dplyr::across(
          dplyr::all_of(value_cols),
          ~ scales::rescale(., to = value_range)
        )) %>%
        dplyr::ungroup()

      data_matrix <- data_normalized %>%
        dplyr::select(dplyr::all_of(value_cols)) %>%
        as.matrix()
      rownames(data_matrix) <- data_normalized[[feature_col]]

      # ---- 3. circos 初始化 --------------------------------------------------
      circos.clear()
      group_levels <- levels(data_normalized[[group_col]])
      n_groups <- length(group_levels)

      gaps_vector <- rep(gap_between, n_groups)
      gaps_vector[n_groups] <- gap_degree

      track_margin <- if (n_tracks > 5 || will_cluster) {
        c(0.005, 0.005)
      } else {
        c(0.01, 0.01)
      }

      circos.par(
        start.degree = start_degree,
        gap.after = gaps_vector,
        track.margin = track_margin,
        cell.padding = c(0, 0, 0, 0)
      )

      circos.heatmap.initialize(
        data_matrix,
        split = data_normalized[[group_col]]
      )

      # ---- 4. 热图 -----------------------------------------------------------
      default_heatmap_args <- list(
        mat = data_matrix,
        split = data_normalized[[group_col]],
        col = heatmap_col,
        track.height = plot_params$heatmap_track_height,
        cluster = FALSE,
        bg.border = plot_params$heatmap_bg_border,
        bg.lwd = plot_params$heatmap_bg_lwd,
        cell.border = plot_params$heatmap_cell_border,
        cell.lwd = plot_params$heatmap_cell_lwd,
        rownames.side = ifelse(show_rownames, "outside", "none"),
        rownames.cex = plot_params$heatmap_rownames_cex
      )
      final_heatmap_args <- utils::modifyList(default_heatmap_args, user_args)
      do.call(circos.heatmap, final_heatmap_args)
      clustered <- isTRUE(final_heatmap_args$cluster)
      row_order_info <- .extract_row_order(
        data_normalized,
        group_col,
        clustered
      )

      # ---- 5. 列名轨道 -------------------------------------------------------
      if (show_colnames) {
        circos.track(
          track.index = get.current.track.index(),
          bg.border = NA,
          panel.fun = function(x, y) {
            if (CELL_META$sector.numeric.index == n_groups) {
              cn <- colnames(data_matrix)
              n <- length(cn)
              cell_h <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
              y_coords <- seq(
                CELL_META$cell.ylim[1] + cell_h / 2,
                CELL_META$cell.ylim[2] - cell_h / 2,
                length.out = n
              )
              for (i in seq_len(n)) {
                circos.lines(
                  c(
                    CELL_META$cell.xlim[2],
                    CELL_META$cell.xlim[2] + convert_x(1, "mm")
                  ),
                  c(y_coords[i], y_coords[i]),
                  col = "black",
                  lwd = 1
                )
              }
              circos.text(
                rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"),
                y_coords,
                cn,
                cex = plot_params$colnames_cex,
                adj = c(0, 0.5),
                facing = "inside"
              )
            }
          }
        )
      }

      # ---- 6. 附加轨道 -------------------------------------------------------
      if (length(track_configs) > 0) {
        for (i in seq_along(track_configs)) {
          config <- track_configs[[i]]
          required_fields <- c("column", "type", "height", "label")
          missing_fields <- setdiff(required_fields, names(config))
          if (length(missing_fields) > 0) {
            warning(sprintf(
              "Track %d missing: %s",
              i,
              paste(missing_fields, collapse = ", ")
            ))
            next
          }
          if (!config$column %in% colnames(data_normalized)) {
            warning(sprintf("Track column '%s' not found", config$column))
            next
          }

          track_data <- .process_track_data(data_normalized, config)
          data_normalized[[config$column]] <- track_data

          val_range <- range(track_data, na.rm = TRUE)
          range_span <- diff(val_range)
          if (range_span == 0) {
            val_range <- c(val_range[1] - 0.1, val_range[2] + 0.1)
          } else if (config$type %in% c("bars", "lollipop", "area")) {
            if (val_range[1] > 0) {
              val_range[1] <- 0
            }
            val_range[2] <- val_range[2] + range_span * 0.15
          } else {
            val_range[1] <- val_range[1] - range_span * 0.1
            val_range[2] <- val_range[2] + range_span * 0.1
          }

          tryCatch(
            circos.track(
              ylim = val_range,
              track.height = config$height,
              bg.border = "grey90",
              panel.fun = function(x, y) {
                sector_data <- data_normalized[
                  data_normalized[[group_col]] == CELL_META$sector.index,
                ]
                .draw_track(
                  sector_data,
                  config,
                  group_col,
                  n_groups,
                  row_order_info,
                  plot_params
                )
              }
            ),
            error = function(e) {
              warning(sprintf("Track %d failed: %s", i, e$message))
            }
          )
        }
      }

      # ---- 7. 分组标签 -------------------------------------------------------
      if (group_label_params$show) {
        circos.track(
          ylim = c(0, 1),
          track.height = group_label_params$track_height,
          bg.border = NA,
          bg.col = adjustcolor(
            group_colors[group_levels],
            alpha.f = group_label_params$bg_alpha
          ),
          panel.fun = function(x, y) {
            circos.text(
              CELL_META$xcenter,
              0.5,
              CELL_META$sector.index,
              facing = group_label_params$facing,
              niceFacing = group_label_params$niceFacing,
              cex = group_label_params$cex,
              font = group_label_params$font,
              col = group_label_params$col,
              adj = c(0.5, 0.5)
            )
          }
        )
      }

      # ---- 8. 连线 -----------------------------------------------------------
      if (!is.null(links_data) && nrow(links_data) > 0) {
        if (all(c("from", "to") %in% colnames(links_data))) {
          feature_mapping <- .calculate_feature_positions(
            data_normalized,
            group_col,
            feature_col,
            row_order_info
          )
          links_drawn <- 0
          for (i in seq_len(nrow(links_data))) {
            from_meta <- feature_mapping[
              feature_mapping$feature == links_data$from[i],
            ]
            to_meta <- feature_mapping[
              feature_mapping$feature == links_data$to[i],
            ]
            if (nrow(from_meta) == 1 && nrow(to_meta) == 1) {
              link_col <- if ("color" %in% colnames(links_data)) {
                links_data$color[i]
              } else {
                "#00000050"
              }
              tryCatch(
                {
                  circos.link(
                    sector.index1 = as.character(from_meta$sector),
                    point1 = from_meta$x_pos,
                    sector.index2 = as.character(to_meta$sector),
                    point2 = to_meta$x_pos,
                    col = link_col,
                    lwd = 1.5,
                    directional = 1,
                    arr.length = 0.2
                  )
                  links_drawn <- links_drawn + 1
                },
                error = function(e) NULL
              )
            }
          }
          if (links_drawn > 0) message(sprintf("Drew %d links", links_drawn))
        }
      }

      # ---- 9. 图例 -----------------------------------------------------------
      if (show_legend) {
        tryCatch(
          .add_heatmap_legend(
            heatmap_col,
            value_range,
            title = legend_params$title,
            x = legend_params$x,
            y = legend_params$y,
            width = legend_params$width,
            height = legend_params$height
          ),
          error = function(e) warning(sprintf("Legend failed: %s", e$message))
        )
      }

      # ---- 10. 返回元数据 ---------------------------------------------------
      list(
        data_matrix = data_matrix,
        data_normalized = data_normalized,
        group_colors = group_colors,
        group_levels = group_levels,
        n_groups = n_groups,
        clustered = clustered,
        row_order_info = row_order_info,
        feature_mapping = if (!is.null(links_data)) {
          .calculate_feature_positions(
            data_normalized,
            group_col,
            feature_col,
            row_order_info
          )
        } else {
          NULL
        },
        plot_params = plot_params,
        group_label_params = group_label_params,
        track_configs = track_configs
      )
    },
    error = function(e) {
      # 仅在出错时清除（修复：原脚本 finally 无条件 clear 会擦掉正确绘制的图）
      tryCatch(circos.clear(), error = function(e2) NULL)
      stop(sprintf("plot_advanced_circos Error: %s", e$message))
    }
    # 注意：移除了 finally { circos.clear() }
    # circlize 的状态由 circos.clear() 在下次调用前主动重置即可
  )

  return(invisible(result))
}


# ==============================================================================
# 9. 模拟数据测试（复杂场景）
# ==============================================================================
# 场景：scRNA-seq 差异基因分析结果 + 细胞类型分组
#   - 6 种细胞类型（Group）
#   - 60 个基因（Gene），每组 10 个
#   - 4 条热图值列：avg_log2FC, pct.1, pct.2, specificity
#   - 3 条附加轨道：
#       * 柱状图（-log10 p-value）
#       * 棒棒糖图（基因表达量 avg_exp）
#       * 散点图（模块评分 module_score）
#   - 若干 links（基因间相关性连线）
# ==============================================================================

set.seed(42)

n_groups <- 6
n_per_group <- 10
n_genes <- n_groups * n_per_group

cell_types <- c("CD4_T", "CD8_T", "NK", "B_cell", "Monocyte", "DC")
genes <- paste0("Gene_", sprintf("%02d", seq_len(n_genes)))

# 构造基础数据框
sim_data <- data.frame(
  Group = rep(cell_types, each = n_per_group),
  Gene = genes,
  # 热图值：log2FC（组间有差异）
  avg_log2FC = c(
    rnorm(10, mean = 2.0, sd = 0.5), # CD4_T  高表达
    rnorm(10, mean = 1.5, sd = 0.6), # CD8_T
    rnorm(10, mean = 0.5, sd = 0.8), # NK
    rnorm(10, mean = -0.5, sd = 0.7), # B_cell
    rnorm(10, mean = -1.0, sd = 0.6), # Monocyte
    rnorm(10, mean = -1.8, sd = 0.5) # DC  低表达
  ),
  pct.1 = runif(n_genes, 0.3, 0.95),
  pct.2 = runif(n_genes, 0.05, 0.5),
  specificity = runif(n_genes, 0.1, 1.0),
  # 附加轨道数据
  pvalue = runif(n_genes, 1e-20, 0.05), # 显著性
  avg_exp = rgamma(n_genes, shape = 2, rate = 0.5),
  module_score = rnorm(n_genes, mean = 0.1, sd = 0.3),
  stringsAsFactors = FALSE
)

# 轨道配置
track_configs_demo <- list(
  list(
    column = "pvalue",
    type = "pvalue", # 自动做 -log10 变换 + 分级着色
    height = 0.09,
    label = "-log10(pval)",
    transform = "log10", # .process_track_data 会做 -log10
    col_low = "grey75",
    col_mid = "#FFA500",
    col_high = "#E63946",
    p1 = 2,
    p2 = 4 # -log10 阈值
  ),
  list(
    column = "avg_exp",
    type = "lollipop",
    height = 0.09,
    label = "Avg Expression",
    color = "#4DBBD5"
  ),
  list(
    column = "module_score",
    type = "points",
    height = 0.08,
    label = "Module Score",
    color = "#E64B35"
  )
)

# 连线数据（随机选 8 对基因建立有向连接）
from_genes <- sample(genes[1:30], 8)
to_genes <- sample(genes[31:60], 8)
link_colors <- paste0(sample(
  c("#E64B3580", "#4DBBD580", "#00A08780"),
  8,
  replace = TRUE
))
links_demo <- data.frame(
  from = from_genes,
  to = to_genes,
  color = link_colors,
  stringsAsFactors = FALSE
)

# 自定义分组配色
group_palette_demo <- c(
  CD4_T = "#E64B35",
  CD8_T = "#4DBBD5",
  NK = "#00A087",
  B_cell = "#3C5488",
  Monocyte = "#F39B7F",
  DC = "#8491B4"
)

# ---- 执行绘图 ----------------------------------------------------------------
message("\n===== 开始绘制环状热图测试 =====\n")

res <- plot_advanced_circos(
  data = sim_data,
  group_col = "Group",
  feature_col = "Gene",
  value_cols = c("avg_log2FC", "pct.1", "pct.2", "specificity"),
  group_palette = group_palette_demo,
  heatmap_col = NULL, # 使用默认 RdBu
  track_configs = track_configs_demo,
  links_data = links_demo,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gap_degree = 80,
  gap_between = 2,
  start_degree = 135,
  value_range = c(0, 1),
  show_legend = TRUE,
  legend_params = list(title = "log2FC (scaled)", x = 0.87, y = 0.82),
  auto_adjust_heights = TRUE
)

message("\n===== 绘图完成 =====")
message(sprintf("Groups: %s", paste(res$group_levels, collapse = ", ")))
message(sprintf(
  "Features: %d | Heatmap cols: %d | Tracks: %d",
  nrow(res$data_matrix),
  ncol(res$data_matrix),
  length(res$track_configs)
))
