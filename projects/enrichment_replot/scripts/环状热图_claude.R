library(circlize)
library(RColorBrewer)
library(dplyr)
library(scales)

# ==============================================================================
# Helper Functions (模块化)
# ==============================================================================
#' Validate and Adjust Track Heights
#'
#' 确保所有轨道加上热图不会超过可用空间
#' @keywords internal
.validate_track_heights <- function(
  track_configs,
  heatmap_height = 0.2,
  show_colnames = TRUE,
  show_group_labels = TRUE
) {
  # 计算固定占用的空间
  fixed_space <- heatmap_height # 热图
  if (show_colnames) {
    fixed_space <- fixed_space + 0.05
  } # 列名轨道（估算）
  if (show_group_labels) {
    fixed_space <- fixed_space + 0.08
  } # 分组标签轨道

  # 计算轨道总高度
  track_heights <- sapply(track_configs, function(x) x$height %||% 0.1)
  total_track_height <- sum(track_heights)

  # 可用空间（保守估计，留20%缓冲）
  available_space <- 0.8 - fixed_space

  if (total_track_height > available_space) {
    warning(sprintf(
      "Track heights (%.2f) exceed available space (%.2f). Auto-adjusting...",
      total_track_height,
      available_space
    ))

    # 按比例缩小所有轨道
    scale_factor <- available_space / total_track_height * 0.95 # 再留5%安全边际

    for (i in seq_along(track_configs)) {
      original_height <- track_configs[[i]]$height
      new_height <- original_height * scale_factor
      track_configs[[i]]$height <- new_height

      message(sprintf(
        "  Track %d (%s): %.3f -> %.3f",
        i,
        track_configs[[i]]$label,
        original_height,
        new_height
      ))
    }
  }

  return(track_configs)
}

#' Calculate Optimal Track Heights
#'
#' 根据轨道数量智能分配高度
#' @keywords internal
.suggest_track_heights <- function(n_tracks, total_available = 0.5) {
  if (n_tracks == 0) {
    return(numeric(0))
  }

  # 根据轨道数量使用不同策略
  if (n_tracks <= 3) {
    # 少量轨道：可以给较大空间
    base_height <- 0.12
  } else if (n_tracks <= 5) {
    # 中等数量：均匀分配
    base_height <- total_available / n_tracks * 0.9
  } else {
    # 大量轨道：压缩高度
    base_height <- total_available / n_tracks * 0.85
  }

  # 最小高度限制
  min_height <- 0.04
  base_height <- max(base_height, min_height)

  return(rep(base_height, n_tracks))
}

#' Validate Input Data
#' @keywords internal
.validate_circos_input <- function(data, group_col, feature_col, value_cols) {
  # Check data frame
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }
  if (nrow(data) == 0) {
    stop("'data' is empty")
  }

  # Check required columns
  missing_cols <- setdiff(c(group_col, feature_col, value_cols), colnames(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing columns in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # Check for duplicates
  if (anyDuplicated(data[[feature_col]])) {
    stop(sprintf("Duplicate values found in '%s' column", feature_col))
  }

  # Check value columns are numeric
  non_numeric <- value_cols[!sapply(data[value_cols], is.numeric)]
  if (length(non_numeric) > 0) {
    stop(sprintf(
      "Non-numeric value columns: %s",
      paste(non_numeric, collapse = ", ")
    ))
  }

  invisible(TRUE)
}


#' Prepare Group Colors (增强版 - 支持多种输入方式)
#' @keywords internal
.prepare_group_colors <- function(data, group_col, group_palette = NULL) {
  data[[group_col]] <- as.factor(data[[group_col]])
  group_levels <- levels(data[[group_col]])
  n_groups <- length(group_levels)

  # ========================================================================
  # 情况1：用户提供了命名向量（最高优先级）
  # ========================================================================
  if (
    !is.null(group_palette) &&
      is.character(group_palette) &&
      !is.null(names(group_palette))
  ) {
    # 这是用户自定义的命名颜色向量
    # 例如：c("CD4_T" = "#FF0000", "CD8_T" = "#00FF00")

    data[[group_col]] <- factor(
      data[[group_col]],
      levels = names(group_palette)
    )
    missing_levels <- setdiff(group_levels, names(group_palette))

    if (length(missing_levels) > 0) {
      warning(sprintf(
        "Groups without colors in palette: %s. Using default colors for them.",
        paste(missing_levels, collapse = ", ")
      ))

      # 为缺失的分组添加默认颜色
      extra_colors <- get_palette_colors("npg", length(missing_levels))
      names(extra_colors) <- missing_levels
      group_palette <- c(group_palette, extra_colors)
    }

    message("Using user-defined color palette (named vector)")
    return(list(data = data, colors = group_palette))
  }

  # ========================================================================
  # 情况2：用户提供了调色板名称（字符串，无名称）
  # ========================================================================
  if (
    !is.null(group_palette) &&
      is.character(group_palette) &&
      length(group_palette) == 1
  ) {
    # 这是内置调色板名称
    # 例如："npg", "lancet", "jco"

    colors <- get_palette_colors(group_palette, n_groups)
    group_colors <- setNames(colors, group_levels)

    message(sprintf("Using built-in palette: '%s'", group_palette))
    return(list(data = data, colors = group_colors))
  }

  # ========================================================================
  # 情况3：用户提供了无名称的颜色向量
  # ========================================================================
  if (
    !is.null(group_palette) &&
      is.character(group_palette) &&
      length(group_palette) > 1
  ) {
    # 这是用户自定义的颜色向量（无名称）
    # 例如：c("#FF0000", "#00FF00", "#0000FF")

    if (length(group_palette) < n_groups) {
      warning(sprintf(
        "Provided %d colors but have %d groups. Interpolating...",
        length(group_palette),
        n_groups
      ))
      colors <- colorRampPalette(group_palette)(n_groups)
    } else {
      colors <- group_palette[1:n_groups]
    }

    group_colors <- setNames(colors, group_levels)

    message("Using user-defined color vector")
    return(list(data = data, colors = group_colors))
  }

  # ========================================================================
  # 情况4：NULL - 使用默认配色
  # ========================================================================
  if (n_groups <= 12) {
    colors <- RColorBrewer::brewer.pal(max(3, n_groups), "Set3")[seq_len(
      n_groups
    )]
  } else {
    colors <- scales::hue_pal()(n_groups)
  }

  group_colors <- setNames(colors, group_levels)

  message("Using default color palette")
  return(list(data = data, colors = group_colors))
}

#' Prepare Heatmap Colors (增强版 - 支持多种输入)
#' @keywords internal
.prepare_heatmap_colors <- function(
  heatmap_col = NULL,
  value_range = c(0, 0.5)
) {
  # ========================================================================
  # 情况1：用户提供了 colorRamp2 函数（最高优先级）
  # ========================================================================
  if (is.function(heatmap_col)) {
    message("Using user-defined color mapping function")
    return(heatmap_col)
  }

  # ========================================================================
  # 情况2：用户提供了调色板名称
  # ========================================================================
  if (is.character(heatmap_col) && length(heatmap_col) == 1) {
    colors <- get_palette_colors(heatmap_col, 11)
    heatmap_col <- circlize::colorRamp2(
      seq(value_range[1], value_range[2], length.out = 11),
      colors
    )

    message(sprintf("Using built-in heatmap palette: '%s'", heatmap_col))
    return(heatmap_col)
  }

  # ========================================================================
  # 情况3：用户提供了颜色向量
  # ========================================================================
  if (is.character(heatmap_col) && length(heatmap_col) > 1) {
    n_colors <- length(heatmap_col)
    heatmap_col <- circlize::colorRamp2(
      seq(value_range[1], value_range[2], length.out = n_colors),
      heatmap_col
    )

    message(sprintf("Using user-defined color vector (%d colors)", n_colors))
    return(heatmap_col)
  }

  # ========================================================================
  # 情况4：NULL - 使用默认
  # ========================================================================
  colors <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
  heatmap_col <- circlize::colorRamp2(
    seq(value_range[1], value_range[2], length.out = length(colors)),
    colors
  )

  message("Using default heatmap colors (RdBu)")
  return(heatmap_col)
}


# ==============================================================================
# 主函数保持不变，但删除之前的配色预处理部分
# ==============================================================================

plot_advanced_circos <- function(
  data,
  group_col = "Group",
  feature_col = "Gene",
  value_cols,
  group_palette = NULL, # 可以是：调色板名、命名向量、颜色向量、NULL
  heatmap_col = NULL, # 可以是：调色板名、颜色向量、函数、NULL
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

  tryCatch(
    {
      # ========================================================================
      # 1. Input Validation
      # ========================================================================
      .validate_circos_input(data, group_col, feature_col, value_cols)

      # Merge plot parameters
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

      # Group label parameters
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

      # ========================================================================
      # 1.1 检查聚类
      # ========================================================================
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

        if (show_legend) {
          show_legend <- FALSE
        }
      }

      # ========================================================================
      # 1.2 自适应字体大小
      # ========================================================================
      n_tracks <- length(track_configs)

      if (group_label_params$auto_adjust && is.null(group_label_params$cex)) {
        if (n_tracks == 0) {
          group_label_params$cex <- 1.2
        } else if (n_tracks <= 2) {
          group_label_params$cex <- 1.1
        } else if (n_tracks <= 4) {
          group_label_params$cex <- 1.0
        } else if (n_tracks <= 6) {
          group_label_params$cex <- 0.9
        } else {
          group_label_params$cex <- 0.8
        }

        message(sprintf(
          "Auto-adjusted group label size to %.1f (based on %d tracks)",
          group_label_params$cex,
          n_tracks
        ))
      }

      # ========================================================================
      # 1.3 验证轨道高度
      # ========================================================================
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

      # ========================================================================
      # 2. Data Preparation (使用增强的配色函数)
      # ========================================================================
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

      # ========================================================================
      # 3-11. 其余部分保持不变
      # ========================================================================

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

      if (show_colnames) {
        circos.track(
          track.index = get.current.track.index(),
          bg.border = NA,
          panel.fun = function(x, y) {
            if (CELL_META$sector.numeric.index == n_groups) {
              cn <- colnames(data_matrix)
              n <- length(cn)
              cell_height <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) /
                n
              y_coords <- seq(
                CELL_META$cell.ylim[1] + cell_height / 2,
                CELL_META$cell.ylim[2] - cell_height / 2,
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
          } else {
            if (config$type %in% c("bars", "lollipop", "area")) {
              if (val_range[1] > 0) {
                val_range[1] <- 0
              }
              val_range[2] <- val_range[2] + range_span * 0.15
            } else {
              val_range[1] <- val_range[1] - range_span * 0.1
              val_range[2] <- val_range[2] + range_span * 0.1
            }
          }

          tryCatch(
            {
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
              )
            },
            error = function(e) {
              warning(sprintf("Track %d failed: %s", i, e$message))
            }
          )
        }
      }

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

      if (!is.null(links_data) && nrow(links_data) > 0) {
        required_link_cols <- c("from", "to")
        if (all(required_link_cols %in% colnames(links_data))) {
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

      if (show_legend) {
        tryCatch(
          {
            .add_heatmap_legend(
              heatmap_col,
              value_range,
              title = legend_params$title,
              x = legend_params$x,
              y = legend_params$y,
              width = legend_params$width,
              height = legend_params$height
            )
          },
          error = function(e) {
            warning(sprintf("Legend failed: %s", e$message))
          }
        )
      }

      result <- list(
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

      return(invisible(result))
    },
    error = function(e) {
      tryCatch(circos.clear(), error = function(e2) NULL)
      stop(sprintf("Error: %s", e$message))
    },
    finally = {
      circos.clear()
    }
  )
}
