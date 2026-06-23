# ==============================================================================
# Advanced Circular Heatmap Plotting Function
# ==============================================================================
# Version: 2.0
# Date: 2025-02-14
# Author: Claude & User
# Description: Create advanced circular heatmaps with multiple tracks,
#              gene-gene links, clustering, and flexible color schemes
# ==============================================================================

# Required packages
required_packages <- c("circlize", "dplyr", "scales", "RColorBrewer")
missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  warning(sprintf(
    "Missing required packages: %s\nPlease install with: install.packages(c(%s))",
    paste(missing_packages, collapse = ", "),
    paste(sprintf("'%s'", missing_packages), collapse = ", ")
  ))
}

# Load required packages
suppressPackageStartupMessages({
  library(circlize)
  library(dplyr)
  library(scales)
  library(RColorBrewer)
})

# Check for optional packages
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  message(
    "Note: ComplexHeatmap package not found. Legend functionality will be limited."
  )
  message("Install with: BiocManager::install('ComplexHeatmap')")
}

# ==============================================================================
# Color Palette System
# ==============================================================================

#' Get Universal Color Palette
#'
#' Provides multiple scientific publication-grade color schemes
#' @param palette_name Name of the color palette
#' @param n Number of colors needed
#' @return Character vector of colors
#' @export
get_palette_colors <- function(palette_name = "npg", n = NULL) {
  palettes <- list(
    # Nature Publishing Group
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

    # AAAS (Science)
    aaas = c(
      "#3B4992",
      "#EE0000",
      "#008B45",
      "#631879",
      "#008280",
      "#BB0021",
      "#5F559B",
      "#A20056",
      "#808180",
      "#1B1919"
    ),

    # Lancet
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

    # NEJM
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

    # JAMA
    jama = c(
      "#374E55",
      "#DF8F44",
      "#00A1D5",
      "#B24745",
      "#79AF97",
      "#6A6599",
      "#80796B"
    ),

    # JCO (Journal of Clinical Oncology)
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

    # D3 Category 10
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

    # Tableau 10
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
    ),

    # Sequential - Blues
    blues = RColorBrewer::brewer.pal(9, "Blues"),

    # Sequential - Reds
    reds = RColorBrewer::brewer.pal(9, "Reds"),

    # Diverging - Red-Blue
    rdbu = rev(RColorBrewer::brewer.pal(11, "RdBu")),

    # Qualitative - Set1
    set1 = RColorBrewer::brewer.pal(9, "Set1"),

    # Qualitative - Set2
    set2 = RColorBrewer::brewer.pal(8, "Set2"),

    # Qualitative - Set3
    set3 = RColorBrewer::brewer.pal(12, "Set3"),

    # Pastel colors
    pastel = c(
      "#FBB4AE",
      "#B3CDE3",
      "#CCEBC5",
      "#DECBE4",
      "#FED9A6",
      "#FFFFCC",
      "#E5D8BD",
      "#FDDAEC",
      "#F2F2F2"
    )
  )

  selected_palette <- palettes[[tolower(palette_name)]]

  if (is.null(selected_palette)) {
    warning(sprintf("Palette '%s' not found. Using 'npg'", palette_name))
    selected_palette <- palettes$npg
  }

  if (is.null(n)) {
    return(selected_palette)
  }

  # If more colors needed than available, use interpolation
  if (n > length(selected_palette)) {
    return(colorRampPalette(selected_palette)(n))
  } else {
    return(selected_palette[1:n])
  }
}

#' List Available Palettes
#' @export
list_color_palettes <- function() {
  cat("Available color palettes:\n\n")
  cat("Scientific Journals:\n")
  cat("  - npg      : Nature Publishing Group\n")
  cat("  - aaas     : Science (AAAS)\n")
  cat("  - lancet   : The Lancet\n")
  cat("  - nejm     : New England Journal of Medicine\n")
  cat("  - jama     : JAMA\n")
  cat("  - jco      : Journal of Clinical Oncology\n\n")

  cat("General Purpose:\n")
  cat("  - d3       : D3.js Category 10\n")
  cat("  - tableau10: Tableau 10\n")
  cat("  - set1/2/3 : ColorBrewer qualitative\n")
  cat("  - pastel   : Pastel colors\n\n")

  cat("Sequential:\n")
  cat("  - blues    : Blue gradient\n")
  cat("  - reds     : Red gradient\n\n")

  cat("Diverging:\n")
  cat("  - rdbu     : Red-Blue diverging\n")
}

#' Preview Color Palette
#' @param palette_name Name of the palette
#' @param n Number of colors
#' @export
preview_palette <- function(palette_name = "npg", n = NULL) {
  colors <- get_palette_colors(palette_name, n)

  n_colors <- length(colors)
  par(mar = c(0, 0, 2, 0))
  plot(
    1:n_colors,
    rep(1, n_colors),
    pch = 15,
    cex = 8,
    col = colors,
    xlim = c(0.5, n_colors + 0.5),
    ylim = c(0.5, 1.5),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = sprintf("Palette: %s (%d colors)", palette_name, n_colors)
  )

  # Add color labels
  text(1:n_colors, rep(0.7, n_colors), colors, srt = 45, adj = 1, cex = 0.7)
}

# ==============================================================================
# Helper Functions
# ==============================================================================

# Null coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

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

#' Prepare Group Colors
#' @keywords internal
.prepare_group_colors <- function(data, group_col, group_palette = NULL) {
  data[[group_col]] <- as.factor(data[[group_col]])
  group_levels <- levels(data[[group_col]])
  n_groups <- length(group_levels)

  # Case 1: Named color vector (highest priority)
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
      warning(sprintf(
        "Groups without colors in palette: %s. Using default colors for them.",
        paste(missing_levels, collapse = ", ")
      ))

      extra_colors <- get_palette_colors("npg", length(missing_levels))
      names(extra_colors) <- missing_levels
      group_palette <- c(group_palette, extra_colors)
    }

    message("Using user-defined color palette (named vector)")
    return(list(data = data, colors = group_palette))
  }

  # Case 2: Palette name
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

  # Case 3: Unnamed color vector
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
      colors <- group_palette[1:n_groups]
    }

    group_colors <- setNames(colors, group_levels)

    message("Using user-defined color vector")
    return(list(data = data, colors = group_colors))
  }

  # Case 4: NULL - use default
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

#' Prepare Heatmap Colors
#' @keywords internal
.prepare_heatmap_colors <- function(
  heatmap_col = NULL,
  value_range = c(0, 0.5)
) {
  # Case 1: User-provided function
  if (is.function(heatmap_col)) {
    message("Using user-defined color mapping function")
    return(heatmap_col)
  }

  # Case 2: Palette name
  if (is.character(heatmap_col) && length(heatmap_col) == 1) {
    colors <- get_palette_colors(heatmap_col, 11)
    heatmap_col <- circlize::colorRamp2(
      seq(value_range[1], value_range[2], length.out = 11),
      colors
    )

    message(sprintf(
      "Using built-in heatmap palette: '%s'",
      deparse(substitute(heatmap_col))
    ))
    return(heatmap_col)
  }

  # Case 3: Color vector
  if (is.character(heatmap_col) && length(heatmap_col) > 1) {
    n_colors <- length(heatmap_col)
    heatmap_col <- circlize::colorRamp2(
      seq(value_range[1], value_range[2], length.out = n_colors),
      heatmap_col
    )

    message(sprintf("Using user-defined color vector (%d colors)", n_colors))
    return(heatmap_col)
  }

  # Case 4: NULL - use default
  colors <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
  heatmap_col <- circlize::colorRamp2(
    seq(value_range[1], value_range[2], length.out = length(colors)),
    colors
  )

  message("Using default heatmap colors (RdBu)")
  return(heatmap_col)
}

#' Validate and Adjust Track Heights
#' @keywords internal
.validate_track_heights <- function(
  track_configs,
  heatmap_height = 0.2,
  show_colnames = TRUE,
  show_group_labels = TRUE
) {
  # Calculate fixed space
  fixed_space <- heatmap_height
  if (show_colnames) {
    fixed_space <- fixed_space + 0.05
  }
  if (show_group_labels) {
    fixed_space <- fixed_space + 0.08
  }

  # Calculate total track height
  track_heights <- sapply(track_configs, function(x) x$height %||% 0.1)
  total_track_height <- sum(track_heights)

  # Available space (conservative estimate, leave 20% buffer)
  available_space <- 0.8 - fixed_space

  if (total_track_height > available_space) {
    warning(sprintf(
      "Track heights (%.2f) exceed available space (%.2f). Auto-adjusting...",
      total_track_height,
      available_space
    ))

    # Scale down proportionally
    scale_factor <- available_space / total_track_height * 0.95

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

#' Process Track Data with Scaling
#' @keywords internal
.process_track_data <- function(data, config) {
  track_data <- data[[config$column]]

  # First handle outliers (before scaling)
  if (isTRUE(config$cap_outliers)) {
    q_low <- quantile(track_data, 0.01, na.rm = TRUE)
    q_high <- quantile(track_data, 0.99, na.rm = TRUE)
    track_data <- pmin(pmax(track_data, q_low), q_high)
  }

  # Then apply scaling
  if (!is.null(config$scale)) {
    track_data <- switch(
      config$scale,
      zscore = as.numeric(scale(track_data)),
      minmax = scales::rescale(track_data, to = c(0, 1)),
      log2 = log2(pmax(track_data, 0) + 1),
      log10 = log10(pmax(track_data, 0) + 1),
      track_data
    )
  }

  return(track_data)
}

#' Extract Global Row Order After Clustering
#' @keywords internal
.extract_row_order <- function(data, group_col, clustered) {
  if (!clustered) {
    # Non-clustering mode: return simple order info without new_index
    row_order_list <- data %>%
      dplyr::group_by(!!sym(group_col)) %>%
      dplyr::mutate(
        sector = as.character(!!sym(group_col)),
        original_index = dplyr::row_number()
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(sector, original_index)

    return(row_order_list)
  } else {
    # Clustering mode: extract reordering info
    order_df <- data %>%
      dplyr::group_by(!!sym(group_col)) %>%
      dplyr::mutate(
        sector = as.character(!!sym(group_col)),
        original_index = dplyr::row_number()
      ) %>%
      dplyr::ungroup()

    # Try to get clustering order from circlize
    sector_names <- get.all.sector.index()

    for (sector in sector_names) {
      tryCatch(
        {
          cell_order <- get.cell.meta.data(
            "row_order",
            sector.index = sector,
            track.index = 1
          )

          sector_mask <- order_df$sector == sector
          order_df$new_index[sector_mask] <- cell_order
        },
        error = function(e) {
          warning(sprintf(
            "Could not extract row order for sector '%s'",
            sector
          ))
          # Fallback: use original order
          sector_mask <- order_df$sector == sector
          order_df$new_index[sector_mask] <- order_df$original_index[
            sector_mask
          ]
        }
      )
    }

    return(order_df)
  }
}

#' Calculate Feature Positions After Clustering
#' @keywords internal
.calculate_feature_positions <- function(
  data,
  group_col,
  feature_col,
  row_order_df
) {
  feature_mapping <- data %>%
    dplyr::mutate(
      original_index = dplyr::row_number(),
      sector = as.character(!!sym(group_col))
    ) %>%
    dplyr::select(
      feature = !!sym(feature_col),
      sector,
      original_index
    )

  # If clustering info available, use it
  if ("new_index" %in% colnames(row_order_df)) {
    feature_mapping <- feature_mapping %>%
      dplyr::left_join(
        row_order_df %>%
          dplyr::select(sector, original_index, new_index),
        by = c("sector", "original_index")
      ) %>%
      dplyr::group_by(sector) %>%
      dplyr::arrange(new_index) %>%
      dplyr::mutate(x_pos = dplyr::row_number() - 0.5) %>%
      dplyr::ungroup()
  } else {
    # No clustering: use original order
    feature_mapping <- feature_mapping %>%
      dplyr::group_by(sector) %>%
      dplyr::mutate(x_pos = dplyr::row_number() - 0.5) %>%
      dplyr::ungroup()
  }

  return(feature_mapping)
}

#' Draw Single Track
#' @keywords internal
.draw_track <- function(
  sector_data,
  config,
  group_col,
  n_groups,
  row_order_info,
  plot_params
) {
  current_sector <- CELL_META$sector.index

  # Check if new_index column exists
  if (!is.null(row_order_info) && "new_index" %in% colnames(row_order_info)) {
    # Clustering mode: reorder data
    sector_order_info <- row_order_info %>%
      dplyr::filter(sector == current_sector) %>%
      dplyr::arrange(new_index)

    col_data <- sector_data[[config$column]][sector_order_info$original_index]
  } else {
    # Non-clustering mode: use original order
    col_data <- sector_data[[config$column]]
  }

  n_points <- length(col_data)
  x_pos <- CELL_META$xlim[1] +
    (CELL_META$xlim[2] - CELL_META$xlim[1]) *
      (seq_len(n_points) - 0.5) /
      n_points

  point_cols <- if (is.function(config$color)) {
    config$color(col_data)
  } else {
    rep(config$color, length.out = n_points)
  }

  # Draw based on type
  switch(
    config$type,
    points = {
      circos.points(
        x_pos,
        col_data,
        pch = config$pch %||% 16,
        cex = config$cex %||% plot_params$track_point_cex,
        col = point_cols
      )
    },
    bars = {
      circos.barplot(
        col_data,
        pos = x_pos,
        col = point_cols,
        bar_width = config$bar_width %||% 0.6,
        border = config$border %||% NA
      )
    },
    lines = {
      circos.lines(
        x_pos,
        col_data,
        col = point_cols,
        lwd = config$lwd %||% plot_params$track_line_lwd
      )
    },
    area = {
      baseline <- config$baseline %||% 0
      circos.lines(
        x_pos,
        col_data,
        col = point_cols,
        area = TRUE,
        baseline = baseline,
        border = NA
      )
    },
    lollipop = {
      baseline <- config$baseline %||% 0
      circos.segments(
        x0 = x_pos,
        y0 = baseline,
        x1 = x_pos,
        y1 = col_data,
        col = point_cols,
        lwd = config$lwd %||% 1
      )
      circos.points(
        x_pos,
        col_data,
        pch = config$pch %||% 16,
        cex = config$cex %||% plot_params$track_point_cex,
        col = point_cols
      )
    },
    boxplot = {
      circos.boxplot(
        col_data,
        pos = CELL_META$xcenter,
        col = point_cols[1],
        border = "black",
        outline = FALSE
      )
    },
    violin = {
      circos.violin(
        col_data,
        pos = CELL_META$xcenter,
        col = point_cols[1],
        border = NA
      )
    },
    warning(sprintf("Unknown track type: %s", config$type))
  )

  # Add track label
  if (CELL_META$sector.numeric.index == n_groups) {
    val_range <- CELL_META$ylim
    mid_y <- mean(val_range)
    circos.lines(
      c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")),
      c(mid_y, mid_y),
      col = "black",
      lwd = 1
    )
    circos.text(
      CELL_META$cell.xlim[2] + convert_x(1.5, "mm"),
      mid_y,
      config$label,
      cex = plot_params$track_label_cex,
      adj = c(0, 0.5),
      facing = "inside"
    )
  }
}

#' Add Color Legend
#' @keywords internal
.add_heatmap_legend <- function(
  heatmap_col,
  value_range,
  title = "Value",
  x = 0.85,
  y = 0.85,
  width = 0.08,
  height = 0.5
) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    warning("ComplexHeatmap package required for legend")
    return(invisible(NULL))
  }

  # Create legend
  lgd <- ComplexHeatmap::Legend(
    col_fun = heatmap_col,
    title = title,
    at = seq(value_range[1], value_range[2], length.out = 5),
    labels_gp = grid::gpar(fontsize = 8),
    title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    grid_height = grid::unit(4, "mm"),
    grid_width = grid::unit(4, "mm")
  )

  # Draw legend
  ComplexHeatmap::draw(
    lgd,
    x = grid::unit(x, "npc"),
    y = grid::unit(y, "npc"),
    just = c("left", "top")
  )
}

# ==============================================================================
# Main Function
# ==============================================================================

#' Advanced Circular Heatmap with Enhanced Styling
#'
#' @param data Data frame with features (rows) and groups
#' @param group_col Column name for grouping variable
#' @param feature_col Column name for feature identifiers
#' @param value_cols Character vector of numeric columns for heatmap
#' @param group_palette Named color vector OR palette name (e.g., "npg", "lancet")
#' @param heatmap_col Color mapping function OR palette name for heatmap
#' @param track_configs List of track configurations
#' @param links_data Data frame with 'from', 'to', and optional 'color' columns
#' @param show_rownames Logical. Show feature names?
#' @param show_colnames Logical. Show column names?
#' @param gap_degree Numeric. Gap after last sector
#' @param gap_between Numeric. Gap between other sectors
#' @param start_degree Numeric. Starting angle
#' @param value_range Numeric vector. Range for data normalization
#' @param plot_params List of plotting parameters
#' @param group_label_params List of group label styling parameters
#' @param show_legend Logical. Show heatmap legend?
#' @param legend_params List of legend parameters
#' @param auto_adjust_heights Logical. Auto-adjust track heights?
#' @param ... Additional arguments passed to circos.heatmap()
#'
#' @export
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
      # 1.1 Check clustering
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
      # 1.2 Auto-adjust font size
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
      # 1.3 Validate track heights
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
      # 2. Data Preparation
      # ========================================================================
      color_result <- .prepare_group_colors(data, group_col, group_palette)
      data <- color_result$data
      group_colors <- color_result$colors

      heatmap_col <- .prepare_heatmap_colors(heatmap_col, value_range)
      heatmap_col <- local({
        raw_col <- heatmap_col
        function(x) {
          result <- raw_col(x)
          result[is.na(x)] <- "white"
          result
        }
      })

      data_normalized <- data %>%
        dplyr::group_by(!!sym(group_col)) %>%
        dplyr::mutate(dplyr::across(
          dplyr::all_of(value_cols),
          ~ if (all(is.na(.))) . else scales::rescale(., to = value_range)
        )) %>%
        dplyr::ungroup()

      data_matrix <- data_normalized %>%
        dplyr::select(dplyr::all_of(value_cols)) %>%
        as.matrix()
      rownames(data_matrix) <- data_normalized[[feature_col]]

      # ========================================================================
      # 3. Initialize Circos Canvas
      # ========================================================================
      circos.clear()

      group_levels <- levels(data_normalized[[group_col]])
      n_groups <- length(group_levels)

      gaps_vector <- rep(gap_between, n_groups)
      gaps_vector[n_groups] <- gap_degree

      track_margin <- if (n_tracks > 5 || will_cluster) {
        c(0.003, 0.003)
      } else if (n_tracks > 2) {
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

      # ========================================================================
      # 4. Draw Main Heatmap
      # ========================================================================
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

      # ========================================================================
      # 5. Extract Row Order
      # ========================================================================
      row_order_info <- .extract_row_order(
        data_normalized,
        group_col,
        clustered
      )

      # ========================================================================
      # 6. Add Column Names Track
      # ========================================================================
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

      # ========================================================================
      # 7. Add External Tracks
      # ========================================================================
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

      # ========================================================================
      # 8. Add Group Labels Track
      # ========================================================================
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

      # ========================================================================
      # 9. Add Links
      # ========================================================================
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

      # ========================================================================
      # 10. Add Legend
      # ========================================================================
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

      # ========================================================================
      # 11. Return Information
      # ========================================================================
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

# ==============================================================================
# Utility Function
# ==============================================================================

#' Preview Feature Positions
#' @param result Output from plot_advanced_circos()
#' @export
preview_feature_positions <- function(result) {
  if (is.null(result$feature_mapping)) {
    stop("No feature mapping available. Did you provide links_data?")
  }

  cat("Feature Positions:\n")
  cat(strrep("=", 80), "\n\n")

  if (
    !is.null(result$row_order_info) &&
      "new_index" %in% colnames(result$row_order_info)
  ) {
    cat("Note: Clustering was enabled. Positions reflect reordered layout.\n\n")
  }

  print(
    result$feature_mapping %>%
      dplyr::arrange(sector, x_pos) %>%
      dplyr::select(feature, sector, x_pos)
  )

  invisible(result$feature_mapping)
}

# ==============================================================================
# Print welcome message
# ==============================================================================
message("Advanced Circos Plot Functions Loaded Successfully!")
message("Functions available:")
message("  - plot_advanced_circos()")
message("  - get_palette_colors()")
message("  - list_color_palettes()")
message("  - preview_palette()")
message("  - preview_feature_positions()")
message("\nType ?plot_advanced_circos for help (if documentation installed)")
message("Or see the user manual for detailed instructions.")
