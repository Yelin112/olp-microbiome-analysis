library(ggplot2)
library(scales)
library(grDevices)
library(paletteer)

# ==============================================================================
# 1. Built-in Palette Library
# ==============================================================================
.builtin_palettes <- list(
  NPG = c(
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
  AAAS = c(
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
  RWB = c("#313695", "#FFFFFF", "#A50026"),
  Viridis = c("#440154", "#21908C", "#FDE725")
)

# ==============================================================================
# 2. Global Palette Manager (Singleton Pattern)
# ==============================================================================
.palette_manager <- new.env(parent = emptyenv())
.palette_manager$custom <- list()
.palette_manager$default <- "NPG"

#' Validate Color Vector
#'
#' @param colors Character vector of colors
#' @return Logical vector indicating valid colors
#' @keywords internal
.is_valid_color <- function(colors) {
  # Check hex colors: #RGB, #RGBA, #RRGGBB, #RRGGBBAA
  is_hex <- grepl("^#([0-9A-Fa-f]{3}){1,2}([0-9A-Fa-f]{2})?$", colors)

  # Check R named colors
  is_named <- colors %in% grDevices::colors()

  # Check valid rgb/rgba format (optional, for completeness)
  is_rgb <- grepl("^rgba?\\(", colors)

  return(is_hex | is_named | is_rgb)
}

#' 存入自定义色板 (支持单个或批量列表导入)
#'
#' @param name_or_list 字符型(单个画板名) 或 命名列表(批量画板)。
#' @param colors 字符型向量。颜色值。当使用批量列表导入时，此参数可省略。
#' @param overwrite 逻辑值。是否允许覆盖已存在的同名色板？
#' @export
set_palette <- function(name_or_list, colors = NULL, overwrite = FALSE) {
  # ==========================================
  # 模式 1：批量导入模式 (传入的是 List)
  # ==========================================
  if (is.list(name_or_list)) {
    pal_list <- name_or_list

    # 校验：List 必须有完整的名字
    if (is.null(names(pal_list)) || any(names(pal_list) == "")) {
      stop(
        "批量导入时，List 必须为每个元素命名 (例如: list('PaletteA' = c(...)))"
      )
    }

    # 遍历 List，递归调用自身逐个注册
    for (p_name in names(pal_list)) {
      set_palette(
        name_or_list = p_name,
        colors = pal_list[[p_name]],
        overwrite = overwrite
      )
    }

    message(sprintf("🎉 成功批量导入 %d 个画板！", length(pal_list)))
    return(invisible(pal_list))
  }

  # ==========================================
  # 模式 2：单个导入模式 (传入的是 字符串)
  # ==========================================
  name <- name_or_list

  # 输入合法性校验
  if (!is.character(name) || length(name) != 1 || nchar(name) == 0) {
    stop("'name_or_list' 必须是单个非空字符串或一个命名列表(List)")
  }
  if (!is.character(colors) || length(colors) == 0) {
    stop("'colors' 必须是非空颜色向量")
  }

  # 校验 Hex 格式或 R 原生颜色名 (.is_valid_color 为之前定义的内部函数)
  invalid_colors <- !.is_valid_color(colors)
  if (any(invalid_colors)) {
    invalid_idx <- which(invalid_colors)
    warning(sprintf(
      "画板 '%s' 中存在可能无效的颜色: %s\n有效格式: Hex (#FF0000), 或原生名 (red)",
      name,
      paste(colors[invalid_idx], collapse = ", ")
    ))
  }

  # 查重与覆盖逻辑
  if (name %in% names(.palette_manager$custom) && !overwrite) {
    stop(sprintf("画板 '%s' 已存在。如需替换请设置 overwrite = TRUE", name))
  }

  # 存入全局环境
  .palette_manager$custom[[name]] <- colors
  message(sprintf("✓ 画板 '%s' 已保存 (包含 %d 种颜色)", name, length(colors)))
  invisible(colors)
}

#' List Available Palettes
#'
#' @param type Character. "all", "builtin", or "custom"
#' @return Character vector of palette names
#' @export
list_palettes <- function(type = c("all", "builtin", "custom")) {
  type <- match.arg(type)

  builtin <- names(.builtin_palettes)
  custom <- names(.palette_manager$custom)

  switch(type, all = c(builtin, custom), builtin = builtin, custom = custom)
}

#' Remove Custom Palette
#'
#' @param name Character. Palette name to remove
#' @export
remove_palette <- function(name) {
  if (!name %in% names(.palette_manager$custom)) {
    warning(sprintf("Palette '%s' not found", name))
    return(invisible(FALSE))
  }
  .palette_manager$custom[[name]] <- NULL
  message(sprintf("✓ Palette '%s' removed", name))
  invisible(TRUE)
}

#' Preview Color Palette
#'
#' @param palette Character or vector. Palette name or color vector
#' @param n Integer. Number of colors to show
#' @export
preview_palette <- function(palette = "NPG", n = NULL) {
  cols <- get_colors(palette, n = n, type = "discrete")

  n_colors <- length(cols)
  labels <- if (!is.null(names(cols))) names(cols) else seq_along(cols)

  # Create preview plot
  df <- data.frame(
    x = seq_along(cols),
    y = 1,
    label = labels,
    color = cols
  )

  p <- ggplot(df, aes(x = x, y = y, fill = color)) +
    geom_tile(color = "white", linewidth = 2) +
    geom_text(aes(label = label), size = 3.5, fontface = "bold") +
    scale_fill_identity() +
    labs(
      title = sprintf(
        "Palette: %s (%d colors)",
        if (is.character(palette) && length(palette) == 1) {
          palette
        } else {
          "Custom"
        },
        n_colors
      )
    ) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

  print(p)
  invisible(cols)
}
# ==============================================================================
# 3. 万能颜色提取引擎 (Universal Color Extraction Engine)
# ==============================================================================
#' 提取颜色引擎
#'
#' @param palette 字符型或向量。色板名称或颜色向量
#' @param n 整型。需要的颜色数量
#' @param type 字符型。"discrete"(离散) 或 "continuous"(连续)
#' @param warn 逻辑值。是否显示插值警告？
#' @return 颜色字符向量
#' @export
get_colors <- function(
  palette = "NPG",
  n = NULL,
  type = c("discrete", "continuous"),
  warn = TRUE
) {
  type <- match.arg(type)

  # 1. 如果用户直接传入了颜色向量
  if (length(palette) > 1 || (!is.character(palette))) {
    cols <- as.character(palette)
  } else {
    # 按照优先级搜索: 自定义 -> 内置 -> RColorBrewer -> paletteer
    cols <- NULL

    # 优先检测自定义字典
    if (palette %in% names(.palette_manager$custom)) {
      cols <- .palette_manager$custom[[palette]]
    } else if (palette %in% names(.builtin_palettes)) {
      # 检测内置
      cols <- .builtin_palettes[[palette]]
    } else if (
      requireNamespace("RColorBrewer", quietly = TRUE) &&
        palette %in% rownames(RColorBrewer::brewer.pal.info)
    ) {
      # 检测 RColorBrewer
      max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
      cols <- RColorBrewer::brewer.pal(max_n, palette)
    } else if (
      requireNamespace("paletteer", quietly = TRUE) &&
        grepl("::", palette)
    ) {
      # 检测外部生态 Paletteer (格式: "包名::色板名")
      # 【核心升级】：双重捕获机制，彻底打通离散与连续画板
      cols <- tryCatch(
        {
          # 尝试一：离散型画板 (Discrete)
          as.character(paletteer::paletteer_d(palette))
        },
        error = function(e1) {
          # 尝试二：如果离散型报错，尝试作为连续型画板提取 (Continuous)
          tryCatch(
            {
              # 对于连续色板，默认提取 100 个颜色作为基准过渡池
              as.character(paletteer::paletteer_c(palette, n = 100))
            },
            error = function(e2) {
              # 如果两边都找不到，再弹出警告并返回 NULL
              warning(sprintf(
                "无法从 paletteer 加载色板 '%s', 请检查包名和色板名是否拼写正确",
                palette
              ))
              return(NULL)
            }
          )
        }
      )
    }

    # 终极回退机制
    if (is.null(cols)) {
      if (warn) {
        warning(sprintf(
          "未找到色板 '%s', 已回退至默认色板 '%s'",
          palette,
          .palette_manager$default
        ))
      }
      cols <- .builtin_palettes[[.palette_manager$default]]
    }
  }

  # 2. 如果不需要特定数量，直接返回全组颜色
  if (is.null(n)) {
    return(cols)
  }

  # 3. 严格映射模式 (带有 names 的离散变量)
  if (!is.null(names(cols)) && type == "discrete") {
    if (n != length(cols) && warn) {
      warning(sprintf(
        "命名色板包含 %d 种颜色，但请求了 %d 种。为保证映射严格性，将返回全组颜色。",
        length(cols),
        n
      ))
    }
    return(cols)
  }

  # 清除 names 以便灵活取用
  cols <- unname(cols)

  # 4. 连续型或颜色不足时的插值处理
  if (type == "continuous" || n > length(cols)) {
    # 只有在离散型且颜色不够时，才提示插值警告
    if (n > length(cols) && type == "discrete" && warn) {
      message(sprintf(
        "ℹ 颜色不足，正在将色板从 %d 色平滑插值至 %d 色",
        length(cols),
        n
      ))
    }
    return(colorRampPalette(cols)(n))
  }

  # 颜色充足时的离散型截取
  return(cols[seq_len(n)])
}

# ==============================================================================
# 4. ggplot2 Scale Functions
# ==============================================================================
#' Universal Discrete Fill Scale
#'
#' @param palette Character. Palette name
#' @param alpha Numeric. Transparency (0-1)
#' @param ... Additional arguments passed to scale_*_manual or discrete_scale
#' @export
scale_fill_univ_d <- function(palette = "NPG", alpha = 1, ...) {
  cols <- get_colors(palette, type = "discrete", warn = FALSE)

  # Named vector = strict mapping mode
  if (!is.null(names(cols))) {
    if (alpha < 1) {
      # 【修复】：显式使用 scales::alpha
      cols <- setNames(scales::alpha(cols, alpha), names(cols))
    }
    return(scale_fill_manual(values = cols, ...))
  }

  # Unnamed vector = flexible palette mode
  pal_func <- function(n) {
    # 【修复】：显式使用 scales::alpha
    scales::alpha(
      get_colors(palette, n = n, type = "discrete", warn = TRUE),
      alpha
    )
  }
  discrete_scale("fill", "univ_d", palette = pal_func, ...)
}

#' Universal Discrete Color Scale
#' @rdname scale_fill_univ_d
#' @export
scale_color_univ_d <- function(palette = "NPG", alpha = 1, ...) {
  cols <- get_colors(palette, type = "discrete", warn = FALSE)

  if (!is.null(names(cols))) {
    if (alpha < 1) {
      # 【修复】：显式使用 scales::alpha
      cols <- setNames(scales::alpha(cols, alpha), names(cols))
    }
    return(scale_color_manual(values = cols, ...))
  }

  pal_func <- function(n) {
    # 【修复】：显式使用 scales::alpha
    scales::alpha(
      get_colors(palette, n = n, type = "discrete", warn = TRUE),
      alpha
    )
  }
  discrete_scale("colour", "univ_d", palette = pal_func, ...)
}

#' Universal Continuous Fill Scale
#'
#' @param palette Character. Palette name
#' @param alpha Numeric. Transparency (0-1)
#' @param reverse Logical. Reverse color order?
#' @param ... Additional arguments passed to scale_fill_gradientn
#' @export
scale_fill_univ_c <- function(
  palette = "RWB",
  alpha = 1,
  reverse = FALSE,
  ...
) {
  cols <- get_colors(palette, type = "continuous")
  if (reverse) {
    cols <- rev(cols)
  }

  # 【修复】：正确应用 alpha 参数
  if (alpha < 1) {
    cols <- scales::alpha(cols, alpha)
  }

  scale_fill_gradientn(colours = cols, ...)
}

#' Universal Continuous Color Scale
#' @rdname scale_fill_univ_c
#' @export
scale_color_univ_c <- function(
  palette = "RWB",
  alpha = 1,
  reverse = FALSE,
  ...
) {
  cols <- get_colors(palette, type = "continuous")
  if (reverse) {
    cols <- rev(cols)
  }

  # 【修复】：正确应用 alpha 参数
  if (alpha < 1) {
    cols <- scales::alpha(cols, alpha)
  }

  scale_color_gradientn(colours = cols, ...)
}

# ==============================================================================
# 5. Utility Functions
# ==============================================================================
#' Set Default Palette
#'
#' @param palette Character. Default palette name
#' @export
set_default_palette <- function(palette) {
  if (!palette %in% list_palettes("all")) {
    stop(sprintf("Palette '%s' not found", palette))
  }
  .palette_manager$default <- palette
  message(sprintf("✓ Default palette set to '%s'", palette))
}

#' Get Current Default Palette
#' @export
get_default_palette <- function() {
  .palette_manager$default
}


# 1. 保存脚本
# 将上面的代码保存为 "plot_advanced_circos.R"

# 2. 在 R 中加载

data <- data.frame(
  Gene = paste0("Gene_", 1:50),
  CellType = rep(c("T_cell", "B_cell"), each = 25),
  Sample_1 = runif(50, 0, 100),
  Sample_2 = runif(50, 0, 100),
  Sample_3 = runif(50, 0, 100)
)

plot_advanced_circos(
  data = data,
  group_col = "CellType",
  value_cols = c("Sample_1", "Sample_2", "Sample_3")
)


# ==============================================================================
# 综合测试函数：复杂数据 + 所有轨道类型
# ==============================================================================

#' 生成复杂的多组学测试数据
#' @export
generate_complex_test_data <- function() {
  set.seed(42)

  # 定义细胞类型
  cell_types <- c(
    "CD4_T",
    "CD8_T",
    "B_cell",
    "NK_cell",
    "Monocyte",
    "Macrophage"
  )
  n_celltypes <- length(cell_types)
  n_genes_per_type <- 30 # 每个细胞类型30个特征基因
  n_samples <- 8 # 8个样本

  n_total <- n_celltypes * n_genes_per_type

  # 基础数据框
  data <- data.frame(
    Gene = paste0("Gene_", sprintf("%03d", 1:n_total)),
    CellType = rep(cell_types, each = n_genes_per_type),
    stringsAsFactors = FALSE
  )

  # 生成表达矩阵（模拟单细胞RNA-seq数据）
  for (i in 1:n_samples) {
    col_name <- sprintf("Sample_%02d", i)

    # 每个细胞类型有不同的表达模式
    expr <- numeric(n_total)
    for (j in 1:n_celltypes) {
      idx_start <- (j - 1) * n_genes_per_type + 1
      idx_end <- j * n_genes_per_type

      # 细胞类型特异性基因表达
      base_expr <- rnorm(n_genes_per_type, mean = 50 + j * 10, sd = 20)
      # 添加样本间变异
      sample_variation <- rnorm(n_genes_per_type, mean = 0, sd = 5)

      expr[idx_start:idx_end] <- pmax(base_expr + sample_variation, 0)
    }

    data[[col_name]] <- expr
  }

  # ========================================================================
  # 添加多个额外的数据维度（用于不同类型的轨道）
  # ========================================================================

  # 1. 差异表达分析结果
  data$LogFC <- rnorm(n_total, mean = 0, sd = 2.5)
  data$LogFC[sample(n_total, 20)] <- rnorm(20, mean = 5, sd = 1) # 上调基因
  data$LogFC[sample(n_total, 20)] <- rnorm(20, mean = -5, sd = 1) # 下调基因

  # 2. 统计显著性
  data$Pvalue <- runif(n_total, 0.0001, 0.5)
  data$Pvalue[abs(data$LogFC) > 2] <- runif(
    sum(abs(data$LogFC) > 2),
    0.0001,
    0.01
  )
  data$NegLog10P <- -log10(data$Pvalue)

  # 3. 基因重要性评分（随机森林等）
  data$Importance <- abs(rnorm(n_total, mean = 50, sd = 20))
  data$Importance[abs(data$LogFC) > 3] <- data$Importance[abs(data$LogFC) > 3] +
    30

  # 4. 保守性分数（进化保守性）
  data$Conservation <- rbeta(n_total, 2, 5) * 100

  # 5. 染色质可及性（ATAC-seq）
  data$Chromatin_Access <- abs(rnorm(n_total, mean = 40, sd = 15))

  # 6. 甲基化水平
  data$Methylation <- rbeta(n_total, 2, 2) * 100

  # 7. 拷贝数变异（CNV）
  data$CNV_Score <- rnorm(n_total, mean = 0, sd = 0.3)
  data$CNV_Score[sample(n_total, 15)] <- rnorm(15, mean = 1.5, sd = 0.5) # 扩增
  data$CNV_Score[sample(n_total, 15)] <- rnorm(15, mean = -1.5, sd = 0.5) # 缺失

  # 8. 蛋白质丰度
  data$Protein_Abundance <- abs(rnorm(n_total, mean = 30, sd = 12))

  # 9. 细胞周期评分
  data$CellCycle_Score <- runif(n_total, 0, 100)

  # 10. 代谢活性
  data$Metabolic_Activity <- abs(rnorm(n_total, mean = 45, sd = 18))

  # 返回结果
  list(
    data = data,
    cell_types = cell_types,
    sample_cols = grep("^Sample_", colnames(data), value = TRUE),
    n_genes = n_total,
    n_samples = n_samples
  )
}

#' 生成基因-基因相互作用网络
#' @export
generate_gene_interactions <- function(data, n_links = 40) {
  set.seed(123)

  genes <- unique(data$Gene)
  cell_types <- unique(data$CellType)

  links <- data.frame(
    from = character(n_links),
    to = character(n_links),
    color = character(n_links),
    stringsAsFactors = FALSE
  )

  for (i in 1:n_links) {
    # 随机选择两个基因
    from_gene <- sample(genes, 1)
    to_gene <- sample(genes, 1)

    # 确保不是同一个基因
    while (from_gene == to_gene) {
      to_gene <- sample(genes, 1)
    }

    # 获取基因所在的细胞类型
    from_type <- data$CellType[data$Gene == from_gene][1]
    to_type <- data$CellType[data$Gene == to_gene][1]

    # 根据是否跨细胞类型设置颜色
    if (from_type == to_type) {
      # 同一细胞类型内：蓝色
      link_color <- "#4DBBD580"
    } else {
      # 跨细胞类型：红色
      link_color <- "#DC000080"
    }

    links$from[i] <- from_gene
    links$to[i] <- to_gene
    links$color[i] <- link_color
  }

  return(links)
}

#' 完整测试：所有轨道类型 + 复杂数据
#' @export
test_all_track_types <- function() {
  cat(strrep("=", 80), "\n")
  cat("COMPREHENSIVE TEST: All Track Types + Complex Data\n")
  cat(strrep("=", 80), "\n\n")

  # ========================================================================
  # 1. 生成数据
  # ========================================================================
  cat("Step 1: Generating complex multi-omics data...\n")
  test_data <- generate_complex_test_data()
  data <- test_data$data

  cat(sprintf("  - Features: %d genes\n", test_data$n_genes))
  cat(sprintf("  - Groups: %d cell types\n", length(test_data$cell_types)))
  cat(sprintf("  - Samples: %d\n", test_data$n_samples))
  cat(sprintf("  - Additional dimensions: 10\n\n"))

  # ========================================================================
  # 2. 生成连线
  # ========================================================================
  cat("Step 2: Generating gene-gene interaction network...\n")
  links <- generate_gene_interactions(data, n_links = 40)
  cat(sprintf("  - Total links: %d\n", nrow(links)))
  cat(sprintf(
    "  - Intra-type links: %d (blue)\n",
    sum(grepl("#4DBBD5", links$color))
  ))
  cat(sprintf(
    "  - Inter-type links: %d (red)\n\n",
    sum(grepl("#DC0000", links$color))
  ))

  # ========================================================================
  # 3. 配置轨道 - 测试所有7种类型
  # ========================================================================
  cat("Step 3: Configuring tracks (7 types)...\n")

  track_configs <- list(
    # 轨道1: BARS - 差异倍数
    list(
      column = "LogFC",
      type = "bars",
      height = 0.09,
      label = "Log2FC",
      color = function(x) {
        ifelse(
          x > 2,
          "#DC0000",
          ifelse(x > 0, "#F39B7F", ifelse(x > -2, "#3C5488", "#2166AC"))
        )
      },
      scale = NULL,
      cap_outliers = TRUE,
      bar_width = 0.8
    ),

    # 轨道2: POINTS - 显著性
    list(
      column = "NegLog10P",
      type = "points",
      height = 0.07,
      label = "-log10P",
      color = colorRamp2(
        c(0, 1, 2, 3, 4),
        c("#CCCCCC", "#FFA500", "#FF6347", "#DC143C", "#8B0000")
      ),
      scale = NULL,
      pch = 16,
      cex = 0.9
    ),

    # 轨道3: LINES - 染色质可及性
    list(
      column = "Chromatin_Access",
      type = "lines",
      height = 0.06,
      label = "ATAC",
      color = "#00A087",
      scale = "minmax",
      lwd = 2
    ),

    # 轨道4: AREA - 基因重要性
    list(
      column = "Importance",
      type = "area",
      height = 0.06,
      label = "Importance",
      color = "#9467BD60",
      scale = "minmax",
      baseline = 0
    ),

    # 轨道5: LOLLIPOP - 甲基化水平
    list(
      column = "Methylation",
      type = "lollipop",
      height = 0.07,
      label = "Methyl",
      color = function(x) ifelse(x > 60, "#E64B35", "#4DBBD5"),
      scale = NULL,
      lwd = 1.2,
      pch = 16,
      cex = 0.8
    ),

    # 轨道6: BOXPLOT - CNV分布（每组一个箱子）
    list(
      column = "CNV_Score",
      type = "boxplot",
      height = 0.08,
      label = "CNV",
      color = "#8491B4",
      scale = NULL
    ),

    # 轨道7: VIOLIN - 蛋白质丰度分布
    list(
      column = "Protein_Abundance",
      type = "violin",
      height = 0.08,
      label = "Protein",
      color = "#F39B7F",
      scale = "zscore"
    )
  )

  for (i in seq_along(track_configs)) {
    cat(sprintf(
      "  Track %d: %-10s - %s (height=%.2f)\n",
      i,
      track_configs[[i]]$type,
      track_configs[[i]]$label,
      track_configs[[i]]$height
    ))
  }
  cat("\n")

  # ========================================================================
  # 4. 自定义配色
  # ========================================================================
  cat("Step 4: Setting up color schemes...\n")

  celltype_colors <- c(
    "CD4_T" = "#E64B35",
    "CD8_T" = "#4DBBD5",
    "B_cell" = "#00A087",
    "NK_cell" = "#3C5488",
    "Monocyte" = "#F39B7F",
    "Macrophage" = "#8491B4"
  )

  cat("  - Group palette: NPG-inspired (6 colors)\n")
  cat("  - Heatmap palette: Red-Blue diverging\n\n")

  # ========================================================================
  # 5. 绘图
  # ========================================================================
  cat("Step 5: Generating circular heatmap...\n\n")

  result <- plot_advanced_circos(
    # 数据
    data = data,
    group_col = "CellType",
    feature_col = "Gene",
    value_cols = test_data$sample_cols,

    # 配色
    group_palette = celltype_colors,
    heatmap_col = colorRamp2(
      seq(0, 0.5, length.out = 11),
      rev(RColorBrewer::brewer.pal(11, "RdBu"))
    ),

    # 轨道和连线
    track_configs = track_configs,
    links_data = links,

    # 显示控制
    show_rownames = FALSE, # 基因太多，不显示名称
    show_colnames = TRUE,
    show_legend = FALSE, # 7个轨道，空间紧张

    # 布局
    gap_degree = 40,
    gap_between = 2,
    start_degree = 90,
    value_range = c(0, 0.5),

    # 分组标签样式（轨道多，字体自动缩小）
    group_label_params = list(
      cex = NULL, # 自动调整
      font = 2,
      col = "white",
      bg_alpha = 0.6
    ),

    # 绘图参数
    plot_params = list(
      heatmap_track_height = 0.16, # 压缩热图为轨道腾出空间
      heatmap_rownames_cex = 0.4,
      track_label_cex = 0.85,
      group_label_cex = 1.0
    ),

    # 不聚类（保持原始顺序）
    cluster = FALSE
  )

  # ========================================================================
  # 6. 输出统计信息
  # ========================================================================
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("TEST COMPLETED SUCCESSFULLY!\n")
  cat(strrep("=", 80), "\n\n")

  cat("Summary:\n")
  cat(sprintf("  - Clustered: %s\n", result$clustered))
  cat(sprintf("  - Total tracks: %d\n", length(track_configs)))
  cat(sprintf(
    "  - Track types: %s\n",
    paste(unique(sapply(track_configs, function(x) x$type)), collapse = ", ")
  ))
  cat(sprintf("  - Links drawn: %d\n", nrow(links)))
  cat(sprintf(
    "  - Total height: %.2f\n",
    sum(sapply(track_configs, function(x) x$height)) + 0.16 + 0.08
  ))

  cat("\nColors used:\n")
  for (ct in names(celltype_colors)) {
    cat(sprintf("  - %-12s: %s\n", ct, celltype_colors[ct]))
  }

  return(invisible(result))
}

#' 测试：聚类模式 + 所有轨道类型
#' @export
test_all_tracks_with_clustering <- function() {
  cat(strrep("=", 80), "\n")
  cat("ADVANCED TEST: All Track Types + Clustering\n")
  cat(strrep("=", 80), "\n\n")

  # 生成数据
  test_data <- generate_complex_test_data()
  data <- test_data$data
  links <- generate_gene_interactions(data, n_links = 30)

  # 轨道配置（数量减少以适应聚类）
  track_configs <- list(
    list(
      column = "LogFC",
      type = "bars",
      height = 0.08,
      label = "Log2FC",
      color = function(x) ifelse(x > 0, "#DC0000", "#3C5488"),
      cap_outliers = TRUE
    ),
    list(
      column = "NegLog10P",
      type = "lollipop",
      height = 0.07,
      label = "-log10P",
      color = colorRamp2(c(0, 2, 4), c("grey", "orange", "red")),
      lwd = 1.2
    ),
    list(
      column = "Chromatin_Access",
      type = "lines",
      height = 0.06,
      label = "ATAC",
      color = "#00A087",
      scale = "minmax",
      lwd = 2
    ),
    list(
      column = "Importance",
      type = "area",
      height = 0.06,
      label = "Importance",
      color = "#9467BD60",
      scale = "minmax"
    ),
    list(
      column = "Protein_Abundance",
      type = "points",
      height = 0.06,
      label = "Protein",
      color = "#F39B7F",
      scale = "zscore",
      pch = 16
    )
  )

  cat("Configuration:\n")
  cat(sprintf("  - Tracks: %d\n", length(track_configs)))
  cat("  - Clustering: Enabled\n")
  cat("  - Dendrogram: Inside\n\n")

  result <- plot_advanced_circos(
    data = data,
    group_col = "CellType",
    feature_col = "Gene",
    value_cols = test_data$sample_cols,
    group_palette = "npg",
    heatmap_col = "rdbu",
    track_configs = track_configs,
    links_data = links,
    show_rownames = FALSE,
    show_colnames = TRUE,
    gap_degree = 45,

    # 启用聚类
    cluster = TRUE,
    dend.side = "inside",
    dend.track.height = 0.12,

    plot_params = list(
      heatmap_track_height = 0.14 # 为聚类树和轨道腾出空间
    )
  )

  cat("\nClustering results:\n")
  cat(sprintf("  - Clustered: %s\n", result$clustered))
  cat("  - Links automatically synchronized to new order\n")
  cat("  - Tracks automatically synchronized to new order\n\n")

  return(invisible(result))
}

#' 测试：极限情况（最多轨道）
#' @export
test_maximum_tracks <- function() {
  cat(strrep("=", 80), "\n")
  cat("EXTREME TEST: Maximum Number of Tracks\n")
  cat(strrep("=", 80), "\n\n")

  # 生成数据
  test_data <- generate_complex_test_data()
  data <- test_data$data

  # 10个轨道！
  track_configs <- list(
    list(
      column = "LogFC",
      type = "bars",
      height = 0.04,
      label = "FC",
      color = function(x) ifelse(x > 0, "red", "blue")
    ),
    list(
      column = "NegLog10P",
      type = "points",
      height = 0.04,
      label = "P",
      color = "#DC0000",
      pch = 16
    ),
    list(
      column = "Importance",
      type = "lines",
      height = 0.04,
      label = "Imp",
      color = "#4DBBD5",
      lwd = 1.5
    ),
    list(
      column = "Conservation",
      type = "area",
      height = 0.04,
      label = "Cons",
      color = "#00A08750",
      scale = "minmax"
    ),
    list(
      column = "Chromatin_Access",
      type = "lollipop",
      height = 0.04,
      label = "ATAC",
      color = "#3C5488",
      lwd = 1
    ),
    list(
      column = "Methylation",
      type = "bars",
      height = 0.04,
      label = "Meth",
      color = "#F39B7F"
    ),
    list(
      column = "CNV_Score",
      type = "points",
      height = 0.04,
      label = "CNV",
      color = "#8491B4",
      pch = 16
    ),
    list(
      column = "Protein_Abundance",
      type = "lines",
      height = 0.04,
      label = "Prot",
      color = "#E64B35",
      lwd = 1.5
    ),
    list(
      column = "CellCycle_Score",
      type = "area",
      height = 0.04,
      label = "Cycle",
      color = "#91D1C250",
      scale = "minmax"
    ),
    list(
      column = "Metabolic_Activity",
      type = "lollipop",
      height = 0.04,
      label = "Metab",
      color = "#DC0000",
      lwd = 1
    )
  )

  cat(sprintf("Attempting to draw %d tracks...\n\n", length(track_configs)))

  result <- plot_advanced_circos(
    data = data,
    group_col = "CellType",
    feature_col = "Gene",
    value_cols = test_data$sample_cols,
    group_palette = "jco",
    heatmap_col = "rdbu",
    track_configs = track_configs,
    show_rownames = FALSE,
    show_colnames = FALSE, # 关闭列名节省空间
    gap_degree = 30,

    group_label_params = list(
      cex = 0.7, # 小字体
      track_height = 0.06
    ),

    plot_params = list(
      heatmap_track_height = 0.12, # 极度压缩
      track_label_cex = 0.7
    ),

    auto_adjust_heights = TRUE # 自动调整
  )

  cat("\nExtreme test completed!\n")
  cat(sprintf("  - Successfully drew %d tracks\n", length(track_configs)))

  return(invisible(result))
}

# ==============================================================================
# 运行所有测试
# ==============================================================================

#' 运行所有测试并保存结果
#' @export
run_all_comprehensive_tests <- function(
  save_plots = TRUE,
  output_dir = "circos_tests"
) {
  if (save_plots) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    cat(sprintf("Plots will be saved to: %s\n\n", output_dir))
  }

  tests <- list(
    list(name = "all_track_types", func = test_all_track_types),
    list(
      name = "all_tracks_clustering",
      func = test_all_tracks_with_clustering
    ),
    list(name = "maximum_tracks", func = test_maximum_tracks)
  )

  results <- list()

  for (i in seq_along(tests)) {
    test <- tests[[i]]
    cat("\n")
    cat(strrep("#", 80), "\n")
    cat(sprintf("Running Test %d/%d: %s\n", i, length(tests), test$name))
    cat(strrep("#", 80), "\n\n")

    if (save_plots) {
      pdf(
        file.path(output_dir, sprintf("%s.pdf", test$name)),
        width = 14,
        height = 14
      )
    }

    results[[test$name]] <- test$func()

    if (save_plots) {
      dev.off()
      cat(sprintf("\n✓ Saved: %s.pdf\n", test$name))
    }

    cat("\n")
  }

  cat(strrep("=", 80), "\n")
  cat("ALL TESTS COMPLETED!\n")
  cat(strrep("=", 80), "\n\n")

  if (save_plots) {
    cat(sprintf("All plots saved to: %s/\n", output_dir))
  }

  return(invisible(results))
}


# 方式1: 运行单个测试
test_all_track_types() # 7种轨道类型
test_all_tracks_with_clustering() # 聚类 + 5种轨道
test_maximum_tracks() # 极限测试：10个轨道

# 方式2: 运行所有测试并保存PDF
results <- run_all_comprehensive_tests(save_plots = TRUE)

# 方式3: 只运行测试不保存
results <- run_all_comprehensive_tests(save_plots = FALSE)

install.packages("pulsar")
install.packages(
  "D:\\Downloads\\pulsar_0.3.11.tar.gz",
  repos = NULL,
  type = "source"
)
remotes::install_github("zdk123/SpiecEasi")
devtools::install_github("Jiawang1209/ggNetView")


# 网络分析调试 -----------------------------------------------------------------

source('Scripts\\tidyplot_themes.R')

source('调试用.r')


library(ggplot2)

# ============================================================================
# 准备：加载主题管理脚本（如果要使用自定义主题）

# ============================================================================
# 准备示例数据
# ============================================================================
set.seed(123)
abundance <- data.frame(
  Bacteroides = c(rnorm(15, 0.3, 0.05), rnorm(15, 0.2, 0.05)),
  Firmicutes = c(rnorm(15, 0.4, 0.08), rnorm(15, 0.5, 0.08)),
  Proteobacteria = c(rnorm(15, 0.1, 0.03), rnorm(15, 0.15, 0.03))
)
groups <- factor(rep(c("健康", "IBD"), each = 15))

# ============================================================================
# 示例1：基础用法（默认配色）
# ============================================================================
result1 <- compare_species_abundance_v2(
  abundance_matrix = abundance,
  grouping_factor = groups,
  species_name = "Bacteroides",
  plot_type = "boxplot",
  add_points = TRUE,
  point_type = "jitter"
)
result1$plot

# ============================================================================
# 示例2：自定义配色（Nature期刊风格）
# ============================================================================
result2 <- compare_species_abundance_v2(
  abundance_matrix = abundance,
  grouping_factor = groups,
  species_name = "Bacteroides",
  plot_title = "拟杆菌属丰度分析",
  plot_type = "boxplot",
  filled = FALSE,
  border_width = 1.5,
  add_points = TRUE,
  point_type = "jitter",
  color = c("#E64B35", "#4DBBD5") # Nature配色
)
result2$plot

# ============================================================================
# 示例3：自定义配色（JAMA医学期刊风格）
# ============================================================================
result3 <- compare_species_abundance_v2(
  abundance_matrix = abundance,
  grouping_factor = groups,
  species_name = "Firmicutes",
  plot_type = "violin",
  filled = TRUE,
  main_alpha = 0.6,
  add_points = TRUE,
  point_type = "beeswarm",
  color_palette = c("#374E55", "#DF8F44") # JAMA配色
)
result3$plot

# ============================================================================
# 示例4：自定义配色（鲜艳对比色）
# ============================================================================
result4 <- compare_species_abundance_v2(
  abundance_matrix = abundance,
  grouping_factor = groups,
  species_name = "Proteobacteria",
  plot_type = "barplot",
  filled = TRUE,
  add_points = TRUE,
  point_type = "jitter",
  color_palette = c("#FF6B6B", "#4ECDC4") # 鲜艳配色
)
result4$plot

# ============================================================================
# 示例5：色盲友好配色
# ============================================================================
result5 <- compare_species_abundance_v2(
  abundance_matrix = abundance,
  grouping_factor = groups,
  species_name = "Bacteroides",
  plot_type = "violin",
  filled = TRUE,
  add_points = TRUE,
  point_type = "jitter",
  point_shape = 21,
  point_color = "white",
  color_palette = c("#E69F00", "#56B4E9") # 色盲友好
)
result5$plot

# ============================================================================
# 示例6：灰度配色（适合黑白打印）
# ============================================================================
result6 <- compare_species_abundance_v2(
  abundance_matrix = abundance,
  grouping_factor = groups,
  species_name = "Firmicutes",
  plot_type = "boxplot",
  filled = TRUE,
  add_points = TRUE,
  point_type = "beeswarm",
  color_palette = c("#404040", "#BFBFBF") # 灰度配色
)
result6$plot

# ============================================================================
# 示例7：不填充的箱线图 + 深色边框 + 自定义配色
# ============================================================================
result7 <- compare_species_abundance_v2(
  abundance_matrix = abundance,
  grouping_factor = groups,
  species_name = "Proteobacteria",
  plot_type = "boxplot",
  filled = FALSE,
  border_width = 2,
  add_points = TRUE,
  point_type = "beeswarm",
  point_size = 3,
  point_color = "#2C3E50",
  color_palette = c("#3498DB", "#E74C3C") # 蓝红配色
)
result7$plot

# ============================================================================
# 示例8：多组比较 + 自定义配色
# ============================================================================
multi_groups <- factor(rep(c("A", "B", "C"), each = 10))
abundance_multi <- data.frame(
  species1 = c(rnorm(10, 0.2, 0.05), rnorm(10, 0.3, 0.05), rnorm(10, 0.4, 0.05))
)

result8 <- compare_species_abundance_v2(
  abundance_matrix = abundance_multi,
  grouping_factor = multi_groups,
  species_name = "species1",
  plot_type = "violin",
  filled = TRUE,
  add_points = TRUE,
  point_type = "jitter",
  color_palette = c("#E64B35", "#4DBBD5", "#00A087") # 三色配色
)
result8$plot

# ============================================================================
# 示例9：批量分析 + 统一配色
# ============================================================================
batch_results <- batch_compare_species(
  abundance_matrix = abundance,
  grouping_factor = groups,
  p_adjust_method = "BH",
  return_plots = TRUE,
  plot_type = "boxplot",
  filled = FALSE,
  border_width = 1.5,
  add_points = TRUE,
  point_type = "jitter",
  point_shape = 21,
  point_size = 2.5,
  point_color = "white",
  color_palette = c("#BC3C29", "#0072B5") # NEJM配色
)

# 查看结果
print(batch_results$results)
batch_results$plots$Bacteroides

# ============================================================================
# 示例10：条形图 + 淡色填充 + 自定义配色
# ============================================================================
result10 <- compare_species_abundance_v2(
  abundance_matrix = abundance,
  grouping_factor = groups,
  species_name = "Bacteroides",
  plot_type = "barplot",
  filled = FALSE, # 淡色填充（alpha=0.1）
  add_points = TRUE,
  point_type = "jitter",
  color_palette = c("#3B4992", "#EE0000") # Science配色
)
result10$plot

# ============================================================================
# 示例11：查看完整结果
# ============================================================================
result <- compare_species_abundance_v2(
  abundance_matrix = abundance,
  grouping_factor = groups,
  species_name = "Bacteroides",
  test_method = "wilcox",
  p_adjust_method = "BH",
  color_palette = c("#E64B35", "#4DBBD5")
)

# 查看所有结果
print(result$raw_pvalue)
print(result$adjusted_pvalue)
print(result$effect_size)
print(result$summary_stats)
result$plot


set.seed(123)
abundance <- data.frame(
  Bacteroides = c(rnorm(15, 0.3, 0.05), rnorm(15, 0.2, 0.05)),
  Firmicutes = c(rnorm(15, 0.4, 0.08), rnorm(15, 0.5, 0.08)),
  Proteobacteria = c(rnorm(15, 0.1, 0.03), rnorm(15, 0.15, 0.03))
)
groups <- factor(rep(c("健康", "IBD"), each = 15))

# ============================================================================
# 示例1：添加星号显著性标记（推荐）
# ============================================================================
result1 <- compare_species_abundance_v2(
  abundance_matrix = abundance,
  grouping_factor = groups,
  species_name = "Bacteroides",
  plot_type = "barplot",
  filled = TRUE,
  main_alpha = 0.9,
  add_points = TRUE,
  point_type = "jitter",
  add_sig_label = "asterisk", # 添加星号标记 ***, **, *, ns
  sig_label_size = 12, # 星号大小
  sig_hide_ns = TRUE # 隐藏非显著性标记
)
result1$plot
