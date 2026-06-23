# ==============================================================================
# registry.R: plot_enrich() 主函数 + 可视化类型注册表
# 创建: 2026-04-24
#
# 职责:
#   - 维护已注册的可视化类型（type → 绘图函数）
#   - 提供 plot_enrich() 统一调用入口
#   - 提供 register_enrich_plot() 扩展接口
# ==============================================================================


# ------------------------------------------------------------------------------
# 内部注册表：存储 type → list(fn, description, required_cols)
# ------------------------------------------------------------------------------

.ENRICH_PLOT_REGISTRY <- new.env(parent = emptyenv())


# ------------------------------------------------------------------------------
# register_enrich_plot(): 注册一个可视化类型
# ------------------------------------------------------------------------------

#' 注册富集分析可视化类型
#'
#' 将自定义绘图函数注册到 plot_enrich() 的调度系统中，
#' 注册后可直接通过 plot_enrich(data, type = "your_type") 调用。
#'
#' @param type     类型名称字符串，如 "bubble"、"bar"、"my_style"
#' @param fn       绘图函数，签名必须为 fn(data, ...)，返回 ggplot 对象
#' @param description 简短描述（用于 list_enrich_plots() 展示）
#' @param overwrite 是否覆盖已注册的同名类型，默认 FALSE
#'
#' @examples
#' \dontrun{
#' my_plot_fn <- function(data, ...) {
#'   std <- as_enrich_std(data)
#'   ggplot2::ggplot(std, ggplot2::aes(x = gene_count, y = term_name)) +
#'     ggplot2::geom_col()
#' }
#' register_enrich_plot("my_style", my_plot_fn, "自定义条形图")
#' plot_enrich(kegg_result, type = "my_style")
#' }
#'
#' @export
register_enrich_plot <- function(type, fn, description = "", overwrite = FALSE) {
  if (!is.character(type) || length(type) != 1 || nchar(type) == 0) {
    stop("type 必须是长度为 1 的非空字符串")
  }
  if (!is.function(fn)) {
    stop("fn 必须是函数")
  }
  if (exists(type, envir = .ENRICH_PLOT_REGISTRY, inherits = FALSE) && !overwrite) {
    stop("类型 '", type, "' 已注册。使用 overwrite = TRUE 强制覆盖，",
         "或先调用 unregister_enrich_plot('", type, "')")
  }

  assign(type, list(fn = fn, description = description), envir = .ENRICH_PLOT_REGISTRY)
  invisible(type)
}


#' 注销一个已注册的可视化类型
#' @export
unregister_enrich_plot <- function(type) {
  if (!exists(type, envir = .ENRICH_PLOT_REGISTRY, inherits = FALSE)) {
    warning("类型 '", type, "' 未注册，无需注销")
    return(invisible(NULL))
  }
  rm(list = type, envir = .ENRICH_PLOT_REGISTRY)
  invisible(type)
}


#' 列出所有已注册的可视化类型
#' @export
list_enrich_plots <- function() {
  types <- ls(.ENRICH_PLOT_REGISTRY)
  if (length(types) == 0) {
    cat("当前没有已注册的可视化类型。\n")
    return(invisible(character(0)))
  }
  cat("已注册的可视化类型:\n")
  for (t in sort(types)) {
    desc <- .ENRICH_PLOT_REGISTRY[[t]]$description
    cat(sprintf("  %-15s %s\n", t, desc))
  }
  invisible(types)
}


# ------------------------------------------------------------------------------
# plot_enrich(): 统一绘图入口
# ------------------------------------------------------------------------------

#' 富集分析可视化（统一入口）
#'
#' 根据 type 参数选择可视化样式，将富集分析结果绘制为图形。
#' 支持的 type 取决于已注册的绘图函数，使用 list_enrich_plots() 查看。
#'
#' @param data 富集分析结果，接受：
#'   - EnrichStd 对象（由 as_enrich_std() 转换）
#'   - enrichResult / gseaResult（clusterProfiler 直接输出）
#'   - data.frame（标准列名或配合 col_map）
#' @param type 可视化类型，默认 "bubble"。
#'   使用 list_enrich_plots() 查看所有可用类型。
#' @param col_map 列名映射（data 为非标准 data.frame 时使用），
#'   见 as_enrich_std() 的 col_map 参数说明
#' @param ... 传递给对应绘图函数的其他参数
#'   （如 top_n、filter_padj、x_var、color_high 等）
#'
#' @return ggplot 对象
#'
#' @examples
#' \dontrun{
#' source("core/std_format.R")
#' source("core/registry.R")
#' # 加载样式（会自动注册）
#' source("styles/bubble.R")
#' source("styles/bar.R")
#'
#' # 查看可用类型
#' list_enrich_plots()
#'
#' # 使用默认气泡图
#' plot_enrich(kegg_result)
#'
#' # 切换样式
#' plot_enrich(kegg_result, type = "bar", top_n = 20)
#'
#' # 传递样式特定参数
#' plot_enrich(kegg_result, type = "bubble",
#'             x_var = "rich_factor", color_high = "#1a5276")
#' }
#'
#' @export
plot_enrich <- function(data, type = "bubble", col_map = NULL, ...) {
  # 检查类型是否已注册
  if (!exists(type, envir = .ENRICH_PLOT_REGISTRY, inherits = FALSE)) {
    available <- sort(ls(.ENRICH_PLOT_REGISTRY))
    if (length(available) == 0) {
      stop("没有已注册的可视化类型。请先 source 对应的样式文件，",
           "如 source('styles/bubble.R')")
    }
    stop(
      "未知的可视化类型: '", type, "'\n",
      "可用类型: ", paste(available, collapse = ", "), "\n",
      "使用 list_enrich_plots() 查看详情"
    )
  }

  # 调用对应绘图函数
  plot_fn <- .ENRICH_PLOT_REGISTRY[[type]]$fn
  plot_fn(data, col_map = col_map, ...)
}
