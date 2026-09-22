# ============================================================================
# helpers.R — 跨项目共享的操作符和工具函数
# 版本: 1.0
# ============================================================================
#
# 用法: source("R/utils/helpers.R")
# 依赖: 无
#
# 设计原则:
#   - 本文件零依赖，可安全作为任何脚本的第一个 source
#   - 只放极小、极通用的工具函数
#   - 永远不在这里定义调色板、主题、绘图逻辑

# ---- 空值合并操作符 ----

#' 空值合并（a %||% b）
#'
#' 如果 a 为 NULL 则返回 b，否则返回 a。
#' 类似 SQL COALESCE 或 JS ?? 操作符。
#'
#' @param a 主值
#' @param b 回退值
#'
#' @examples
#' NULL %||% "default"   # → "default"
#' "hello" %||% "world"  # → "hello"
`%||%` <- function(a, b) if (is.null(a)) b else a


# ---- 反向 %in% ----

#' 反向匹配（x %ni% table）
#'
#' 等价于 !x %in% table，但更可读。
#'
#' @param x 要查找的值
#' @param table 查找表
#'
#' @examples
#' "z" %ni% letters[1:10]  # → TRUE
`%ni%` <- Negate(`%in%`)


# ---- 安全 source ----

#' 文件存在时才 source
#'
#' 避免因路径错误导致的 invisible 失败。
#'
#' @param path 文件路径
#' @param verbose 是否打印信息
#'
#' @examples
#' source_if_exists("optional_config.R")
source_if_exists <- function(path, verbose = FALSE) {
  if (file.exists(path)) {
    if (verbose) {
      message("→ 加载: ", path)
    }
    source(path)
    return(invisible(TRUE))
  } else {
    if (verbose) {
      message("⊘ 跳过 (不存在): ", path)
    }
    return(invisible(FALSE))
  }
}
