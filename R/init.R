# ============================================================================
# init.R — 统一入口：定位仓库根目录，加载 utils 与 lib 函数
# ============================================================================
#
# 目的: 本地 Windows 与 Linux 服务器共用同一套脚本，脚本不再写 "../../R/..."
#       这类依赖目录深度的相对路径，也不依赖 sys.frame(1)$ofile（Rscript 下失效）。
#
# 用法（任何脚本的开头）:
#   source(file.path(Sys.getenv("OLP_ROOT"), "R/init.R"))
#   load_utils()                    # helpers + palette_system
#   load_utils(theme = TRUE)        # 另加 theme_system
#   load_lib("compare_plot")        # 加载 lib/compare_plot/ 下的函数
#
# 每台机器在 ~/.Renviron 中设置一次:
#   OLP_ROOT=/path/to/分析          # Windows 例: OLP_ROOT=E:/打工人/Zeng/2023-3 OLP Microbiome/分析
#
# 未设置 OLP_ROOT 时，会依次尝试: git rev-parse → 从当前目录向上查找。
# 若在未设置 OLP_ROOT 的情况下 source 本文件，请用 source("<相对路径>/R/init.R")，
# 加载后根目录会被自动解析并写入 OLP_ROOT。
#
# 依赖: 无（仅用 base R；git 可选）

# 仓库根目录的标志文件
.OLP_MARKER <- file.path("R", "utils", "palette_system.R")

.olp_is_root <- function(dir) {
  nzchar(dir) && file.exists(file.path(dir, .OLP_MARKER))
}

# 定位仓库根目录；找不到则报错
.olp_find_root <- function() {
  # 1. 环境变量
  env <- Sys.getenv("OLP_ROOT", unset = "")
  if (.olp_is_root(env)) {
    return(normalizePath(env, winslash = "/", mustWork = TRUE))
  }
  if (nzchar(env)) {
    warning("OLP_ROOT='", env, "' 下未找到 ", .OLP_MARKER, "，尝试自动查找")
  }

  # 2. git
  top <- tryCatch(
    suppressWarnings(system2(
      "git",
      c("rev-parse", "--show-toplevel"),
      stdout = TRUE,
      stderr = FALSE
    )),
    error = function(e) character(0)
  )
  if (length(top) == 1 && .olp_is_root(top)) {
    return(normalizePath(top, winslash = "/", mustWork = TRUE))
  }

  # 3. 从当前工作目录向上查找
  dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    if (.olp_is_root(dir)) {
      return(dir)
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) {
      break
    }
    dir <- parent
  }

  stop(
    "找不到项目根目录。请在 ~/.Renviron 中设置 OLP_ROOT=<仓库根目录>，",
    "或在仓库内运行脚本。"
  )
}

OLP_ROOT <- .olp_find_root()
Sys.setenv(OLP_ROOT = OLP_ROOT)

#' 拼接相对于仓库根目录的路径
#' @examples olp_path("projects", "16S", "data")
olp_path <- function(...) file.path(OLP_ROOT, ...)

#' 加载 utils 基础工具
#'
#' 始终加载 helpers.R 与 palette_system.R；其余按需。
#'
#' @param theme  是否加载 theme_system.R
#' @param panel  是否加载 panel_fix.R
#' @param pptx   是否加载 export_pptx.R（自动带上 panel_fix.R；缺 officer/rvg 时
#'   导出函数会给出提示并跳过，不会中断脚本）
load_utils <- function(theme = FALSE, panel = FALSE, pptx = FALSE) {
  u <- function(f) source(olp_path("R", "utils", f), chdir = FALSE)
  u("helpers.R")
  u("palette_system.R")
  if (theme) {
    u("theme_system.R")
  }
  if (panel || pptx) {
    u("panel_fix.R")
  }
  if (pptx) {
    u("export_pptx.R")
  }
  invisible(TRUE)
}

#' 加载 lib 中的函数
#'
#' @param name  lib 下的分类目录名（如 "abundance_processing"、"compare_plot"），
#'   或单个函数文件名（不含 .R，如 "beta_calc"）。目录会加载其中除 test_ 开头
#'   以外的全部 .R 文件。
#' @return 不可见地返回已加载的文件路径。
load_lib <- function(name) {
  lib <- olp_path("R", "r_functions", "lib")
  dir <- file.path(lib, name)

  if (dir.exists(dir)) {
    files <- list.files(dir, pattern = "\\.R$", full.names = TRUE)
    files <- files[!grepl("^test_", basename(files))]
  } else {
    files <- Sys.glob(file.path(lib, "*", paste0(name, ".R")))
  }

  if (length(files) == 0) {
    stop("lib 中找不到 '", name, "'（既不是分类目录，也不是函数文件名）")
  }

  # lib 函数依赖 utils；未加载时补加载
  if (!exists("get_colors", mode = "function") || !exists("%||%", mode = "function")) {
    load_utils()
  }
  for (f in files) {
    source(f)
  }
  invisible(files)
}
