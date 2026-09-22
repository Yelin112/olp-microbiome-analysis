# ============================================================================
# check_conventions.R — 检查 R 代码是否符合项目规范（多机可用 / 配色主题统一）
# ============================================================================
#
# 用法:
#   Rscript R/check_conventions.R                       # 检查 dev/ 和 lib/
#   Rscript R/check_conventions.R R/r_functions/dev/xxx # 检查指定文件或目录
#   Rscript R/check_conventions.R --strict path         # WARN 也算失败（晋升 lib 前用）
#   source("R/check_conventions.R"); check_conventions("path")
#
# 级别:
#   ERROR — 换机器/服务器就会坏，必须修（绝对路径、sys.frame ofile、../ 相对 source、setwd 字面路径）
#   WARN  — 违反配色/主题约定，晋升到 lib 前应处理（写死十六进制颜色、自定义 %||%、
#           绘图函数缺 palette / theme_use 参数、palette 默认值不是 NULL）
#
# 说明: 基于文本模式的粗检查，会有少量误报；注释和 roxygen 行已跳过。
#       个别确需保留的写法，在该行末尾加 `# convention-ok` 即可豁免。
#       只读，不修改任何文件。退出码: 有 ERROR 为 1（--strict 时有 WARN 也为 1）。
#       utils/、archive/、_legacy/、tests/ 目录及 test*/example* 文件不检查颜色和签名规则。

.is_comment <- function(x) grepl("^\\s*#", x)

# 每条规则: id, level, 正则(perl), 说明, 生效范围（"all" / "func" = 仅函数文件）
.RULES <- list(
  list("abs-path", "ERROR", "[\"'][A-Za-z]:[/\\\\]{1,2}|[\"']/(home|Users|data|mnt)/",
       "写死的绝对路径（换机器即失效），改用 olp_path() 或相对 OLP_ROOT", "all"),
  list("ofile", "ERROR", "sys\\.frame\\(1\\)\\$ofile",
       "sys.frame(1)$ofile 在 Rscript/RStudio 中失效，改用 init.R", "all"),
  list("rel-source", "ERROR", "source\\(\\s*[\"'](\\.\\./|\\./)+[^\"']*[\"']",
       "依赖工作目录的相对 source，改用 load_utils()/load_lib()", "all"),
  list("setwd-literal", "ERROR", "setwd\\(\\s*[\"']",
       "setwd() 字面路径，改用 file.path(OLP_ROOT, ...)", "all"),
  list("hex-color", "WARN", "[\"']#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?[\"']",
       "写死的十六进制颜色，颜色应来自 palette_system（get_colors/scale_*_pub_*）", "func"),
  list("own-or-op", "WARN", "^`%\\|\\|%`\\s*<-",
       "自定义 %||%，改用 helpers.R 的", "func")
)

.norm <- function(x) gsub("\\\\", "/", x)

.check_file <- function(f) {
  lines <- tryCatch(readLines(f, warn = FALSE, encoding = "UTF-8"), error = function(e) character(0))
  if (length(lines) == 0) return(NULL)
  path <- .norm(f)
  is_func <- !grepl("^(test|example|Fortest)", basename(f), ignore.case = TRUE) &&
    !grepl("(^|/)(tests?|_legacy|archive|utils)/", path)
  skip <- .is_comment(lines) | grepl("#\\s*convention-ok", lines)
  out <- list()
  add <- function(i, lvl, id, msg) {
    out[[length(out) + 1]] <<- data.frame(file = f, line = i, level = lvl, rule = id,
                                          msg = msg, stringsAsFactors = FALSE)
  }
  for (r in .RULES) {
    if (r[[5]] == "func" && !is_func) next
    hit <- which(!skip & grepl(r[[3]], lines, perl = TRUE))
    for (i in hit) add(i, r[[2]], r[[1]], r[[4]])
  }
  # 绘图函数签名检查（仅函数文件；代码中出现 ggplot( 调用）
  code <- lines[!skip]
  if (is_func && any(grepl("(^|[^_.A-Za-z0-9])ggplot2?(::ggplot)?\\(|ggplot2::ggplot\\(", code))) {
    body <- paste(code, collapse = "\n")
    if (!grepl("palette\\s*=", body)) {
      add(1, "WARN", "no-palette-arg", "绘图函数应有 palette 参数（对接 palette_system）")
    } else if (grepl("palette\\s*=\\s*[\"']", body) && !grepl("palette\\s*=\\s*NULL", body)) {
      add(1, "WARN", "palette-default", "palette 默认值建议为 NULL，让全局默认生效")
    }
    if (!grepl("theme_use\\s*=", body)) {
      add(1, "WARN", "no-theme-arg", "绘图函数应有 theme_use 参数（对接 theme_system）")
    }
  }
  if (length(out) == 0) NULL else do.call(rbind, out)
}

check_conventions <- function(paths = NULL, strict = FALSE, root = NULL) {
  if (is.null(root)) {
    root <- Sys.getenv("OLP_ROOT", unset = "")
    if (!nzchar(root)) {
      root <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE),
                       error = function(e) "")
    }
    root <- .norm(root)
  }
  if (is.null(paths)) {
    if (!nzchar(root)) stop("找不到仓库根目录：请设置 OLP_ROOT，或显式传入路径")
    paths <- file.path(root, "R", "r_functions", c("dev", "lib"))
  }
  files <- unlist(lapply(paths, function(p) {
    if (dir.exists(p)) list.files(p, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
    else if (file.exists(p)) p else character(0)
  }))
  files <- files[!grepl("(^|/)(archive|_legacy)/", .norm(files))]
  parts <- lapply(files, .check_file)
  parts <- parts[!vapply(parts, is.null, TRUE)]
  res <- if (length(parts)) do.call(rbind, parts) else NULL

  cat(sprintf("检查 %d 个文件\n", length(files)))
  if (is.null(res) || nrow(res) == 0) {
    cat("未发现问题。\n")
    return(invisible(TRUE))
  }
  res <- res[order(res$level, res$file, res$line), ]
  for (i in seq_len(nrow(res))) {
    r <- res[i, ]
    shown <- if (nzchar(root)) sub(paste0("^", root, "/?"), "", .norm(r$file)) else .norm(r$file)
    cat(sprintf("[%s] %s:%d  (%s) %s\n", r$level, shown, r$line, r$rule, r$msg))
  }
  n_err <- sum(res$level == "ERROR"); n_warn <- sum(res$level == "WARN")
  cat(sprintf("\n%d ERROR, %d WARN\n", n_err, n_warn))
  invisible(n_err == 0 && (!strict || n_warn == 0))
}

if (sys.nframe() == 0 && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  strict <- "--strict" %in% args
  paths <- setdiff(args, "--strict")
  ok <- check_conventions(if (length(paths)) paths else NULL, strict = strict)
  quit(status = if (isTRUE(ok)) 0 else 1)
}
