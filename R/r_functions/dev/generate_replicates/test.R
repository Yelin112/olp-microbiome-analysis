# ============================================================================
# 测试文件：generate_replicates
# 创建日期：2026-08-30
#
# 用法：source() 本文件，自动运行全部测试并打印 pass/fail 汇总。
# ============================================================================

FUNC_DIR <- "e:/打工人/Zeng/2023-3 OLP Microbiome/分析/R/r_functions/dev/generate_replicates"
source(file.path(FUNC_DIR, "function.R"))

# ---- 测试框架 ----

pass_count <- 0
fail_count <- 0

run_test <- function(name, expr) {
  cat("▶ ", name, "... ")
  tryCatch(
    {
      expr
      cat("✅ 通过\n")
      pass_count <<- pass_count + 1
    },
    error = function(e) {
      cat("❌ 失败:", conditionMessage(e), "\n")
      fail_count <<- fail_count + 1
    }
  )
}

# 静默调用（吞掉安全阀 warning），但保留返回值
quiet <- function(expr) suppressWarnings(expr)

# ============================================================================
# 基础功能
# ============================================================================

run_test("continuous: 正常值", {
  r <- quiet(generate_replicates(100, type = "continuous", n = 5, seed = 1))
  stopifnot(length(r) == 5, all(r > 0))
})

run_test("continuous: 负值（应 warning，不 stop）", {
  r <- quiet(generate_replicates(-10, type = "continuous", n = 3, seed = 1))
  stopifnot(length(r) == 3)
})

run_test("continuous: 0（应返回 n 个 0）", {
  r <- quiet(generate_replicates(0, type = "continuous", n = 3))
  stopifnot(all(r == 0))
})

run_test("proportion: 0~1 输入", {
  r <- quiet(generate_replicates(0.75, type = "proportion", n = 5, seed = 1, scale = "0-1"))
  stopifnot(length(r) == 5, all(r > 0 & r < 1))
})

run_test("proportion: 0~100 输入", {
  r <- quiet(generate_replicates(75, type = "proportion", n = 5, seed = 1, scale = "0-100"))
  stopifnot(length(r) == 5, all(r > 0 & r < 100))
})

run_test("proportion: auto 检测 >1 时按百分比处理", {
  r <- quiet(generate_replicates(75, type = "proportion", n = 5, seed = 1, scale = "auto"))
  stopifnot(all(r > 0 & r < 100))
})

run_test("proportion: 边界值 0 应 stop", {
  ok <- tryCatch({
    quiet(generate_replicates(0, type = "proportion", n = 3, scale = "0-1"))
    FALSE
  }, error = function(e) TRUE)
  stopifnot(ok)
})

run_test("proportion: 边界值 1（100%）应 stop", {
  ok <- tryCatch({
    quiet(generate_replicates(1, type = "proportion", n = 3, scale = "0-1"))
    FALSE
  }, error = function(e) TRUE)
  stopifnot(ok)
})

run_test("abundance: 正常值", {
  r <- quiet(generate_replicates(0.1, type = "abundance", n = 5, seed = 1))
  stopifnot(length(r) == 5, all(r >= 0 & r <= 1))
})

run_test("abundance: 极小值", {
  r <- quiet(generate_replicates(0.001, type = "abundance", n = 5, seed = 1))
  stopifnot(length(r) == 5, all(r >= 0))
})

run_test("abundance: 0（应返回 n 个 0）", {
  r <- quiet(generate_replicates(0, type = "abundance", n = 3))
  stopifnot(all(r == 0))
})

run_test("count: 整数", {
  r <- quiet(generate_replicates(50, type = "count", n = 5, seed = 1))
  stopifnot(length(r) == 5, all(r >= 0))
})

run_test("count: 非整数（不再强制 round，允许非整数 mu）", {
  r <- quiet(generate_replicates(50.5, type = "count", n = 5, seed = 1))
  stopifnot(length(r) == 5, all(r >= 0))
})

# ============================================================================
# 输入格式
# ============================================================================

run_test("单值输入 -> 返回向量", {
  r <- quiet(generate_replicates(0.5, type = "continuous", n = 4))
  stopifnot(is.numeric(r), is.null(dim(r)), length(r) == 4)
})

run_test("向量输入 -> 返回矩阵（行=数据点，列=重复）", {
  r <- quiet(generate_replicates(c(0.12, 0.08, 0.15), type = "abundance", n = 3, seed = 1))
  stopifnot(is.matrix(r), nrow(r) == 3, ncol(r) == 3)
})

run_test("data.frame + col -> 返回扩增后数据框", {
  df <- data.frame(sample = c("A", "B"), expression = c(10, 20))
  out <- quiet(generate_replicates(df, type = "continuous", n = 3, col = "expression", seed = 1))
  stopifnot(all(c("expression_rep1", "expression_rep2", "expression_rep3") %in% names(out)))
  stopifnot(nrow(out) == 2, "expression" %in% names(out))
})

# ============================================================================
# 参数行为
# ============================================================================

run_test("seed 固定时结果可重复", {
  r1 <- quiet(generate_replicates(0.5, type = "continuous", n = 5, seed = 42))
  r2 <- quiet(generate_replicates(0.5, type = "continuous", n = 5, seed = 42))
  stopifnot(identical(as.numeric(r1), as.numeric(r2)))
})

run_test("seed 不污染全局随机状态", {
  set.seed(999)
  before <- get(".Random.seed", envir = .GlobalEnv)
  quiet(generate_replicates(0.5, type = "continuous", n = 5, seed = 1))
  after <- get(".Random.seed", envir = .GlobalEnv)
  stopifnot(identical(before, after))
})

run_test("clip=TRUE 裁剪越界值到 [0,1]", {
  r <- quiet(generate_replicates(0.99, type = "abundance", n = 50, cv = 0.5, clip = TRUE, seed = 1))
  stopifnot(all(r >= 0 & r <= 1))
})

run_test("clip=FALSE 允许越界（仅 warning）", {
  r <- suppressWarnings(generate_replicates(0.99, type = "abundance", n = 50, cv = 0.5, clip = FALSE, seed = 1))
  stopifnot(is.numeric(r))
})

run_test("n=1 时向量输入仍返回矩阵（类型稳定）", {
  r <- quiet(generate_replicates(c(0.1, 0.2), type = "continuous", n = 1))
  stopifnot(is.matrix(r), nrow(r) == 2, ncol(r) == 1)
})

# ============================================================================
# 错误处理
# ============================================================================

run_test("proportion 传入 0 或 1 应 stop（已覆盖于上方）", {
  stopifnot(TRUE)
})

run_test("data.frame 缺少 col 参数应 stop", {
  df <- data.frame(x = 1:3)
  ok <- tryCatch({ generate_replicates(df, type = "continuous", n = 3); FALSE },
                  error = function(e) TRUE)
  stopifnot(ok)
})

# ============================================================================
# 安全阀：每次调用都必须 warning + 打标记
# ============================================================================

run_test("安全阀：调用产生 warning", {
  triggered <- FALSE
  withCallingHandlers(
    generate_replicates(0.5, type = "continuous", n = 3),
    warning = function(w) { triggered <<- TRUE; invokeRestart("muffleWarning") }
  )
  stopifnot(triggered)
})

run_test("安全阀：返回值带 synthetic attribute", {
  r <- quiet(generate_replicates(0.5, type = "continuous", n = 3))
  stopifnot(isTRUE(attr(r, "synthetic")))
  stopifnot(attr(r, "synthetic_type") == "continuous")
})

# ============================================================================
# 汇总
# ============================================================================

cat("\n========================================\n")
cat("通过:", pass_count, " / 失败:", fail_count, "\n")
cat("========================================\n")
