# ============================================================================
# 可视化函数标准模板 — 新绘图函数的开发起点
# 版本: 1.0
# 最后更新: 2026-06-27
# ============================================================================
#
# 使用方法:
#   1. 复制本文件到 projects/<项目>/scripts/ 或 R/r_functions/dev/
#   2. 替换所有 "my_function" 为你自己的函数名
#   3. 按需增删参数
#   4. 填入业务逻辑
#
# 设计原则:
#   - 参数命名: 与 compare_plot 系列保持一致 (palette / theme_use / add_stat)
#   - 配色: 永远通过 scale_fill_pub_d() 或 get_colors() 获取，不在函数内定义调色板
#   - 主题: 永远通过 theme_use 参数传入，不在函数内硬编码 theme_classic()
#   - 统计: 如涉及假设检验，统一用 add_stat + stat_label 参数对
#   - 返回: 推荐返回 list(plot = ..., data = ..., stats = ...)，保留元数据

# ============================================================================
# 0. 依赖声明
# ============================================================================

library(ggplot2)
library(dplyr)

# 工具层：不在函数文件里写死路径，也不拼 "../../utils"。
# 由调用方先 source R/init.R 并 load_utils()；这里只做兜底（依赖 OLP_ROOT 环境变量）。
# 不要用 sys.frame(1)$ofile（Rscript 下失效），不要写 E:/ 或 /home/ 绝对路径。
if (!exists("%||%") || !exists("get_colors")) {
  .init <- file.path(Sys.getenv("OLP_ROOT"), "R", "init.R")
  if (!file.exists(.init)) {
    stop("未找到 utils 函数。请先 source R/init.R 并调用 load_utils()，或设置 OLP_ROOT 环境变量。")
  }
  source(.init)
  load_utils(theme = TRUE)   # helpers + palette_system + theme_system
  rm(.init)
}

# ============================================================================
# 1. 函数签名
# ============================================================================

#' 组间比较可视化（示例）
#'
#' @description 一段话描述这个函数做什么。
#'
#' @param data        data.frame，包含分组列和数值列。
#' @param group_col   分组列名，字符串。
#' @param value_col   数值列名，字符串。
#' @param palette     配色方案。支持:
#'   \itemize{
#'     \item 字符串 — 系统注册的色板名 (如 \code{"NPG"}, \code{"JAMA"})
#'     \item 颜色向量 — 自定义颜色 (如 \code{c("#FF0000", "#00FF00")})
#'     \item NULL — 使用全局默认色板
#'   }
#'   详见 \code{?palette_system}。
#' @param theme_use   ggplot2 主题函数。
#'   内置: \code{theme_nature()}, \code{theme_cell()}, \code{theme_pub_base()},
#'   \code{theme_pub_stat()}。也可传入任意 ggplot2 主题 (如 \code{theme_bw})。
#'   NULL 时使用 ggplot2 默认主题。
#' @param add_stat    统计检验方法。\code{"none"} / \code{"t.test"} /
#'   \code{"wilcox.test"} / \code{"anova"} / \code{"kruskal.test"}。
#' @param stat_label  统计标签格式: \code{"p.signif"} (星号) /
#'   \code{"p.format"} ("p = 0.001")。
#' @param comparisons 手动指定比较组 list，如 \code{list(c("A","B"), c("A","C"))}。
#'   NULL 时自动生成全部两两组合。
#' @param xlab, ylab, title  轴标签与标题。NULL 时自动从列名生成。
#' @param verbose     是否打印运行信息，默认 TRUE。
#' @param ...         传递给 theme_use 的额外参数。
#'
#' @return 列表: list(plot = ggplot, data = data.frame, stats = data.frame)
#' @export
#'
#' @examples
#' # 基础用法
#' my_function(iris, "Species", "Sepal.Length")
#'
#' # 自定义配色 + 主题 + 统计检验
#' my_function(iris, "Species", "Sepal.Length",
#'             palette = "JAMA", theme_use = theme_nature(),
#'             add_stat = "anova")
my_function <- function(
    data,
    group_col,
    value_col,
    # ── 视觉参数 ──
    palette    = NULL,
    theme_use  = NULL,
    # ── 统计参数 ──
    add_stat   = c("none", "t.test", "wilcox.test", "anova", "kruskal.test"),
    stat_label = c("p.signif", "p.format"),
    comparisons = NULL,
    hide_ns    = FALSE,
    step_increase = 0.12,
    # ── 标签 ──
    xlab  = NULL,
    ylab  = NULL,
    title = NULL,
    # ── 控制 ──
    verbose = TRUE,
    ...
) {
  # ==========================================================================
  # 阶段 1: 参数校验
  # ==========================================================================

  add_stat   <- match.arg(add_stat)
  stat_label <- match.arg(stat_label)

  if (!is.data.frame(data))
    stop("data 必须是 data.frame")
  if (!group_col %in% names(data))
    stop("列 '", group_col, "' 在数据框中不存在")
  if (!value_col %in% names(data))
    stop("列 '", value_col, "' 在数据框中不存在")

  # ==========================================================================
  # 阶段 2: 数据预处理
  # ==========================================================================

  n_before <- nrow(data)
  data <- data[!is.na(data[[group_col]]) & !is.na(data[[value_col]]), ]
  if (verbose && n_before - nrow(data) > 0) {
    message("已移除 ", n_before - nrow(data), " 行缺失值")
  }
  if (nrow(data) == 0)
    stop("移除缺失值后数据为空")

  if (!is.factor(data[[group_col]])) {
    data[[group_col]] <- factor(data[[group_col]],
                                 levels = unique(data[[group_col]]))
  }

  # ==========================================================================
  # 阶段 3: 初始化 ggplot（统一 aes 映射）
  # ==========================================================================

  p <- ggplot(data, aes(
    x    = .data[[group_col]],
    y    = .data[[value_col]],
    fill = .data[[group_col]]
  ))

  # ==========================================================================
  # 阶段 4: 主体图层
  # ==========================================================================

  p <- p +
    geom_boxplot(                    # ← 替换为你的主图类型
      alpha     = 0.7,
      linewidth = 0.5,
      color     = "black"
    )

  # ==========================================================================
  # 阶段 5: 配色 — 统一对接点
  # ==========================================================================

  p <- p + scale_fill_pub_d(palette)

  # ==========================================================================
  # 阶段 6: 统计检验 — 统一对接点（可选）
  # ==========================================================================

  if (add_stat != "none") {
    # 自动生成两两比较
    comps <- comparisons
    n_groups <- nlevels(data[[group_col]])
    if (is.null(comps) && n_groups > 2 &&
        add_stat %in% c("t.test", "wilcox.test")) {
      comps <- utils::combn(levels(data[[group_col]]), 2, simplify = FALSE)
    }

    p <- p + ggpubr::stat_compare_means(
      comparisons   = comps,
      method        = add_stat,
      label         = if (is.null(comps) && n_groups > 2) "p.format" else stat_label,
      hide.ns       = hide_ns,
      step.increase = step_increase,
      size          = 3.5,
      bracket.size  = 0.5
    )

    # Y 轴扩展（为显著性标注留空间）
    p <- p + scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.15))
    )
  }

  # ==========================================================================
  # 阶段 7: 标签
  # ==========================================================================

  p <- p + labs(
    title = title,
    x     = xlab %||% group_col,
    y     = ylab %||% value_col
  )

  # ==========================================================================
  # 阶段 8: 主题 — 统一对接点
  # ==========================================================================

  if (is.function(theme_use)) {
    p <- p + theme_use(...)
  }

  # ==========================================================================
  # 阶段 9: 返回
  # ==========================================================================

  return(list(
    plot = p,
    data = data
  ))
}


# ============================================================================
# 附录 A: 常用参数组合速查
# ============================================================================

# A1 — 探索性分析（快速出图，不关心美观）
# my_function(data, "group", "value")

# A2 — 发表级（Nature 风格 + NPG 配色 + 统计检验）
# my_function(data, "group", "value",
#             palette = "NPG", theme_use = theme_nature(),
#             add_stat = "wilcox.test")

# A3 — 组内比较（同一个 palette，多组时自动插值）
# my_function(data, "time", "value",
#             palette = c("#E64B35", "#4DBBD5", "#00A087"),
#             theme_use = theme_pub_base(border = TRUE))

# A4 — 多面板
# my_function(data, "group", "value",
#             palette = "JCO", theme_use = theme_pub_stat()) +
#   facet_wrap(~ batch, scales = "free")


# ============================================================================
# 附录 B: 参数命名规范（所有上层函数必须遵守）
# ============================================================================

# ┌──────────────┬─────────────────────────────────────────┐
# │ 参数名        │ 语义                                    │
# ├──────────────┼─────────────────────────────────────────┤
# │ palette      │ 配色方案（名/向量/NULL）                  │
# │ theme_use    │ 主题函数（或 NULL）                       │
# │ add_stat     │ 统计检验方法                             │
# │ stat_label   │ 统计标签格式                             │
# │ comparisons  │ 手动指定比较组                           │
# │ hide_ns      │ 隐藏非显著                               │
# │ data         │ 输入数据框（永远第一个参数）               │
# │ group_col    │ 分组列名（字符串，不是裸名）               │
# │ value_col    │ 数值列名（字符串，不是裸名）               │
# │ xlab/ylab    │ 轴标签（NULL = 自动生成）                 │
# │ title        │ 图标题（NULL = 无标题）                   │
# │ verbose      │ 是否打印运行信息                          │
# │ ...          │ 传递给 theme_use() 的额��参数              │
# └──────────────┴─────────────────────────────────────────┘


# ============================================================================
# 附录 C: 目录布局建议
# ============================================================================

# R/
# ├── utils/
# │   ├── helpers.R              # %||%, %ni%, source_if_exists()
# │   ├── palette_system.R       # 颜色注册中心 + 解析引擎 + scale 函数
# │   ├── theme_system.R         # 主题工厂 + 期刊预设
# │   ├── panel_fix.R            # 固定面板尺寸工具
# │   └── FUNCTION_TEMPLATE.R    # 本文件
# ├── r_functions/
# │   ├── dev/                   # 开发中的函数
# │   ├── lib/                   # 稳定函数（source 目标）
# │   └── archive/               # 旧版本
# └── scratch/                   # 一次性/实验脚本


# ============================================================================
# 附录 D: 从旧函数迁移检查清单
# ============================================================================

# 迁移一个旧函数到新标准时，逐项核对:
#
# [ ] 删除了函数内的调色板列表（改用 get_colors / scale_fill_pub_d）
# [ ] 删除了函数内的 theme_classic() 等硬编码主题（改用 theme_use 参数）
# [ ] palette 参数默认值改为 NULL（而非 "NPG"）——让全局默认生效
# [ ] theme_use 参数默认值改为 NULL（而非 theme_classic）
# [ ] 工具层改用 init.R 兜底加载（见第 0 节），没有 "../../utils" 或绝对路径
# [ ] 运行 Rscript R/check_conventions.R <文件或目录> 无 ERROR
# [ ] 删除了自定义的 %||% 定义（改用 helpers.R 的）
# [ ] 返回格式统一为 list(plot, data, ...)
# [ ] 添加了完整的 @param 文档（特别是 palette 和 theme_use）
