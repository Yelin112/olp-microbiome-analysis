# 分析图表产出规范

> 记录 2026-07 期间在流式细胞术（`projects/flow_cytometry/2026-07/`）和 IHC（`projects/ihc/2026-07/`）两个独立分析中确定下来的图表产出习惯。
> 目的：即使中间隔了很久没碰这套流程，也能照着这份文档把新分析接上同一套风格，不用重新摸索。
> 这些是**目前用下来效果不错的默认做法，不是强制标准**——如果某次分析有更合适的方式，直接换掉即可，不用迁就这里的写法。

> **关于这份文档与 Claude 记忆的关系**：这里记录的规范已经同步保存进 Claude Code 的项目记忆（`feedback_figure_output_conventions.md` / `feedback_stat_comparison_design.md`），所以新开对话时 Claude 会自动带上这些习惯，不需要每次重新说明。但这套规范目前只覆盖了两次分析里实际遇到的情况，**并不完整**，以后大概率还会继续修改、补充或优化。如果你（或以后的 Claude）调整了这里的内容，最好也顺手让 Claude 同步更新一下对应的记忆文件，避免文档和记忆内容不一致。

## 1. 图表导出三格式

每张图默认同时导出三种格式：

| 格式 | 用途 |
|---|---|
| PDF | 矢量图，投稿/印刷 |
| PNG（300dpi） | 位图，快速预览、嵌入文档/聊天 |
| PPTX | 可编辑矢量，PowerPoint 里继续调整文字/颜色/图例位置 |

PPTX 通过 [`R/utils/export_pptx.R`](utils/export_pptx.R) 实现，基于 `officer` + `rvg::dml()`，把 ggplot 渲染成 Office DrawingML，图形和文字在 PPT 里都是可编辑对象（不是插入一张位图）：

```r
source("R/utils/export_pptx.R")

# 单图单页
save_plot_pptx(p, "output/figures/01_xxx.pptx",
                panel_width = 9, panel_height = 7, title = "图标题")

# 多图打包成一份多页 PPTX（每图一页）
save_plots_pptx(named_plot_list, "output/figures/all_metrics.pptx",
                 panel_width = 9, panel_height = 7)
```

`fix_panel = TRUE`（默认）时内部调用 `panel_fix.R` 的 `fix_panel_size()`，保证 PPTX 版本跟同一张图的 PDF/PNG 版本 panel 尺寸一致；如果传入的是 patchwork 拼图对象，改用 `fix_panel = FALSE`（按整图宽高直接导出，不强制统一子面板尺寸）。

## 2. 图片内文字一律用英文

标题、轴标签、图例、`plot_annotation()` 总标题——全部用英文，**即使脚本注释和 `report.md` 正文照常用中文**。这是刻意跟项目"默认中文"习惯反着来的例外，容易被顺手写错，所以单独强调。

```r
# metric_info 的 label 列、compare_plot() 的 title/xlab/ylab、labs(subtitle=...) 都写英文
compare_plot(..., title = "CD8+ T Cell (%)", xlab = "Group", ylab = "CD8+ T Cell (%)")
```

## 3. 组间比较的统计设计

- **比较对默认"各组 vs Control"**，不是相邻链式比较，也不是全两两比较：
  ```r
  key_pairs <- list(c("Control", "OXA"), c("Control", "OXA_Mn"), ...)
  ```
- **整体检验（Kruskal-Wallis）放在 subtitle，不要用 `stat_compare_means()` 画在图内**——后者用 `label.y.npc="top"` 时会跟两两比较的括号重叠：
  ```r
  p <- compare_plot(...) + labs(subtitle = sprintf("Kruskal-Wallis, p = %.4f", kw$p))
  ```
- **统计表里同时存原始 p 值和 BH 校正 p 值**（`p.adjust(method="BH")`），图上的显著性星号是**未校正**的探索性展示，`report.md` 里要写清楚这一点，避免被误读成经过多重检验校正。

## 4. 百分比类指标：Y 轴上限不能超过 100%

细胞占比、阳性区域占比这类指标物理上不可能超过 100%，但堆叠多层比较括号时 `y_expand` 的自动扩展可能把轴顶推到 100 以上（尤其是数据基线本身就偏高的指标）。非百分比指标（染色强度、积分光密度等）不适用这条，应该按数据自动缩放，不能瞎设上限。

做法：先按常规参数出图，检测实际计算出的轴顶是否超过该指标的上限；只有超过时才收紧参数重画并加硬上限。参考两个脚本里的 `build_stat_plot()`：

```r
build_stat_plot <- function(df_m, m_label, comparisons, cap = Inf) {
  make_plot <- function(step, y_exp) {
    compare_plot(df_m, value.var = "value", group.by = "group", strategy = "auto",
      add_stat = "wilcox.test", comparisons = comparisons, stat_label = "p.signif",
      hide_ns = FALSE, step_increase = step, y_expand = y_exp,
      palette = "NPG", theme_use = theme_pub_base, title = m_label, xlab = "Group", ylab = m_label)
  }
  p <- make_plot(0.14, 0.5)
  top <- ggplot_build(p)$layout$panel_params[[1]]$y.range[2]
  if (is.finite(cap) && top > cap) {
    p <- make_plot(0.09, 0.02) +
      scale_y_continuous(limits = c(0, cap), expand = expansion(mult = c(0.05, 0)))
  }
  p
}
```

`cap` 按指标区分（百分比指标传 100，其他传 `Inf` 即跳过检测），一般放在 `metric_info` tibble 的一列里统一管理。

## 5. 输出目录结构（当前默认，非强制）

```
projects/<项目>/<日期>/output/
├── figures/           # 所有图，文件名带数字前缀保证顺序（00_overview, 01_xxx, 02_xxx...）
├── stats/              # CSV 统计结果表（含原始 p 值 + BH 校正 p 值）
└── report.md           # Markdown 分析报告，引用 figures/ 里的图片路径
```

两次分析都用了这个结构且没有问题，可以作为起点；但这只是目前用着顺手的默认布局，不是规定，某次分析如果有更合适的组织方式，直接换。

## 6. 参考实现

完整可运行的参考脚本（包含以上全部约定的实际用法）：
- [`projects/flow_cytometry/2026-07/scripts/cell_type_comparison.R`](../projects/flow_cytometry/2026-07/scripts/cell_type_comparison.R)
- [`projects/ihc/2026-07/scripts/cd8_field_comparison.R`](../projects/ihc/2026-07/scripts/cd8_field_comparison.R)

配色/主题系统本身见 [`R/utils/README.md`](utils/README.md)；`compare_plot()` 函数参数详见 [`R/r_functions/lib/compare_plot/compare_plot_完整文档.md`](r_functions/lib/compare_plot/compare_plot_完整文档.md)。
