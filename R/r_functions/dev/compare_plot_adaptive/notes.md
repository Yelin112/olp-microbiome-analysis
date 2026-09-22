# 开发笔记：compare_plot_adaptive

## 项目信息
- **开始日期**: 2025-04-05
- **前身函数**: compare_plot_optimized.R（已在 lib/compare_plot/ 中独立维护）
- **目标**: 基于样本量的自适应可视化策略，管线式架构

---

## 版本历史
- **v1** (归档于 archive/): 初始版本，手动调整参数（compare_plot_optimized.R）
- **v2** (归档于 strategies/archive/): 自适应版本，每个策略独立生成完整 ggplot
- **v3** (当前): 管线式架构 — 策略只管 geom 层，共享功能由统一管线叠加

---

## v3 架构设计

```
compare_plot_adaptive(data, group_col, value_col, ...)
│
├─ 1. 参数校验
├─ 2. 数据预处理（NA 移除、因子转换）
├─ 3. 策略检测（detect_strategy）
├─ 4. 颜色配置（.get_colors）
├─ 5. 初始化 ggplot（统一 aes 映射）
├─ 6. 背景层（可选: 斑马纹）
├─ 7. 策略图层 ← build_layers_xxx() 只返回 geom 列表
├─ 8. 统计检验层 ← 新增: ggpubr t.test/wilcox/anova/kruskal
├─ 9. Y轴扩展（为显著性标注留空间）
├─ 10. 分面（可选: facet_wrap）
├─ 11. 主题与标签（scale_fill_manual + labs + theme）
└─ 12. 返回 list(plot, strategy, n_per_group, metadata)
```

### 关键设计决策

| 决策 | 理由 |
|---|---|
| 策略函数只返回 geom 层 | 新增策略只需写一个 build_layers_xxx.R，自动继承统计/分面/主题 |
| 管线统一处理共享功能 | 改一处（统计/主题/分面）全部策略生效，消除4×重复 |
| 使用 .data[[col]] 替代 !!rlang::sym() | 不需要额外 rlang 依赖，与 ggplot2 官方推荐一致 |
| detect_strategy 纯函数无副作用 | message 由主函数 verbose 统一控制，函数可单独测试 |
| 返回 list 而非纯 ggplot | 保留元数据供下游代码查询策略/样本量 |

---

## 文件结构

```
compare_plot_adaptive/
├── function.R                      # 管线式主函数（~270行）
├── strategies/
│   ├── detect_strategy.R           # 策略检测（纯逻辑）
│   ├── layers_pure_scatter.R       # n < 5: 散点 + crossbar + errorbar
│   ├── layers_bar_points.R         # 5 ≤ n < 10: 柱状图 + errorbar + 散点
│   ├── layers_boxplot_points.R     # 10 ≤ n < 20: 箱线图/小提琴 + 散点
│   ├── layers_pure_boxplot.R       # n ≥ 20: 纯箱线图/小提琴 + 均值点
│   └── archive/                    # v2 旧策略文件（已废弃）
├── test.R                          # 24 个测试用例
├── test_output/                    # 测试生成的图片
├── example.R                       # 使用示例
└── notes.md                        # 本文件
```

---

## 开发进度

### 已完成 ✅
- [x] Phase 1: 项目结构建立
- [x] Phase 2: 策略检测模块（v3 已移除 rlang 依赖）
- [x] Phase 3-6: 4 个策略图层构建器
- [x] Phase 7: 管线式主函数整合
- [x] Phase 8: 完整测试（24 个用例，含统计检验/分面/主题/背景）
- [x] Phase 9: 使用示例
- [x] Phase 10: v2 旧策略文件归档

### 待解决 / 未来扩展
- [ ] 组内簇状比较（fill.by 参数，需处理 position_dodge）
- [ ] 统计检验在分面模式下的行为验证（ggpubr facet 兼容性）
- [ ] 极端值对箱线图的影响评估
- [ ] 更多策略类型（raincloud plot, beeswarm, sina plot 等）
- [ ] 文档：完整参数手册（vignette）

---

## 新增策略的方法

在管线式架构下，新增策略只需 3 步：

1. 创建 `strategies/layers_<name>.R`，写 `build_layers_<name>()` 返回 geom 列表
2. 在 `function.R` 的 `switch()` 中加一行映射
3. 在参数 `strategy` 的 choices 中加名字

```r
# 示例：新增 raincloud plot 策略
# 1. strategies/layers_raincloud.R
build_layers_raincloud <- function(...) {
  list(
    geom_flat_violin(...),
    geom_boxplot(width = 0.15, ...),
    geom_jitter(...)
  )
}

# 2. function.R switch 中加:
#   raincloud = build_layers_raincloud(...)
```

统计检验、分面、主题、背景、颜色自动继承，无需在策略文件中重复。

---

## 参考

### 对话记录
- Claude chat (2025-04-05): v2 初始开发
- Claude chat (2026-06-26): v3 管线式架构重构

### 相关函数
- `lib/compare_plot/compare_plot_optimized.R`: 稳定版（功能更全，但单体架构）
- `R_FUNCTION_DEV_WORKFLOW.md`: 函数开发流程
