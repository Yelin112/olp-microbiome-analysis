# 开发笔记：compare_plot_adaptive

## 项目信息
- **开始日期**: 2025-04-05
- **前身函数**: compare_plot_optimized.R（已归档至 archive/）
- **目标**: 在 compare_plot 基础上，添加基于样本量的自适应可视化策略自动切换

---

## 版本历史
- **v1** (归档于 archive/): 初始版本，手动调整参数（compare_plot_optimized.R）
- **v2** (当前开发): 自适应版本，自动策略切换

---

## 开发目标

### 核心需求
根据每组样本量自动选择可视化策略：
- n < 5: 纯散点图
- 5 ≤ n < 10: 柱状图 + 散点
- 10 ≤ n < 20: 箱线图/小提琴图 + 散点
- n ≥ 20: 纯箱线图/小提琴图

### 设计原则
1. 以 max_n（最大组样本量）决定策略
2. 用户可自定义阈值
3. 支持手动覆盖
4. 尽量向后兼容（如果可能）

---

## 参考资料

### 对话记录
- Claude chat conversation (2025-04-05)
- 需求文档: R_FUNCTION_DEV_WORKFLOW.md

### 参考图风格
1. **纯散点**: 大散点 (size=3) + 均值误差线 (crossbar, size=1.2)，黑色边框 (stroke=0.5)
2. **柱状图+散点**: 柱宽 0.7，黑色散点，颜色: #4DBBD5, #E64B35, #00A087, #F39B7F
3. **箱线图+散点**: 箱线图 alpha=0.7，小散点 (size=1.5, alpha=0.6)，不重复显示离群点

---

## 当前进度

### 已完成 ✅
- [x] Phase 1: 项目结构建立
- [x] Phase 2: 策略检测模块 (detect_strategy.R)
- [x] Phase 3: 策略1实现 (pure_scatter.R)
- [x] Phase 4: 策略2实现 (bar_points.R)
- [x] Phase 5: 策略3实现 (boxplot_points.R)
- [x] Phase 6: 策略4实现 (pure_boxplot.R)
- [x] Phase 7: 主函数整合 (function.R)
- [x] Phase 8: 完整测试 (test.R) — 12个测试用例
- [x] Phase 9: 使用示例 (example.R) — 6个示例
- [ ] Phase 10: 归档决策（待用户确认后执行）

---

## 技术细节

### 依赖包
```r
library(dplyr)      # 数据处理
library(ggplot2)    # 绘图
library(rlang)      # NSE处理
```

### 关键函数
1. `detect_strategy()`: 样本量检测和策略选择
2. `plot_pure_scatter()`: 纯散点绘图
3. `plot_bar_points()`: 柱状图+散点绘图
4. `plot_boxplot_points()`: 箱线图+散点绘图
5. `plot_pure_boxplot()`: 纯箱线图绘图
6. `compare_plot_adaptive()`: 主函数

---

## 已知问题

### 待解决
- [ ] 缺失值处理策略（删除 vs 提示）
- [ ] 极端值对箱线图的影响
- [ ] 统计检验集成（ggpubr vs ggsignif）
- [ ] 是否需要支持分面（facet）

### 已解决 ✅
- [x] NSE变量处理 (使用 !!sym())
- [x] 策略判断逻辑

---

## 笔记
- 用户希望在VSCode中使用Claude Code插件进行开发
- 遵循 r-function-dev workflow
- 强调不要一次性重写，而是逐步测试
