# 开发笔记：pcoa_plot

## 项目信息
- **开始日期**：2026-04-06
- **归档日期**：2026-04-06
- **归档路径**：`r_functions/lib/pcoa_plot.R`
- **定位**：通用二维散点图，支持密度曲线和箱线图两种边际图样式

---

## 参考图样式

### 样式A（density 边际）
- 主图：散点 + 置信椭圆（填充 + 虚线轮廓）
- 边际：顶部 + 右侧核密度曲线，按组着色
- 关键包：`ggExtra::ggMarginal()`

### 样式B（boxplot 边际）
- 主图：散点，多行斜体统计注释
- 右侧面板：y 轴变量分组箱线图 + 显著性字母
- 顶部面板：x 轴变量分组箱线图（横向）
- 关键包：`aplot::insert_right()` + `aplot::insert_top()`

---

## 当前进度
- [x] function.R — 主函数完成
- [x] test.R — 20 个测试用例全部通过
- [x] 归档至 lib/pcoa_plot.R

---

## 设计决策记录

| 决策点 | 结论 |
|--------|------|
| 命名向量配色 | 有 names 时保留用户映射，无 names 时按分组水平顺序 |
| RColorBrewer 支持 | 识别全部 35 种色板，未知名称给 warning 回退 NPG |
| sig_letters 位置 | 各组 Q75 + 6% range offset 自动计算 |
| boxplot 轴对齐 | 右侧面板共享 y 轴；顶部面板 coord_flip() 共享 x 轴 |

---

## 依赖包
- `ggplot2`（必需）
- `ggExtra`（density 边际）
- `aplot`（boxplot 边际）
- `RColorBrewer`（可选，非 NPG 色板）
