# enrich_clusters_pro — 开发笔记

## 当前状态
流式数据参数优化已完成，进入测试验证阶段。

## 核心逻辑
对 2D 坐标做向心收缩变换：`new_pos = old_pos + α × (center - old_pos)`
- linear 模式：α = strength（常数）
- gaussian 模式：α = strength × exp(−dist² / 2radius²)（随距离衰减）

## 参数速查

| 参数 | logicle 默认 | linear 默认 | 说明 |
|------|------------|------------|------|
| `strength` | 0.35 | 0.50 | 收缩强度，0~1 |
| `radius` | 自动估算 | 自动估算 | 高斯引力半径 |
| `min_events` | 50 | 50 | 稀有群体保护阈值 |

`radius` 自动估算逻辑：有 group_col 时 = 群体重心间平均距离 / 4；
无 group_col 时 = 数据范围 / 8，或 logicle 模式固定 0.5。

## 已解决的流式适配问题

- [x] `radius` 默认值对 logicle 坐标系不适配 → 自动估算 / coord_type 感知
- [x] 无坐标系感知，strength 默认值粗糙 → coord_type 参数 + 自动检测
- [x] 无稀有群体保护 → min_events 参数
- [x] umap/force 适用场景不清晰 → 添加警告 message

## 已知问题 / 待修

### 1. centers 列名假设
按组精准匹配（B 模式）用 `setdiff(colnames, group_col)[1:2]` 取坐标列，
假定非分组的前两列为 x/y。列顺序不同时可能出错。
→ 改进方向：增加 `x_center_col` / `y_center_col` 显式参数

### 2. umap/force 无 coord_type 感知
这两种模式输出坐标与原始坐标系完全无关，没有对应的 rescale 逻辑。
→ 改进方向：增加 `rescale = TRUE` 参数将输出缩放到原始数据范围

### 3. strength 上限无警告
strength > 0.85 时群体会接近坍缩为单点。
→ 改进方向：加 `if (strength > 0.85) warning(...)`

## 待开发功能

- [ ] 与 `palette_system.R` 对接：`palette` 参数 + 内置绘图（`plot = TRUE`）
- [ ] 多维扩展：当前写死 2D（new_x/new_y）
- [ ] strength 上限警告
- [ ] x_center_col / y_center_col 显式参数

## 版本历史

- v0（草稿）：`attract_logic` — 原型，rowwise 循环
- v1：`enrich_clusters` — 向量化，3 种模式
- v2：`enrich_clusters_pro` — 策略模式重构，5 种模式，3 种 centers 匹配
- v3（当前）：流式数据适配 — coord_type / 自动 radius / min_events / umap 警告

原始草稿：`projects/flow_cytometry/scripts/流式优化.r`
