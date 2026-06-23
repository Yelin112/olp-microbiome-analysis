# 开发笔记：beta_calc

## 项目信息
- **开始日期**：2026-04-06
- **开发目录**：`r_functions/dev/beta_calc/`
- **定位**：通用多元距离计算与排序，不限于微生物组

---

## 设计决策（讨论记录）

| 决策点 | 结论 |
|--------|------|
| 数据通用性 | 通用"样本×特征"矩阵为主，phyloseq/OTU 为可选扩展 |
| 输入方向 | 样本为行（与 vegan 一致），OTU 表自动转置 |
| 降维方法 | PCoA、NMDS、PCA 三种 |
| PCA 统计 | 自动用欧氏距离做 Adonis（与 PCA 等价距离一致）|
| 两两比较 | 默认关（`pairwise = FALSE`） |
| 输出对接 | `res$coords` + `res$stat_text` + `res$axis_labs` 直接传给 `pcoa_plot()` |

---

## 函数架构

```
beta_calc()
├── .resolve_input()   解析三种输入方式，对齐样本
├── .calc_dist()       计算距离矩阵（vegan 或 phyloseq）
├── .do_ordination()   PCoA / NMDS / PCA 排序
├── .run_one_stat()    单次统计检验（adonis2 / anosim / mrpp）
└── .format_stat_text() 格式化注释字符串
```

---

## 输入方式优先级

```
ps（phyloseq）> otu（特征×样本）> data（样本×特征）
```

样本对齐逻辑：
- 有行名时：取 `intersect(rownames(data), rownames(meta))`
- 无行名时：假设行顺序一致（给出 message 提示）

---

## 方法-距离对应关系

| method | 距离处理 | 统计用距离 |
|--------|---------|-----------|
| PCoA | `vegdist()` → `cmdscale()` | 同上 |
| NMDS | `vegdist()` → `metaMDS()` | 同上 |
| PCA  | 忽略 dist，直接 `prcomp()` | 欧氏距离（自动）|

PCA + 欧氏距离在数学上等价，统计结果解释合理。

---

## 输出字段说明

| 字段 | 类型 | 说明 |
|------|------|------|
| `coords` | data.frame | 降维坐标 + 分组列，直接传给 `pcoa_plot(data=)` |
| `dist_mat` | dist | 距离矩阵，可复用于其他分析 |
| `stat_table` | data.frame | 整体检验结果，4列：method/statistic_name/statistic/p_value |
| `pair_table` | data.frame | 两两比较，6列：group1/group2/method/statistic_name/statistic/p_value |
| `stat_text` | character | 格式化字符串，直接传给 `pcoa_plot(stat_text=)` |
| `axis_labs` | character | 格式化轴标签，如 "PCoA1 (39.9%)" |
| `stress` | numeric/NULL | NMDS stress，已包含在 stat_text 第一行 |

---

## pcoa_plot() 对接

```r
res <- beta_calc(data = otu_relabund, meta = meta,
                 group = "Group", dist = "bray")

pcoa_plot(res$coords,
          x        = "PC1",
          y        = "PC2",
          group.by = "Group",
          xlab     = res$axis_labs[1],
          ylab     = res$axis_labs[2],
          stat_text = res$stat_text,
          marginal  = "density")
```

---

## 依赖包
- `vegan`（必需：vegdist, adonis2, anosim, mrpp, metaMDS）
- `phyloseq`（可选：unifrac/wunifrac 距离，或 ps 输入）

---

## 当前进度
- [x] function.R — 主函数完成
- [x] test.R — 测试套件（Part 1-9，含对接测试）
- [ ] 运行测试，记录结果
- [ ] 归档至 lib/

---

## 已知限制与待处理
- [ ] `adonis2` 在两两比较时置换次数仅 999，极小样本可能不稳定
- [ ] phyloseq 路径下 UniFrac 未经实际数据测试（需进化树）
- [ ] NMDS 多次运行结果可能因 stress 略有浮动（seed 参数可控制）
