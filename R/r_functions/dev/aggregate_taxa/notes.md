# aggregate_taxa — 开发笔记

## 当前状态

已完成开发、测试 (8/8 通过) 并归档。正式版本位于
`R/r_functions/lib/abundance_processing/aggregate_taxa.R`，
与同目录下的 `filter_abundance.R` / `norm_abundance.R` / `diff_abundance.R` 并列，
共同构成微生物组丰度处理流水线: 聚合 -> 过滤 -> 标准化 -> 差异分析。
本 dev 目录保留作为开发历史记录 (function.R/test.R/notes.md)，
后续修改请直接改 lib 版本，dev 版本不再同步维护。

## 来源

早期版本记录在 `自定义函数.md`（对话存档），当时基于 dplyr + tibble 实现
(rownames_to_column + left_join + group_by/summarise)。仓库里此前从未落地过
真正的 `.R` 文件——这次是照该文档的接口和行为重新实现。

## 与旧版 (md 存档) 的主要差异

1. **去掉 dplyr/tibble 依赖，改用纯 base R**：
   `abundance_processing/` 目录下的另外三个函数
   (`diff_abundance` / `filter_abundance` / `norm_abundance`) 都是纯 base R
   实现，没有引入 tidyverse 依赖。为保持目录风格一致，聚合逻辑改用
   `match()` 做 ID 对齐 + `rowsum()` 做分组求和，替代原来的
   `left_join()` + `group_by() %>% summarise()`。
2. **`rowsum()` 替代 group_by/summarise**：直接对矩阵按分组向量求和，
   比 dplyr 管道更快，也不需要中间的 rownames_to_column 转换。
3. 接口 (参数名、默认值、na_action/fill_na_with_rank/keep_all_abundance
   的语义) 与旧文档保持一致，未引入 breaking change。

## 参数速查

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `target_rank` | 必填 | 聚合目标分类等级，必须是 tax_table 的列名 |
| `na_action` | "unclassified" | "unclassified" 保留缺失值行 (重命名) / "remove" 直接丢弃 |
| `fill_na_with_rank` | NULL | 用上一级分类名智能填充缺失，如 "Unclassified_Firmicutes" |
| `keep_all_abundance` | TRUE | TRUE=保留丰度表全部特征 (无注释的记为 Unclassified)；FALSE=只保留交集 |

## 已知限制

- 要求 `abund_table` 与 `tax_table` 都有行名 (Feature ID)；两者均无法自动生成
  占位 ID（不同于 `diff_abundance` 等函数会自动补 `Feature_N`），因为聚合
  的核心就是靠 ID 匹配，自动补全的占位符无法与分类表对齐，容易掩盖真实的
  匹配失败问题，所以这里选择直接报错而不是静默生成。
- 未对超大特征数 (如百万级 ASV) 做专门性能优化，`rowsum()` 已经是
  base R 里效率较高的分组聚合方式，暂无需要进一步优化的证据。

## 测试覆盖

`test.R` 覆盖：基本聚合、total 丰度守恒、na_action 两种模式、
fill_na_with_rank 智能填充 (含双重缺失 Unassigned 场景)、
keep_all_abundance 两种模式、target_rank 不存在报错、
无匹配 ID 报错、与 filter_abundance/norm_abundance 的链式调用。

## 下一步

- [x] 归档为 `R/r_functions/lib/abundance_processing/aggregate_taxa.R`
      (归档版 sanity check 通过，与 dev 版行为一致)
- [ ] 更新 `filter_abundance.R`/`norm_abundance.R`/`diff_abundance.R` 的顶部
      注释或 README，补充 aggregate_taxa 在流水线中的位置说明（可选）
