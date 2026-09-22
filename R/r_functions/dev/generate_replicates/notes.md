# generate_replicates — 开发笔记

创建日期：2026-08-30
基于原始交接文档整理，已按 code review 结论修改。

## 安全警告（必读）

本函数从单个真实测量值生成噪声注入的伪重复，**不是真实的独立重复**。
如果这些伪重复被用来跑假设检验（t 检验等）或写进论文正式统计结论，等同于编造数据。

安全阀实现：
- 每次调用 `generate_replicates()` 都会 `warning()`，明确说明返回值是合成数据、不得用于假设检验。
- 返回值（向量/矩阵/data.frame）都带 attribute：
  - `attr(x, "synthetic") == TRUE`
  - `attr(x, "synthetic_type")` — 数据类型
  - `attr(x, "synthetic_warning")` — 完整警示文本
  - 下游代码/审查可以用这个 attribute 识别数据来源。已知局限：`cbind()`/`dplyr` 等后续处理可能丢失 attribute，这是 base R attribute 机制的固有限制，不是本函数能完全保证的。

## 与原始交接文档的差异（及原因）

1. **proportion 类型新增 `scale` 参数**（`"auto"` / `"0-1"` / `"0-100"`，默认 `"auto"`）。
   原文档的"value > 1 就当百分比"启发式无法区分"0~1 比例"和"0~100 百分比"在数值 ≤1
   时的情况（例如 0.5 可能是 50%，也可能是抑制率 0.5%）。默认行为保留原文档的
   auto 启发式（向后兼容），但触发时会 warning 提示这个局限，并建议用户显式传入
   `scale` 参数消歧义。

2. **count 类型不再强制 `round()` 原始值**。`rnbinom()` 的 `mu` 参数本身不要求整数，
   强制取整只会丢失精度，没有必要。

3. **`seed` 参数不再污染全局随机状态**。原文档未提及这一点，但直接 `set.seed()`
   会改变调用脚本后续所有随机行为。实现上用 `on.exit()` 保存/恢复调用前的
   `.Random.seed`。

4. **向量输入返回值类型统一为矩阵**，包括 `n = 1` 的情况（不会 drop 成向量），
   避免下游代码因为返回类型不稳定（有时向量有时矩阵）而出 bug。

5. **"warning 但不中断"的表述做了澄清**：多数情况 warning 不中断，唯一的例外是
   proportion 遇到 0 或 1（logit 无法计算）时 stop——这一点原文档表格标题和内容
   本身有点自相矛盾，这里以"设计决策记录"里的版本为准。

## 归档计划

- 开发目录：`R/r_functions/dev/generate_replicates/`（当前）
- 归档目标：`R/r_functions/lib/data_simulation/generate_replicates.R`
  （新建 `data_simulation` 分类，因为原文档提到后续可能扩展"多次运行 + bootstrap"
  的相关函数，适合放在同一个分类下）

## 后续扩展（暂不开发，记录备用）

某些分析工具每次只能跑出一个结果，后续可能需要支持"多次真实运行 + bootstrap 汇总"
的模式（这是获取真实重复的正确做法，优先于本工具的噪声注入方式）。
