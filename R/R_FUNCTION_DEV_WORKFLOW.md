# R 函数开发流程（dev → lib → archive）

> **最后更新**: 2026-09-19
> 本文档是通用流程说明。具体某个函数的设计决策、版本历史和待办，写在它自己的 `dev/<函数名>/notes.md` 里，不在这里。
> 旧版本文档（2025-04）是 `compare_plot_adaptive` 的一次性开发计划，已被本文取代。

---

## 1. 三个区域

```
R/r_functions/
├── dev/        # 开发区：正在写、正在改、还没定型的函数
├── lib/        # 稳定区：项目脚本 source 的目标，改动需谨慎
└── archive/    # 归档区：被替换掉的旧版本，只备份不使用
```

与 `R/utils/` 的分工：

| 位置 | 放什么 | 例子 |
|---|---|---|
| `R/utils/` | 与生物学无关、跨项目通用的绘图基础设施 | 配色、主题、面板固定、PPTX 导出 |
| `R/r_functions/lib/` | 领域分析/绘图函数（可能依赖 utils） | 丰度过滤、beta 多样性、网络分析、组间比较图 |

判断标准：函数里出现"物种/OTU/样本/丰度"等领域概念，放 `lib/`；只是在管 ggplot2 外观或导出，放 `utils/`。

---

## 2. dev/ 中的函数目录

每个函数一个文件夹，`dev/<function_name>/`：

```
dev/<function_name>/
├── function.R      # 函数本体（必需）
├── test.R          # 测试脚本（必需）
├── notes.md        # 开发笔记（必需）
├── example.R       # 使用示例（可选）
└── test_output/    # 测试产生的图（可选，PNG 不入库）
```

现有函数的实际写法：

- **单文件函数**：`aggregate_taxa`、`beta_calc`、`pcoa_plot`、`generate_replicates`、`enrich_clusters`
- **带子模块的函数**：`compare_plot_adaptive`（`strategies/` 下每个策略一个文件）、`enrich_plot`（`R/registry.R` + `R/styles/*.R`，入口是 `load_enrich_plot.R`）。函数变复杂、需要按"可插拔单元"拆分时采用这种结构。

### notes.md 建议包含

- 定位与目标（这个函数解决什么问题、不解决什么）
- 设计决策表（决策点 → 结论 → 理由）
- 函数架构图（主函数 + 内部辅助函数的调用关系）
- 与其他函数的对接方式（输入/输出数据结构）
- 版本历史、当前进度、已知问题

`beta_calc/notes.md` 和 `compare_plot_adaptive/notes.md` 是比较完整的范例。

---

## 3. 编码约定

### 3.1 签名与配色/主题

绘图函数必须包含 `palette` 和 `theme_use` 参数，配色通过 `palette_system.R` 获取，不在函数内定义调色板列表。详见 [utils/FUNCTION_TEMPLATE.R](utils/FUNCTION_TEMPLATE.R) 和 [utils/README.md](utils/README.md)。

```r
# 配色
scale_fill_pub_d(palette)            # 或 get_colors(palette, n)
# 主题
if (is.function(theme_use)) p <- p + theme_use()
```

### 3.1a 多机可用（本地 Windows + Linux 服务器）

- 函数文件里不写死绝对路径（`E:/`、`/home/`），不用 `sys.frame(1)$ofile`，不用 `../` 相对 `source()`。
- 工具层（`%||%`、`get_colors` 等）由调用方 `source(R/init.R)` + `load_utils()` 提供；函数文件里只保留依赖 `OLP_ROOT` 的兜底加载，写法见 [utils/FUNCTION_TEMPLATE.R](utils/FUNCTION_TEMPLATE.R) 第 0 节。
- 缺可选包（如 `officer`/`rvg`）时降级并提示，不中断；新依赖登记到 `R/check_deps.R`。

### 3.2 内部函数

私有辅助函数以 `.` 开头，并带函数名前缀避免污染全局命名空间，例如 `.norm_clr`、`.diff_wilcox`、`.grf_count`。只有主函数和确实需要被外部调用的函数不加点。

### 3.3 依赖处理

- 必需的包在文件头 `library()`；可选的重型包（`DESeq2`、`edgeR`、`metagenomeSeq`、`ANCOMBC` 等）在用到的分支里用 `requireNamespace()` 检查，缺失时给出清晰报错，不要让整个文件加载失败。
- 依赖 `utils/` 的函数，若调用方没有先 source，需要自行兜底加载（参考 `compare_plot_optimized.R` 文件头的 `if (!exists("get_colors"))` 写法）。

### 3.4 函数间对接

函数之间通过结构化返回值衔接，而不是让用户手动整理。例如 `beta_calc()` 返回 `coords`、`stat_text`、`axis_labs`，直接传给 `pcoa_plot()`。新增函数时先确定输入输出结构，再写实现。

### 3.5 注释

Roxygen2 风格（`#'`），至少写清 `@param`、`@return` 和一个可运行的 `@examples`。

---

## 4. 开发流程

### 各阶段的规范要求

| 阶段 | 位置 | 要求 |
|---|---|---|
| scratch | `R/scratch/` | 一次性探索，不强制；打算继续开发的，移到 dev 前先清掉 ERROR |
| dev | `dev/<name>/` | `check_conventions.R` 无 ERROR（路径类问题，换机器就坏） |
| lib | `lib/<category>/` | 无 ERROR、无 WARN（`--strict`）；确需保留的写法行末加 `# convention-ok` 并在 notes.md 说明 |

检查命令（只读，不改文件；默认检查 dev 和 lib，也可指定文件或目录）：

```bash
Rscript R/check_conventions.R R/r_functions/dev/<function_name>
Rscript R/check_conventions.R --strict R/r_functions/dev/<function_name>   # 晋升前
```

存量代码里有一批历史违规（如 `enrich_plot` 的样式文件、`abundance_processing` 的诊断图写死颜色），不要求一次改完，只要求新增和被修改的部分合规。

### 步骤 1：建目录

```
R/r_functions/dev/<function_name>/
```

创建 `function.R`、`test.R`、`notes.md`。在 `notes.md` 里先写清目标和设计决策，再动手写代码。

### 步骤 2：小步开发

- 一次只加一个功能，写完立刻在 `test.R` 里验证，不要一次性重写。
- 测试脚本用 `source(file.path(Sys.getenv("OLP_ROOT"), "R/init.R"))` 加载工具层和被测函数（dev 里的 `function.R` 用 `source(file.path(olp_path("R", "r_functions", "dev", "<name>"), "function.R"))`），不要写 `../` 相对路径或绝对路径，这样本地和服务器都能直接运行。
- 测试产出的图存到 `test_output/`。

### 步骤 3：晋升到 lib/

满足以下条件再晋升：

- [ ] `test.R` 全部通过
- [ ] `Rscript R/check_conventions.R --strict R/r_functions/dev/<function_name>` 无 ERROR、无 WARN
- [ ] 新增的依赖包已登记到 `R/check_deps.R`
- [ ] 参数校验和错误信息完整
- [ ] 至少在一个真实项目里用过
- [ ] `notes.md` 里没有阻塞性的已知问题

晋升操作：

1. 把 `function.R` 复制到 `lib/<category>/<function_name>.R`。`<category>` 按功能归类，已有：`abundance_processing`、`beta_calc`、`compare_plot`、`data_simulation`、`network_analysis`、`pcoa_plot`；没有合适分类就新建。
2. 需要给使用者看的说明，写成 `lib/<category>/<函数名>使用说明.md`（`beta_calc`、`pcoa_plot` 已有）。
3. **晋升后 dev 里不再保留 `function.R` 副本**，避免 dev 和 lib 两份代码各自演化。`notes.md` 和 `test.R` 可以留在 dev 作为开发记录。

### 步骤 4：修改已在 lib/ 的函数

1. 先把当前 lib 版本备份到 `archive/`，文件名加版本和日期：`<name>_v1_YYYYMMDD.R`。
2. 在 `dev/` 里改并测试（新建或复用同名 dev 文件夹）。
3. 通过后覆盖 `lib/` 中的文件。
4. 破坏性改动（改参数名、改返回结构）要检查 `projects/` 和 `collaborators/` 里有哪些脚本调用了它。

---

## 5. 在项目脚本中使用

分析脚本（`projects/<项目>/scripts/`）从 lib source：

```r
source(file.path(Sys.getenv("OLP_ROOT"), "R/init.R"))
load_utils(theme = TRUE)     # helpers + palette（+ theme；pptx = TRUE 加 PPTX 导出）
load_lib("compare_plot")     # lib 下的分类目录，或单个函数文件名（如 "beta_calc"）
```

`init.R` 负责定位仓库根目录（`OLP_ROOT` 环境变量 → git → 向上查找），脚本不再依赖所在目录深度，本地和 Linux 服务器通用。多机协作规则见 [CLAUDE.md](../CLAUDE.md) 的"多机协作"一节。新函数如果依赖 utils，靠调用方先 `load_utils()`，不要再自己拼路径 source。

**只 source `lib/`，不 source `dev/`。** 需要试用 dev 里的函数时，在 scratch 或 dev 自己的 `test.R` 里 source，验证后再晋升。

---

## 6. 当前状态与待整理事项

截至 2026-09-19：

| 函数 | dev | lib | 备注 |
|---|---|---|---|
| aggregate_taxa | ✓ | ✓ | 已晋升，dev 副本与 lib 一致，待清理 |
| beta_calc | ✓ | ✓ | 同上 |
| pcoa_plot | ✓ | ✓ | 同上 |
| generate_replicates | ✓ | ✓ | 同上 |
| compare_plot | — | ✓ | 已有 `strategy` 自适应参数 |
| compare_plot_adaptive | ✓ | — | 开发中（v3 管线式），与 `compare_plot` 的关系待决定：合并回去还是独立共存 |
| enrich_plot | ✓ | — | 开发中（注册表 + 多种样式） |
| enrich_clusters | ✓ | — | 开发中 |
| filter/norm/diff_abundance | — | ✓ | 从 `utils/` 迁入，无 dev 记录 |
| NetworkAnalyzer | — | ✓ | 从 `utils/` 迁入；`test_NetworkAnalyzer.R` 目前放在 lib 里，按约定应放 dev 或单独的测试目录 |
