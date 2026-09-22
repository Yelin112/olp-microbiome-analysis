# ============================================================================
# 使用示例：generate_replicates
#
# !! 再次提醒 !!
# 这些示例仅演示 API 用法。生成的伪重复不是真实重复，不能用于论文里的
# 假设检验或统计结论——详见 function.R 顶部的安全警告与 notes.md。
# ============================================================================

FUNC_DIR <- "e:/打工人/Zeng/2023-3 OLP Microbiome/分析/R/r_functions/dev/generate_replicates"
source(file.path(FUNC_DIR, "function.R"))

# ---- 示例 1：表达量（continuous）----
# 场景：qPCR 相对表达量，单次测量，想快速出一张探索性箱线图占位

tpm_value <- 128.4
tpm_reps <- generate_replicates(tpm_value, type = "continuous", n = 3, cv = 0.05, seed = 1)
print(tpm_reps)

# ---- 示例 2：抑制率 / 相对活性（proportion）----
# 场景：裂解酶抑制率，单次实验得到 68%
# 显式指定 scale 避免 0~1 / 0~100 歧义

inhibition_rate <- 68  # 百分比格式
inhibition_reps <- generate_replicates(
  inhibition_rate, type = "proportion", n = 3, scale = "0-100", seed = 1
)
print(inhibition_reps)

# ---- 示例 3：细菌相对丰度（abundance）----
# 场景：某属在宏基因组中的相对丰度向量，右偏、近 0

abundance_vec <- c(0.15, 0.003, 0.42, 0.001)
abundance_reps <- generate_replicates(
  abundance_vec, type = "abundance", n = 3, cv = 0.08, clip = TRUE, seed = 1
)
print(abundance_reps)

# ---- 示例 4：菌落计数（count）----
# 场景：CFU 计数，单次平板计数结果

cfu_count <- 152
cfu_reps <- generate_replicates(cfu_count, type = "count", n = 3, dispersion = 3, seed = 1)
print(cfu_reps)

# ---- 示例 5：data.frame 批量处理 ----

df <- data.frame(
  sample_id  = c("S1", "S2", "S3"),
  expression = c(10.2, 55.8, 3.1)
)
df_expanded <- generate_replicates(df, type = "continuous", n = 3, col = "expression", seed = 1)
print(df_expanded)

# ---- 示例 6：检查安全阀标记 ----
# 任何下游代码都可以通过 attr() 识别一个对象是否含合成数据

r <- generate_replicates(0.5, type = "abundance", n = 3, seed = 1)
cat("是否为合成数据:", attr(r, "synthetic"), "\n")
cat("类型:", attr(r, "synthetic_type"), "\n")
