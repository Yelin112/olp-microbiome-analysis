# 测试文件: aggregate_taxa
# 创建日期: 2026-09-03

source("function.R")

cat("========================================\n")
cat("测试 aggregate_taxa()\n")
cat("========================================\n\n")

# 构造测试数据: 8 个 OTU, 4 个样本
abund <- matrix(
  c(
    10, 20, 5, 0, # OTU1
    8, 15, 3, 2, # OTU2
    0, 5, 10, 8, # OTU3
    3, 0, 2, 1, # OTU4 (Genus 缺失)
    12, 8, 6, 4, # OTU5
    1, 1, 1, 1, # OTU6 (Phylum 也缺失, 用于 Unassigned 场景)
    0, 0, 0, 0, # OTU7 (仅存在于丰度表, tax_table 中没有它)
    7, 7, 7, 7 # OTU8
  ),
  nrow = 8,
  byrow = TRUE,
  dimnames = list(paste0("OTU", 1:8), paste0("S", 1:4))
)

tax <- data.frame(
  Phylum = c(
    "Firmicutes",
    "Firmicutes",
    "Bacteroidetes",
    "Firmicutes",
    "Proteobacteria",
    NA,
    "Firmicutes"
  ),
  Genus = c("Lactobacillus", "Streptococcus", "Bacteroides", NA, "Ecoli", NA, "Clostridium"),
  row.names = paste0("OTU", c(1, 2, 3, 4, 5, 6, 8)) # OTU7 故意不在 tax_table 中
)

# === 测试1: 基本聚合 (Phylum, 默认 unclassified) ===
cat("测试1: 基本聚合到 Phylum\n")
res1 <- aggregate_taxa(abund, tax, target_rank = "Phylum", verbose = FALSE)
stopifnot(is.data.frame(res1))
stopifnot(sum(res1) == sum(abund)) # keep_all_abundance=TRUE 时总丰度不丢失
stopifnot("Unclassified" %in% rownames(res1)) # OTU6 (Phylum=NA) 应归为 Unclassified
stopifnot("Firmicutes" %in% rownames(res1))
cat("  OK - 总丰度守恒:", sum(res1) == sum(abund), "\n")
cat("  OK - 聚合后分类:", paste(rownames(res1), collapse = ", "), "\n\n")

# === 测试2: na_action = "remove" ===
cat("测试2: na_action='remove'\n")
res2 <- aggregate_taxa(
  abund,
  tax,
  target_rank = "Phylum",
  na_action = "remove",
  verbose = FALSE
)
stopifnot(!"Unclassified" %in% rownames(res2))
stopifnot(sum(res2) < sum(abund)) # OTU6 被丢弃, 总丰度应减少
cat("  OK - Unclassified 已移除, 总丰度:", sum(res2), "<", sum(abund), "\n\n")

# === 测试3: fill_na_with_rank 智能填充 ===
cat("测试3: fill_na_with_rank='Phylum' 时 Genus 缺失的智能填充\n")
res3 <- aggregate_taxa(
  abund,
  tax,
  target_rank = "Genus",
  fill_na_with_rank = "Phylum",
  verbose = FALSE
)
# OTU4: Genus=NA, Phylum=Firmicutes -> "Unclassified_Firmicutes"
stopifnot("Unclassified_Firmicutes" %in% rownames(res3))
# OTU6: Genus=NA, Phylum=NA -> "Unassigned"
stopifnot("Unassigned" %in% rownames(res3))
cat("  OK - 智能填充结果:", paste(rownames(res3), collapse = ", "), "\n\n")

# === 测试4: keep_all_abundance = FALSE (仅保留交集) ===
cat("测试4: keep_all_abundance=FALSE\n")
res4 <- suppressWarnings(aggregate_taxa(
  abund,
  tax,
  target_rank = "Phylum",
  keep_all_abundance = FALSE,
  verbose = FALSE
))
stopifnot(sum(res4) == sum(abund[rownames(abund) != "OTU7", ])) # OTU7 应被剔除
cat("  OK - OTU7 (无注释) 已按交集剔除, 总丰度:", sum(res4), "\n\n")

# === 测试5: keep_all_abundance = TRUE 时 OTU7 应计入 Unclassified ===
cat("测试5: keep_all_abundance=TRUE 时未注释特征归入 Unclassified\n")
res5 <- suppressWarnings(aggregate_taxa(
  abund,
  tax,
  target_rank = "Phylum",
  keep_all_abundance = TRUE,
  verbose = FALSE
))
stopifnot(sum(res5) == sum(abund)) # 一个特征都不能丢
cat("  OK - 总丰度守恒 (含未注释特征):", sum(res5) == sum(abund), "\n\n")

# === 测试6: target_rank 不存在应报错 ===
cat("测试6: target_rank 不存在时应报错\n")
err_caught <- tryCatch(
  {
    aggregate_taxa(abund, tax, target_rank = "NotAColumn", verbose = FALSE)
    FALSE
  },
  error = function(e) TRUE
)
stopifnot(err_caught)
cat("  OK - 正确报错\n\n")

# === 测试7: 无任何匹配 ID 时应报错 ===
cat("测试7: abund_table 与 tax_table 无匹配 ID 时应报错\n")
abund_no_match <- abund
rownames(abund_no_match) <- paste0("Gene", 1:8)
err_caught2 <- tryCatch(
  {
    aggregate_taxa(abund_no_match, tax, target_rank = "Phylum", verbose = FALSE)
    FALSE
  },
  error = function(e) TRUE
)
stopifnot(err_caught2)
cat("  OK - 正确报错\n\n")

# === 测试8: 与 filter_abundance / norm_abundance 的链式衔接 ===
cat("测试8: 与下游函数链式衔接 (filter_abundance -> norm_abundance)\n")
source("../../lib/abundance_processing/filter_abundance.R")
source("../../lib/abundance_processing/norm_abundance.R")

agg <- aggregate_taxa(abund, tax, target_rank = "Phylum", verbose = FALSE)
filt <- filter_abundance(agg, prev_threshold = 0, abund_threshold = 0, verbose = FALSE)
norm <- norm_abundance(filt$filtered, method = "TSS")
stopifnot(is.data.frame(norm))
stopifnot(all(abs(colSums(norm) - 1) < 1e-9))
cat("  OK - 链式调用成功, TSS 标准化后各样本列和为 1\n\n")

cat("========================================\n")
cat("所有测试通过!\n")
cat("========================================\n")
