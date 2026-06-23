# ==============================================================================
# test_std_format.R
# 测试 as_enrich_std() 的三种输入场景
# ==============================================================================

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../R/std_format.R")

# ------------------------------------------------------------------------------
# 辅助：打印分隔线
# ------------------------------------------------------------------------------
.section <- function(title) {
  cat("\n", strrep("=", 60), "\n", title, "\n", strrep("=", 60), "\n", sep = "")
}
.pass <- function(msg) cat("  [PASS]", msg, "\n")
.fail <- function(msg) cat("  [FAIL]", msg, "\n")

.check <- function(expr, msg) {
  result <- tryCatch(expr, error = function(e) FALSE)
  if (isTRUE(result)) .pass(msg) else .fail(paste(msg, "(got:", result, ")"))
}


# ==============================================================================
# 场景 1: enrichResult 对象（clusterProfiler KEGG）
# ==============================================================================
.section("场景1: enrichResult 对象 (clusterProfiler KEGG)")

if (
  !requireNamespace("clusterProfiler", quietly = TRUE) ||
    !requireNamespace("org.Hs.eg.db", quietly = TRUE)
) {
  cat("  [跳过] 需要安装: clusterProfiler, org.Hs.eg.db\n")
} else {
  library(clusterProfiler)
  library(org.Hs.eg.db)

  # 准备测试基因
  test_genes_sym <- c(
    "NFKB1",
    "TNF",
    "IL1B",
    "IL6",
    "CXCL8",
    "PTGS2",
    "MMP9",
    "TLR4",
    "MAPK14",
    "ICAM1",
    "JUN",
    "STAT3",
    "MYD88",
    "CCL2",
    "IL10",
    "TGFB1",
    "VEGFA",
    "HIF1A",
    "NFKBIA",
    "RELA"
  )
  gene_ids <- bitr(
    test_genes_sym,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  cat("  基因ID转换:", nrow(gene_ids), "个\n")

  kegg_res <- enrichKEGG(
    gene = gene_ids$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  cat("  KEGG富集结果:", nrow(kegg_res), "个通路\n")

  # 转换
  std1 <- as_enrich_std(kegg_res)

  cat("\n--- 验证结果 ---\n")
  .check(inherits(std1, "EnrichStd"), "类标记为 EnrichStd")
  .check(inherits(std1, "data.frame"), "同时继承 data.frame")
  .check("term_id" %in% colnames(std1), "包含 term_id 列")
  .check("term_name" %in% colnames(std1), "包含 term_name 列")
  .check("gene_ids" %in% colnames(std1), "包含 gene_ids 列")
  .check("p_value" %in% colnames(std1), "包含 p_value 列")
  .check("p_adjust" %in% colnames(std1), "包含 p_adjust 列")
  .check("gene_count" %in% colnames(std1), "包含 gene_count 列")
  .check("gene_ratio" %in% colnames(std1), "包含 gene_ratio 列（衍生）")
  .check("rich_factor" %in% colnames(std1), "包含 rich_factor 列（衍生）")
  .check(is.integer(std1$gene_count), "gene_count 为整数类型")
  .check(
    all(std1$p_value >= 0 & std1$p_value <= 1, na.rm = TRUE),
    "p_value 值域 [0, 1]"
  )
  .check(all(std1$gene_ratio > 0, na.rm = TRUE), "gene_ratio 均为正数")
  .check(nrow(std1) == nrow(kegg_res), "行数与原始结果一致")

  cat("\n--- print() 输出 ---\n")
  print(std1)

  cat("\n--- summary() 输出 ---\n")
  summary(std1)
}


# ==============================================================================
# 场景 2: data.frame 输入（自动识别列名）
# ==============================================================================
.section("场景2: data.frame 输入（自动识别列名）")

# 模拟一个使用 clusterProfiler 标准列名的 data.frame
mock_df_standard <- data.frame(
  ID = c("hsa04064", "hsa04668", "hsa04620"),
  Description = c(
    "NF-kappa B signaling pathway",
    "TNF signaling pathway",
    "Toll-like receptor signaling pathway"
  ),
  geneID = c("3551/7124/3553", "7124/6885/1147", "7099/6654/3557"),
  pvalue = c(0.001, 0.005, 0.020),
  p.adjust = c(0.010, 0.025, 0.080),
  Count = c(3L, 3L, 3L),
  BgRatio = c("100/8000", "80/8000", "60/8000"),
  stringsAsFactors = FALSE
)

std2 <- as_enrich_std(mock_df_standard)

cat("\n--- 验证结果 ---\n")
.check(inherits(std2, "EnrichStd"), "类标记为 EnrichStd")
.check(
  std2$term_name[1] == "NF-kappa B signaling pathway",
  "term_name 正确映射"
)
.check(std2$gene_count[1] == 3L, "gene_count 正确")
.check(
  !is.na(std2$rich_factor[1]),
  "rich_factor 已计算（bg_count 来自 BgRatio）"
)

cat("\n预览:\n")
print(std2[, c(
  "term_id",
  "term_name",
  "p_adjust",
  "gene_count",
  "gene_ratio",
  "rich_factor"
)])


# ==============================================================================
# 场景 3: data.frame 输入（非标准列名 + col_map 指定）
# ==============================================================================
.section("场景3: data.frame 输入（非标准列名，使用 col_map）")

# 模拟微生物组分析工具的输出（列名不标准）
mock_df_custom <- data.frame(
  ko_id = c("K00001", "K00002", "K00003"),
  ko_name = c(
    "alcohol dehydrogenase",
    "alcohol dehydrogenase",
    "homoserine dehydrogenase"
  ),
  gene_list = c(
    "gene_A/gene_B/gene_C",
    "gene_D/gene_E",
    "gene_F/gene_G/gene_H"
  ),
  raw_pval = c(0.002, 0.008, 0.040),
  adj_pval = c(0.012, 0.032, 0.100),
  hit_count = c(3L, 2L, 3L),
  bg_genes = c(50L, 40L, 60L),
  extra_info = c("pathway_A", "pathway_B", "pathway_C"),
  stringsAsFactors = FALSE
)

std3 <- as_enrich_std(
  mock_df_custom,
  col_map = list(
    term_id = "ko_id",
    term_name = "ko_name",
    gene_ids = "gene_list",
    p_value = "raw_pval",
    p_adjust = "adj_pval",
    gene_count = "hit_count",
    bg_count = "bg_genes"
  )
)

cat("\n--- 验证结果 ---\n")
.check(inherits(std3, "EnrichStd"), "类标记为 EnrichStd")
.check(std3$term_id[1] == "K00001", "term_id 正确映射")
.check(std3$gene_count[2] == 2L, "gene_count 正确")
.check(!is.na(std3$rich_factor[1]), "rich_factor 已计算")
.check(
  round(std3$rich_factor[1], 4) == round(3 / 50, 4),
  "rich_factor 值正确 (3/50)"
)
.check("extra_info" %in% colnames(std3), "额外列 extra_info 被保留")
.check(std3$extra_info[1] == "pathway_A", "额外列值正确")

cat("\n预览:\n")
print(std3[, c(
  "term_id",
  "term_name",
  "p_adjust",
  "gene_count",
  "rich_factor",
  "extra_info"
)])


# ==============================================================================
# 场景 4: 错误处理测试
# ==============================================================================
.section("场景4: 错误处理")

bad_df <- data.frame(
  my_col = c("a", "b"),
  another = c(1, 2)
)

cat("\n--- 缺少必需列时应报错 ---\n")
tryCatch(
  as_enrich_std(bad_df),
  error = function(e) {
    .pass("正确抛出错误")
    cat("  错误信息:\n")
    cat(
      paste0(
        "  ",
        strsplit(conditionMessage(e), "\n")[[1]][1:3],
        collapse = "\n"
      ),
      "\n"
    )
  }
)

cat("\n--- col_map 指定不存在的列时应报错 ---\n")
tryCatch(
  as_enrich_std(
    mock_df_standard,
    col_map = list(term_name = "nonexistent_col")
  ),
  error = function(e) {
    .pass("正确抛出错误")
    cat("  错误信息:", conditionMessage(e), "\n")
  }
)


# ==============================================================================
# 场景 5: prep_enrich_data() 筛选函数
# ==============================================================================
.section("场景5: prep_enrich_data() 筛选和排序")

std_for_prep <- as_enrich_std(mock_df_standard)

# top_n
prepped_top2 <- prep_enrich_data(std_for_prep, top_n = 2)
.check(nrow(prepped_top2) == 2, "top_n = 2，结果为 2 行")

# filter_padj
prepped_sig <- prep_enrich_data(std_for_prep, filter_padj = 0.05)
.check(all(prepped_sig$p_adjust <= 0.05), "filter_padj = 0.05，所有结果符合")

# 手动指定 terms
prepped_terms <- prep_enrich_data(
  std_for_prep,
  terms = c("NF-kappa B signaling pathway", "TNF signaling pathway")
)
.check(nrow(prepped_terms) == 2, "手动指定 2 个 terms，结果为 2 行")

# order_by gene_count
prepped_order <- prep_enrich_data(
  std_for_prep,
  order_by = "gene_count",
  decreasing = TRUE
)
.check(
  prepped_order$gene_count[1] >= prepped_order$gene_count[nrow(prepped_order)],
  "order_by gene_count 降序正确"
)

cat("\n")


# ==============================================================================
# 场景 6: expand_gene_ids()
# ==============================================================================
.section("场景6: expand_gene_ids() 展开基因列")

expanded <- expand_gene_ids(std_for_prep)
.check(
  nrow(expanded) == sum(std_for_prep$gene_count),
  "展开后行数 = gene_count 之和"
)
.check(!any(grepl("/", expanded$gene_ids)), "展开后 gene_ids 不含 '/'")

cat("\n展开前:", nrow(std_for_prep), "行，展开后:", nrow(expanded), "行\n")
cat("前5行:\n")
print(expanded[1:5, c("term_name", "gene_ids", "p_adjust")], row.names = FALSE)

cat("\n", strrep("=", 60), "\n")
cat("全部测试完成\n")
cat(strrep("=", 60), "\n\n")
