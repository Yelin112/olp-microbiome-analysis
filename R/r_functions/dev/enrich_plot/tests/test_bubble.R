# ==============================================================================
# test_bubble.R
# 测试 plot_enrich_bubble()
# ==============================================================================

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../R/std_format.R")
source("../R/styles/bubble.R")

library(ggplot2)

# ------------------------------------------------------------------------------
# 场景 1: enrichResult 对象直接传入（最常用场景）
# ------------------------------------------------------------------------------
cat("\n=== 场景1: enrichResult 直接传入 ===\n")

if (
  requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)
) {
  library(clusterProfiler)
  library(org.Hs.eg.db)

  test_genes <- c(
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
    test_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  kegg <- enrichKEGG(
    gene = gene_ids$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  cat("KEGG 结果:", nrow(kegg), "个通路\n")

  # 默认参数
  p1 <- plot_enrich_bubble(kegg, top_n = 10)
  ggsave(
    "figures/test_bubble_s1_default.png",
    p1,
    width = 8,
    height = 6,
    dpi = 150
  )
  cat("[PASS] 场景1 默认参数 → test_bubble_s1_default.png\n")
} else {
  cat("[跳过] 需要 clusterProfiler + org.Hs.eg.db\n")
}


# ------------------------------------------------------------------------------
# 场景 2: mock data.frame（标准列名，无需 col_map）
# ------------------------------------------------------------------------------
cat("\n=== 场景2: data.frame 标准列名 ===\n")

mock_df <- data.frame(
  ID = paste0("hsa040", 60:74),
  Description = paste0("Pathway ", LETTERS[1:15]),
  geneID = sapply(3:17, function(n) paste(paste0("g", 1:n), collapse = "/")),
  pvalue = sort(runif(15, 0.0001, 0.05)),
  p.adjust = sort(runif(15, 0.001, 0.20)),
  Count = 3:17,
  BgRatio = paste0(sample(80:200, 15), "/8000"),
  stringsAsFactors = FALSE
)

p2 <- plot_enrich_bubble(mock_df, top_n = 10)
ggsave("figures/test_bubble_s2_mock.png", p2, width = 8, height = 6, dpi = 150)
cat("[PASS] 场景2 mock data.frame → test_bubble_s2_mock.png\n")


# ------------------------------------------------------------------------------
# 场景 3: x_var 切换
# ------------------------------------------------------------------------------
cat("\n=== 场景3: 切换 x_var ===\n")

p3a <- plot_enrich_bubble(mock_df, x_var = "rich_factor", top_n = 10)
p3b <- plot_enrich_bubble(mock_df, x_var = "gene_count", top_n = 10)

ggsave(
  "figures/test_bubble_s3_richfactor.png",
  p3a,
  width = 8,
  height = 6,
  dpi = 150
)
ggsave(
  "figures/test_bubble_s3_genecount.png",
  p3b,
  width = 8,
  height = 6,
  dpi = 150
)
cat("[PASS] 场景3 x_var 切换 → test_bubble_s3_*.png\n")


# ------------------------------------------------------------------------------
# 场景 4: 自定义配色 + 标题
# ------------------------------------------------------------------------------
cat("\n=== 场景4: 自定义配色与标题 ===\n")

p4 <- plot_enrich_bubble(
  mock_df,
  top_n = 10,
  color_low = "#d6eaf8",
  color_high = "#1a5276",
  size_range = c(3, 12),
  title = "KEGG Pathway Enrichment"
)
ggsave(
  "figures/test_bubble_s4_custom.png",
  p4,
  width = 8,
  height = 6,
  dpi = 150
)
cat("[PASS] 场景4 自定义配色 → test_bubble_s4_custom.png\n")


# ------------------------------------------------------------------------------
# 场景 5: 非标准 data.frame + col_map
# ------------------------------------------------------------------------------
cat("\n=== 场景5: 非标准 data.frame + col_map ===\n")

custom_df <- data.frame(
  ko_id = paste0("K000", 1:10),
  ko_name = paste0("KO function ", 1:10),
  gene_list = sapply(2:11, function(n) {
    paste(paste0("OTU", 1:n), collapse = "/")
  }),
  raw_p = sort(runif(10, 0.001, 0.05)),
  adj_p = sort(runif(10, 0.005, 0.20)),
  hits = 2:11,
  bg = sample(40:100, 10),
  stringsAsFactors = FALSE
)

p5 <- plot_enrich_bubble(
  custom_df,
  col_map = list(
    term_id = "ko_id",
    term_name = "ko_name",
    gene_ids = "gene_list",
    p_value = "raw_p",
    p_adjust = "adj_p",
    gene_count = "hits",
    bg_count = "bg"
  ),
  top_n = 8,
  x_var = "rich_factor",
  title = "KO Enrichment (Microbiome)"
)
ggsave(
  "figures/test_bubble_s5_colmap.png",
  p5,
  width = 8,
  height = 6,
  dpi = 150
)
cat("[PASS] 场景5 col_map → test_bubble_s5_colmap.png\n")


# ------------------------------------------------------------------------------
# 场景 6: filter_padj + terms 手动指定
# ------------------------------------------------------------------------------
cat("\n=== 场景6: filter_padj & 手动指定 terms ===\n")

p6 <- plot_enrich_bubble(
  mock_df,
  filter_padj = 0.10,
  terms = c("Pathway A", "Pathway B", "Pathway C", "Pathway D", "Pathway E")
)
ggsave("figures/test_bubble_s6_terms.png", p6, width = 8, height = 5, dpi = 150)
cat("[PASS] 场景6 filter_padj + terms → test_bubble_s6_terms.png\n")


cat("\n所有场景完成，请检查生成的 PNG 文件。\n")
