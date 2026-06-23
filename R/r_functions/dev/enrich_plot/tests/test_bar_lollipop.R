# ==============================================================================
# test_bar_lollipop.R
# ==============================================================================

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../R/std_format.R")
source("../R/registry.R")
source("../R/styles/lollipop.R")
source("../R/styles/bar.R")

library(ggplot2)

# ------------------------------------------------------------------------------
# 构建 mock 数据（单组 + 双组）
# ------------------------------------------------------------------------------

make_mock <- function(n = 12, seed = 1) {
  set.seed(seed)
  data.frame(
    ID = paste0("hsa04", 100 + seq_len(n)),
    Description = paste0("Pathway ", LETTERS[seq_len(n)]),
    geneID = sapply(3:(n + 2), function(k) {
      paste(paste0("g", 1:k), collapse = "/")
    }),
    pvalue = sort(runif(n, 0.0001, 0.05)),
    p.adjust = sort(runif(n, 0.001, 0.15)),
    Count = sample(3:20, n, replace = TRUE),
    BgRatio = paste0(sample(50:200, n, replace = TRUE), "/8000"),
    stringsAsFactors = FALSE
  )
}

mock_up <- make_mock(10, seed = 1)
mock_down <- make_mock(10, seed = 2)
mock_down$Description <- paste0("Down_", mock_down$Description)

# ==============================================================================
# Lollipop 测试
# ==============================================================================
cat("\n=== Lollipop 测试 ===\n")

# 场景1: 单组 lollipop
p_lol_single <- plot_enrich_lollipop(mock_up, top_n = 8)
ggsave(
  "figures/test_lollipop_s1_single.png",
  p_lol_single,
  width = 8,
  height = 6,
  dpi = 150
)
cat("[PASS] lollipop 单组 → test_lollipop_s1_single.png\n")

# 场景2: 双向 lollipop（核心场景）
p_lol_dual <- plot_enrich_lollipop(
  list(Up = mock_up, Down = mock_down),
  top_n = 8,
  group_colors = c('#009dd3', '#f29c98'),
  size_var = "gene_count"
)
ggsave(
  "figures/test_lollipop_s2_dual.png",
  p_lol_dual,
  width = 9,
  height = 7,
  dpi = 150
)
cat("[PASS] lollipop 双向 → test_lollipop_s2_dual.png\n")

# 场景3: 双向 + 自定义标签 + 无大小映射
p_lol_dual2 <- plot_enrich_lollipop(
  list(Up = mock_up, Down = mock_down),
  top_n = 6,
  group_colors = c('#2ecc71', '#e74c3c'),
  group_labels = c("Enriched", "Depleted"),
  size_var = NULL,
  point_size = 5,
  title = "KEGG Pathway Comparison"
)
ggsave(
  "figures/test_lollipop_s3_custom.png",
  p_lol_dual2,
  width = 9,
  height = 6,
  dpi = 150
)
cat("[PASS] lollipop 自定义 → test_lollipop_s3_custom.png\n")

# 场景4: enrichResult 直接传入（如有 clusterProfiler）
if (
  requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)
) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  g <- bitr(
    c(
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
    ),
    "SYMBOL",
    "ENTREZID",
    org.Hs.eg.db
  )
  kegg <- enrichKEGG(g$ENTREZID, "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  if (nrow(kegg) >= 4) {
    p_lol_real <- plot_enrich_lollipop(
      kegg,
      top_n = 10,
      size_var = "gene_count"
    )
    ggsave(
      "figures/test_lollipop_s4_real.png",
      p_lol_real,
      width = 8,
      height = 6,
      dpi = 150
    )
    cat("[PASS] lollipop enrichResult → test_lollipop_s4_real.png\n")
  }
} else {
  cat("[跳过] enrichResult 场景需要 clusterProfiler\n")
}

# ==============================================================================
# Bar 测试
# ==============================================================================
cat("\n=== Bar 测试 ===\n")

# 场景1: 单组 bar
p_bar_single <- plot_enrich_bar(mock_up, top_n = 8)
ggsave(
  "figures/test_bar_s1_single.png",
  p_bar_single,
  width = 8,
  height = 6,
  dpi = 150
)
cat("[PASS] bar 单组 → test_bar_s1_single.png\n")

# 场景2: 双向 bar（核心场景）
p_bar_dual <- plot_enrich_bar(
  list(Up = mock_up, Down = mock_down),
  top_n = 8,
  group_colors = c('#009dd3', '#f29c98')
)
ggsave(
  "figures/test_bar_s2_dual.png",
  p_bar_dual,
  width = 9,
  height = 7,
  dpi = 150
)
cat("[PASS] bar 双向 → test_bar_s2_dual.png\n")

# 场景3: 自定义配色和标签
p_bar_dual2 <- plot_enrich_bar(
  list(Up = mock_up, Down = mock_down),
  top_n = 6,
  group_colors = c('#2ecc71', '#e74c3c'),
  group_labels = c("Enriched", "Depleted"),
  legend_title = "Direction",
  title = "Pathway Enrichment Comparison"
)
ggsave(
  "figures/test_bar_s3_custom.png",
  p_bar_dual2,
  width = 9,
  height = 6,
  dpi = 150
)
cat("[PASS] bar 自定义 → test_bar_s3_custom.png\n")

# ==============================================================================
# 通过 plot_enrich() 主函数调用
# ==============================================================================
cat("\n=== 通过 plot_enrich() 主函数 ===\n")

p_via_main_lol <- plot_enrich(mock_up, type = "lollipop", top_n = 8)
p_via_main_bar <- plot_enrich(
  list(Up = mock_up, Down = mock_down),
  type = "bar",
  top_n = 6
)
ggsave(
  "figures/test_main_lollipop.png",
  p_via_main_lol,
  width = 8,
  height = 6,
  dpi = 150
)
ggsave(
  "figures/test_main_bar_dual.png",
  p_via_main_bar,
  width = 9,
  height = 6,
  dpi = 150
)
cat("[PASS] plot_enrich(type='lollipop') → test_main_lollipop.png\n")
cat("[PASS] plot_enrich(type='bar', dual) → test_main_bar_dual.png\n")

list_enrich_plots()

cat("\n所有场景完成，请检查生成的 PNG 文件。\n")
