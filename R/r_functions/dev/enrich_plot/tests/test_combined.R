# ==============================================================================
# test_combined.R
# 测试 plot_enrich_combined()
# ==============================================================================

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../R/std_format.R")
source("../R/registry.R")
source("../R/styles/combined.R")

library(ggplot2)

# ------------------------------------------------------------------------------
# 构建 mock 数据（模拟 GO + KEGG）
# ------------------------------------------------------------------------------

make_mock_cat <- function(n, cat_name, seed) {
  set.seed(seed)
  data.frame(
    ID = paste0(cat_name, "_", seq_len(n)),
    Description = paste0(cat_name, " term ", seq_len(n)),
    geneID = sapply(seq_len(n), function(k) {
      paste(paste0("Gene", sample(1:30, sample(3:10, 1))), collapse = "/")
    }),
    pvalue = sort(runif(n, 0.0001, 0.04)),
    p.adjust = sort(runif(n, 0.001, 0.10)),
    Count = sample(3:15, n, replace = TRUE),
    BgRatio = paste0(sample(50:300, n, replace = TRUE), "/8000"),
    stringsAsFactors = FALSE
  )
}

mock_bp <- make_mock_cat(8, "BP", seed = 1)
mock_cc <- make_mock_cat(6, "CC", seed = 2)
mock_mf <- make_mock_cat(5, "MF", seed = 3)
mock_kegg <- make_mock_cat(7, "KEGG", seed = 4)


# ==============================================================================
# 场景 1: named list 输入（核心场景）
# ==============================================================================
cat("\n=== 场景1: named list 输入 ===\n")

p1 <- plot_enrich_combined(
  list(BP = mock_bp, CC = mock_cc, MF = mock_mf, KEGG = mock_kegg),
  top_n = 5
)
ggsave(
  "figures/test_combined_s1_default.png",
  p1,
  width = 10,
  height = 8,
  dpi = 150
)
cat("[PASS] 场景1 默认参数 → test_combined_s1_default.png\n")


# ==============================================================================
# 场景 2: 自定义分类顺序和颜色
# ==============================================================================
cat("\n=== 场景2: 自定义顺序和配色 ===\n")

p2 <- plot_enrich_combined(
  list(BP = mock_bp, CC = mock_cc, MF = mock_mf, KEGG = mock_kegg),
  top_n = 4,
  category_order = c("KEGG", "MF", "CC", "BP"), # 自定义顺序（y轴从下到上）
  category_colors = c(
    BP = "#4393C3",
    CC = "#74C476",
    MF = "#FD8D3C",
    KEGG = "#D6604D"
  ),
  title = "GO & KEGG Pathway Enrichment"
)
ggsave(
  "figures/test_combined_s2_custom.png",
  p2,
  width = 10,
  height = 8,
  dpi = 150
)
cat("[PASS] 场景2 自定义 → test_combined_s2_custom.png\n")


# ==============================================================================
# 场景 3: data.frame 输入 + ontology_col
# ==============================================================================
cat("\n=== 场景3: data.frame 输入 + ontology_col ===\n")

# 合并 mock 数据，加分类列
combined_df <- rbind(
  cbind(mock_bp, category = "BP"),
  cbind(mock_cc, category = "CC"),
  cbind(mock_mf, category = "MF"),
  cbind(mock_kegg, category = "KEGG")
)

p3 <- plot_enrich_combined(
  combined_df,
  ontology_col = "category",
  top_n = 4
)
ggsave("figures/test_combined_s3_df.png", p3, width = 10, height = 8, dpi = 150)
cat("[PASS] 场景3 data.frame 输入 → test_combined_s3_df.png\n")


# ==============================================================================
# 场景 4: 关闭基因列表显示
# ==============================================================================
cat("\n=== 场景4: 不显示基因列表 ===\n")

p4 <- plot_enrich_combined(
  list(BP = mock_bp, CC = mock_cc, KEGG = mock_kegg),
  top_n = 5,
  show_gene_ids = FALSE
)
ggsave(
  "figures/test_combined_s4_nogenes.png",
  p4,
  width = 9,
  height = 7,
  dpi = 150
)
cat("[PASS] 场景4 无基因列表 → test_combined_s4_nogenes.png\n")


# ==============================================================================
# 场景 5: 两分类（非 GO/KEGG 场景，模拟微生物组自定义分类）
# ==============================================================================
cat("\n=== 场景5: 自定义两分类（微生物组场景） ===\n")

mock_a <- make_mock_cat(8, "Bacteria", seed = 10)
mock_b <- make_mock_cat(6, "Archaea", seed = 11)

p5 <- plot_enrich_combined(
  list(Bacteria = mock_a, Archaea = mock_b),
  top_n = 5,
  category_colors = c(Bacteria = "#2196F3", Archaea = "#FF9800"),
  title = "Functional Enrichment by Domain"
)
ggsave(
  "figures/test_combined_s5_microbiome.png",
  p5,
  width = 9,
  height = 7,
  dpi = 150
)
cat("[PASS] 场景5 微生物组场景 → test_combined_s5_microbiome.png\n")


# ==============================================================================
# 场景 6: 通过 plot_enrich() 主函数调用
# ==============================================================================
cat("\n=== 场景6: 通过 plot_enrich() 主函数 ===\n")

p6 <- plot_enrich(
  list(BP = mock_bp, CC = mock_cc, MF = mock_mf, KEGG = mock_kegg),
  type = "combined",
  top_n = 4
)
ggsave(
  "figures/test_combined_s6_main.png",
  p6,
  width = 10,
  height = 8,
  dpi = 150
)
cat("[PASS] 场景6 plot_enrich(type='combined') → test_combined_s6_main.png\n")

list_enrich_plots()

cat("\n所有场景完成，请检查生成的 PNG 文件。\n")
