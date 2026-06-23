# ==============================================================================
# test_point_bar.R
# 测试 plot_enrich_point_bar()
# ==============================================================================

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../R/std_format.R")
source("../R/registry.R")
source("../R/styles/point_bar.R")

library(ggplot2)

# ------------------------------------------------------------------------------
# mock 数据
# ------------------------------------------------------------------------------

make_mock_cat <- function(n, cat_name, seed) {
  set.seed(seed)
  data.frame(
    ID = paste0(cat_name, "_", seq_len(n)),
    Description = paste0(cat_name, " term ", seq_len(n)),
    geneID = sapply(seq_len(n), function(k) {
      paste(paste0("Gene", sample(1:30, sample(3:8, 1))), collapse = "/")
    }),
    pvalue = sort(runif(n, 0.0001, 0.04)),
    p.adjust = sort(runif(n, 0.001, 0.05)),
    Count = sample(3:15, n, replace = TRUE),
    BgRatio = paste0(sample(50:300, n, replace = TRUE), "/8000"),
    GeneRatio = paste0(
      sample(3:15, n, replace = TRUE),
      "/",
      sample(150:250, n, replace = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

mock_bp <- make_mock_cat(8, "BP", seed = 1)
mock_cc <- make_mock_cat(6, "CC", seed = 2)
mock_mf <- make_mock_cat(5, "MF", seed = 3)
mock_kegg <- make_mock_cat(7, "KEGG", seed = 4)

# ------------------------------------------------------------------------------
# 场景 1: 四分类（GO + KEGG），默认参数
# ------------------------------------------------------------------------------
cat("\n=== 场景1: GO+KEGG 四分类 ===\n")

p1 <- plot_enrich_point_bar(
  list(BP = mock_bp, CC = mock_cc, MF = mock_mf, KEGG = mock_kegg),
  top_n = 5
)
ggsave(
  "figures/test_pointbar_s1_default.png",
  p1,
  width = 14,
  height = 8,
  dpi = 150
)
cat("[PASS] 场景1 → test_pointbar_s1_default.png\n")


# ------------------------------------------------------------------------------
# 场景 2: 自定义配色 + 两列布局
# ------------------------------------------------------------------------------
cat("\n=== 场景2: 自定义配色 + 两列 facet ===\n")

p2 <- plot_enrich_point_bar(
  list(BP = mock_bp, CC = mock_cc, MF = mock_mf, KEGG = mock_kegg),
  top_n = 5,
  category_colors = c(
    BP = "#1cc7d0",
    CC = "#2dde98",
    MF = "#ff6c5f",
    KEGG = "#3369e7"
  ),
  facet_ncol = 2,
  title = "GO & KEGG Enrichment"
)
ggsave(
  "figures/test_pointbar_s2_custom.png",
  p2,
  width = 12,
  height = 10,
  dpi = 150
)
cat("[PASS] 场景2 → test_pointbar_s2_custom.png\n")


# ------------------------------------------------------------------------------
# 场景 3: 两分类（微生物组自定义场景）
# ------------------------------------------------------------------------------
cat("\n=== 场景3: 两分类（微生物组） ===\n")

mock_a <- make_mock_cat(8, "Bacteria", seed = 10)
mock_b <- make_mock_cat(6, "Archaea", seed = 11)

p3 <- plot_enrich_point_bar(
  list(Bacteria = mock_a, Archaea = mock_b),
  top_n = 6,
  category_colors = c(Bacteria = "#2196F3", Archaea = "#FF9800"),
  title = "Functional Enrichment by Domain"
)
ggsave(
  "figures/test_pointbar_s3_microbiome.png",
  p3,
  width = 10,
  height = 7,
  dpi = 150
)
cat("[PASS] 场景3 → test_pointbar_s3_microbiome.png\n")


# ------------------------------------------------------------------------------
# 场景 4: data.frame 输入 + ontology_col
# ------------------------------------------------------------------------------
cat("\n=== 场景4: data.frame + ontology_col ===\n")

combined_df <- rbind(
  cbind(mock_bp, category = "BP"),
  cbind(mock_cc, category = "CC"),
  cbind(mock_kegg, category = "KEGG")
)

p4 <- plot_enrich_point_bar(
  combined_df,
  ontology_col = "category",
  top_n = 5
)
ggsave("figures/test_pointbar_s4_df.png", p4, width = 12, height = 6, dpi = 150)
cat("[PASS] 场景4 → test_pointbar_s4_df.png\n")


# ------------------------------------------------------------------------------
# 场景 5: 调整样式参数
# ------------------------------------------------------------------------------
cat("\n=== 场景5: 样式参数调整 ===\n")

p5 <- plot_enrich_point_bar(
  list(BP = mock_bp, KEGG = mock_kegg),
  top_n = 6,
  bar_alpha = 0.8,
  point_size = 5,
  line_color = "#333333",
  line_width = 1.2,
  label_size = 3.5,
  x_buffer = 3
)
ggsave(
  "figures/test_pointbar_s5_style.png",
  p5,
  width = 10,
  height = 7,
  dpi = 150
)
cat("[PASS] 场景5 → test_pointbar_s5_style.png\n")


# ------------------------------------------------------------------------------
# 场景 6: 通过 plot_enrich() 主函数调用
# ------------------------------------------------------------------------------
cat("\n=== 场景6: plot_enrich() 主函数 ===\n")

p6 <- plot_enrich(
  list(BP = mock_bp, CC = mock_cc, MF = mock_mf, KEGG = mock_kegg),
  type = "point_bar",
  top_n = 4
)
ggsave(
  "figures/test_pointbar_s6_main.png",
  p6,
  width = 14,
  height = 6,
  dpi = 150
)
cat("[PASS] 场景6 → test_pointbar_s6_main.png\n")

list_enrich_plots()

cat("\n所有场景完成，请检查 tests/figures/ 下的 PNG 文件。\n")
