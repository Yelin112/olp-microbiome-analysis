# ==============================================================================
# test_radial.R
# 测试 plot_enrich_radial()
# ==============================================================================

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../R/std_format.R")
source("../R/registry.R")
source("../R/styles/radial.R")

library(ggplot2)

# ------------------------------------------------------------------------------
# mock 数据
# ------------------------------------------------------------------------------

make_mock_cat <- function(n, cat_name, seed) {
  set.seed(seed)
  data.frame(
    ID = paste0(cat_name, "_", seq_len(n)),
    Description = paste0(cat_name, " pathway term ", seq_len(n)),
    geneID = sapply(seq_len(n), function(k) {
      paste(paste0("Gene", sample(1:30, sample(3:8, 1))), collapse = "/")
    }),
    pvalue = sort(runif(n, 0.0001, 0.04)),
    p.adjust = sort(runif(n, 0.001, 0.05)),
    Count = sample(5:25, n, replace = TRUE),
    BgRatio = paste0(sample(50:300, n, replace = TRUE), "/8000"),
    GeneRatio = paste0(
      sample(5:25, n, replace = TRUE),
      "/",
      sample(150:250, n, replace = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

mock_bp <- make_mock_cat(10, "BP", seed = 1)
mock_cc <- make_mock_cat(8, "CC", seed = 2)
mock_mf <- make_mock_cat(8, "MF", seed = 3)
mock_kegg <- make_mock_cat(9, "KEGG", seed = 4)

# ------------------------------------------------------------------------------
# 场景 1: 四分类（GO + KEGG），默认参数
# ------------------------------------------------------------------------------
cat("\n=== 场景1: GO+KEGG 四分类，默认参数 ===\n")

p1 <- plot_enrich_radial(
  list(BP = mock_bp, CC = mock_cc, MF = mock_mf, KEGG = mock_kegg),
  top_n = 8
)
ggsave(
  "figures/test_radial_s1_default.png",
  p1,
  width = 8,
  height = 8,
  dpi = 150
)
cat("[PASS] 场景1 → test_radial_s1_default.png\n")


# ------------------------------------------------------------------------------
# 场景 2: 两分类（微生物组场景）
# ------------------------------------------------------------------------------
cat("\n=== 场景2: 两分类（微生物组） ===\n")

mock_a <- make_mock_cat(10, "Bacteria", seed = 10)
mock_b <- make_mock_cat(8, "Archaea", seed = 11)

p2 <- plot_enrich_radial(
  list(Bacteria = mock_a, Archaea = mock_b),
  top_n = 8,
  category_colors = c(Bacteria = "#437f64", Archaea = "#d88c51"),
  category_colors_light = c(Bacteria = "#84ae9b", Archaea = "#e0b28e"),
  title = "Functional Enrichment by Domain"
)
ggsave(
  "figures/test_radial_s2_microbiome.png",
  p2,
  width = 8,
  height = 8,
  dpi = 150
)
cat("[PASS] 场景2 → test_radial_s2_microbiome.png\n")


# ------------------------------------------------------------------------------
# 场景 3: 调整 inner_radius 和 inset 范围
# ------------------------------------------------------------------------------
cat("\n=== 场景3: 调整内圈大小 ===\n")

p3 <- plot_enrich_radial(
  list(BP = mock_bp, CC = mock_cc, MF = mock_mf, KEGG = mock_kegg),
  top_n = 6,
  inner_radius = -300,
  inset_left = 0.15,
  inset_right = 0.85,
  inset_bottom = 0.15,
  inset_top = 0.85
)
ggsave(
  "figures/test_radial_s3_radius.png",
  p3,
  width = 8,
  height = 8,
  dpi = 150
)
cat("[PASS] 场景3 → test_radial_s3_radius.png\n")


# ------------------------------------------------------------------------------
# 场景 4: 通过 plot_enrich() 主函数调用
# ------------------------------------------------------------------------------
cat("\n=== 场景4: plot_enrich() 主函数 ===\n")

p4 <- plot_enrich(
  list(BP = mock_bp, CC = mock_cc, MF = mock_mf, KEGG = mock_kegg),
  type = "radial",
  top_n = 6
)
ggsave("figures/test_radial_s4_main.png", p4, width = 8, height = 8, dpi = 150)
cat("[PASS] 场景4 → test_radial_s4_main.png\n")

list_enrich_plots()

cat("\n所有场景完成，请检查 tests/figures/ 下的 PNG 文件。\n")
