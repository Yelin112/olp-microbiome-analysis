# ==============================================================================
# test_scatter.R
# 测试 plot_enrich_scatter()
# ==============================================================================

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../R/std_format.R")
source("../R/registry.R")
source("../R/styles/scatter.R")

library(ggplot2)

# ------------------------------------------------------------------------------
# mock 数据
# ------------------------------------------------------------------------------

make_mock <- function(n, seed = 1) {
  set.seed(seed)
  data.frame(
    ID = paste0("hsa04", 100 + seq_len(n)),
    Description = paste0("Pathway ", LETTERS[seq_len(n)]),
    geneID = sapply(seq_len(n), function(k) {
      paste(paste0("g", sample(1:30, sample(3:10, 1))), collapse = "/")
    }),
    pvalue = sort(runif(n, 1e-5, 0.04)),
    p.adjust = sort(runif(n, 1e-4, 0.05)),
    Count = sample(5:35, n, replace = TRUE),
    BgRatio = paste0(sample(50:300, n, replace = TRUE), "/8000"),
    stringsAsFactors = FALSE
  )
}

mock_up <- make_mock(15, seed = 1)
mock_down <- make_mock(15, seed = 2)
mock_down$Description <- paste0("Down_", mock_down$Description)

# ------------------------------------------------------------------------------
# 场景 1: 默认参数
# ------------------------------------------------------------------------------
cat("\n=== 场景1: 默认参数 ===\n")

p1 <- plot_enrich_scatter(
  list(Up = mock_up, Down = mock_down),
  top_n = 12
)
ggsave(
  "figures/test_scatter_s1_default.png",
  p1,
  width = 9,
  height = 7,
  dpi = 150
)
cat("[PASS] 场景1 → test_scatter_s1_default.png\n")


# ------------------------------------------------------------------------------
# 场景 2: 自定义配色与标签
# ------------------------------------------------------------------------------
cat("\n=== 场景2: 自定义配色 ===\n")

p2 <- plot_enrich_scatter(
  list(Up = mock_up, Down = mock_down),
  top_n = 12,
  group_colors = c(Up = "#009dd3", Down = "#f29c98"),
  group_labels = c(Up = "Enriched", Down = "Depleted"),
  legend_title = "Direction",
  bg_color = "#fafaf5",
  title = "KEGG Pathway Enrichment"
)
ggsave(
  "figures/test_scatter_s2_custom.png",
  p2,
  width = 9,
  height = 7,
  dpi = 150
)
cat("[PASS] 场景2 → test_scatter_s2_custom.png\n")


# ------------------------------------------------------------------------------
# 场景 3: 调整 x 范围和 y 刻度
# ------------------------------------------------------------------------------
cat("\n=== 场景3: 自定义坐标范围 ===\n")

p3 <- plot_enrich_scatter(
  list(Up = mock_up, Down = mock_down),
  top_n = 10,
  x_limits = c(-40, 40),
  y_tick_step = 2,
  size_range = c(2, 10),
  label_size = 2.5
)
ggsave("figures/test_scatter_s3_axis.png", p3, width = 9, height = 7, dpi = 150)
cat("[PASS] 场景3 → test_scatter_s3_axis.png\n")


# ------------------------------------------------------------------------------
# 场景 4: 通过 plot_enrich() 主函数调用
# ------------------------------------------------------------------------------
cat("\n=== 场景4: plot_enrich() 主函数 ===\n")

p4 <- plot_enrich(
  list(Up = mock_up, Down = mock_down),
  type = "scatter",
  top_n = 10
)
ggsave("figures/test_scatter_s4_main.png", p4, width = 9, height = 7, dpi = 150)
cat("[PASS] 场景4 → test_scatter_s4_main.png\n")

list_enrich_plots()

cat("\n所有场景完成，请检查 tests/figures/ 下的 PNG 文件。\n")
