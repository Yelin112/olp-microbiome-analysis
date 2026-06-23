# ==============================================================================
# load_enrich_plot.R: 一键加载富集分析可视化系统
# 用法: source("path/to/enrich_plot/load_enrich_plot.R")
# ==============================================================================

.ep_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)

source(file.path(.ep_dir, "R/std_format.R"))
source(file.path(.ep_dir, "R/registry.R"))
for (.f in list.files(file.path(.ep_dir, "R/styles"), pattern = "\\.R$", full.names = TRUE)) {
  source(.f)
}

rm(.ep_dir, .f)
cat("enrich_plot 系统已加载，可用类型:\n")
list_enrich_plots()
