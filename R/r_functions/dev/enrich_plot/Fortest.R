
# 先测试核心层
setwd("e:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/r_functions/dev/enrich_plot")
source("tests/test_std_format.R")

# 再测试各样式（任意顺序）
source("tests/test_bar_lollipop.R")
source("tests/test_bubble.R")
source("tests/test_combined.R")
source("tests/test_point_bar.R")
source("tests/test_scatter.R")
source("tests/test_radial.R")

# 或用新的一键加载后测试
source("load_enrich_plot.R")
