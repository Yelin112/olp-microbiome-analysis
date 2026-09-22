# ============================================================================
# test_NetworkAnalyzer.R — NetworkAnalyzer 改造后功能测试
# ============================================================================
#
# 运行方式:
#   Rscript R/r_functions/lib/network_analysis/test_NetworkAnalyzer.R   # 在仓库根目录运行
#   （已设置 OLP_ROOT 环境变量时可在任意目录运行）
# 或在 RStudio 中 source 此文件

# ── 加载项目工具 ──────────────────────────────────────────────────────────
source(if (nzchar(Sys.getenv("OLP_ROOT"))) file.path(Sys.getenv("OLP_ROOT"), "R/init.R") else "R/init.R")
load_utils(theme = TRUE)
load_lib("network_analysis")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("  NetworkAnalyzer 改造测试\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# ============================================================================
# TEST 1: 从 abundance 矩阵构建相关网络
# ============================================================================
cat("── TEST 1: 构建相关网络 ──\n")

# 模拟 20 个 OTU x 30 个样本的 count 矩阵
set.seed(42)
otu_names <- paste0("OTU", sprintf("%03d", 1:20))
sample_names <- paste0("Sample", 1:30)
abundance <- matrix(
  rnbinom(20 * 30, mu = 100, size = 2),
  nrow = 20,
  ncol = 30,
  dimnames = list(otu_names, sample_names)
)

na <- NetworkAnalyzer$new(
  data = abundance,
  data_type = "abundance",
  filter_threshold = 0.001
)
cat("✓ 初始化成功\n")

na$build_network(
  method = "correlation",
  cor_method = "spearman",
  cor_threshold = 0.4,
  p_threshold = 0.05,
  use_abs = TRUE
)
cat("✓ 网络构建成功\n")

# ============================================================================
# TEST 2: 模块检测
# ============================================================================
cat("\n── TEST 2: 模块检测 ──\n")

na$detect_modules(method = "louvain")
cat("✓ 模块检测完成\n")

# ============================================================================
# TEST 3: 属性计算 + 节点角色
# ============================================================================
cat("\n── TEST 3: 网络属性 + 节点角色 ──\n")

na$calculate_properties(calculate_roles = TRUE)
cat("✓ 属性计算完成\n")

# 打印摘要
na

# 检查节点属性表
node_tbl <- na$get_node_table()
cat("\n节点属性表 (前 6 行):\n")
print(head(node_tbl))

# 验证 Zi-Pi 角色分类
if ("role" %in% names(node_tbl)) {
  cat("\n角色分布:\n")
  print(table(node_tbl$role))
}

# ============================================================================
# TEST 4: igraph 后端绘图（配色来自 palette_system）
# ============================================================================
cat("\n── TEST 4: igraph 绘图（palette 参数）──\n")

# 先确保有 module 属性
if (!"module" %in% vertex_attr_names(na$network)) {
  V(na$network)$module <- paste0(
    "M",
    sample(1:3, vcount(na$network), replace = TRUE)
  )
}

# 使用默认配色
tryCatch(
  {
    p_igraph <- na$plot_network(
      layout = "fr",
      color_by = "module",
      size_by = "degree",
      method = "igraph",
      palette = "NPG",
      main = "Network (igraph + NPG palette)"
    )
    cat("✓ igraph 绘图成功（NPG 配色）\n")
  },
  error = function(e) {
    cat("⚠ igraph 绘图失败:", e$message, "\n")
  }
)

# 测试直接用颜色向量
tryCatch(
  {
    p_igraph2 <- na$plot_network(
      layout = "fr",
      color_by = "module",
      size_by = "degree",
      method = "igraph",
      palette = c("#E64B35", "#4DBBD5", "#00A087"),
      main = "Network (igraph + custom colors)"
    )
    cat("✓ igraph 绘图成功（自定义颜色向量）\n")
  },
  error = function(e) {
    cat("⚠ igraph 绘图失败:", e$message, "\n")
  }
)

# ============================================================================
# TEST 5: ggraph 后端绘图
# ============================================================================
cat("\n── TEST 5: ggraph 绘图 ──\n")

has_ggraph <- requireNamespace("ggraph", quietly = TRUE)
if (has_ggraph) {
  tryCatch(
    {
      p_ggraph <- na$plot_network(
        method = "ggraph",
        palette = "NPG",
        theme_use = theme_nature
      )
      cat("✓ ggraph 绘图成功（NPG + theme_nature）\n")
    },
    error = function(e) {
      cat("⚠ ggraph 绘图失败:", e$message, "\n")
    }
  )

  # 测试默认主题回退
  tryCatch(
    {
      p_ggraph2 <- na$plot_network(method = "ggraph", palette = "D3")
      cat("✓ ggraph 绘图成功（默认 theme_pub_stat）\n")
    },
    error = function(e) {
      cat("⚠ ggraph 绘图失败:", e$message, "\n")
    }
  )
} else {
  cat("⏭ ggraph 未安装，跳过\n")
}

# ============================================================================
# TEST 6: visNetwork 后端绘图
# ============================================================================
cat("\n── TEST 6: visNetwork 绘图 ──\n")

has_visNetwork <- requireNamespace("visNetwork", quietly = TRUE)
if (has_visNetwork) {
  tryCatch(
    {
      p_vis <- na$plot_network(
        method = "visNetwork",
        palette = "JAMA"
      )
      cat("✓ visNetwork 绘图成功（JAMA 配色）\n")
    },
    error = function(e) {
      cat("⚠ visNetwork 绘图失败:", e$message, "\n")
    }
  )
} else {
  cat("⏭ visNetwork 未安装，跳过\n")
}

# ============================================================================
# TEST 7: 子网络提取
# ============================================================================
cat("\n── TEST 7: 子网络提取 ──\n")

sub_na <- na$subset_network(min_degree = 2)
cat("✓ 子网络提取成功\n")

# ============================================================================
# TEST 8: 网络比较
# ============================================================================
cat("\n── TEST 8: 网络比较 ──\n")

# 用不同阈值构建第二个网络
na2 <- NetworkAnalyzer$new(
  data = abundance,
  data_type = "abundance",
  filter_threshold = 0.001
)
na2$build_network(
  method = "correlation",
  cor_threshold = 0.3, # 更宽松的阈值
  p_threshold = 0.1,
  use_abs = TRUE
)
na2$detect_modules(method = "louvain")
na2$calculate_properties(calculate_roles = TRUE)

comparison <- na$compare_with(na2, permutations = 500)
cat("✓ 网络比较完成\n")
print(comparison)

# ============================================================================
# TEST 9: 对比图 S3 方法（验证配色/主题对接）
# ============================================================================
cat("\n── TEST 9: 对比图 S3 方法 ──\n")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  tryCatch(
    {
      p_comp_props <- plot.NetworkComparison(
        comparison,
        type = "properties",
        palette = "NPG",
        theme_use = theme_pub_stat
      )
      cat("✓ properties 对比图成功\n")
    },
    error = function(e) {
      cat("⚠ properties 对比图失败:", e$message, "\n")
    }
  )

  tryCatch(
    {
      p_comp_roles <- plot.NetworkComparison(
        comparison,
        type = "roles",
        palette = "JAMA",
        theme_use = theme_nature
      )
      cat("✓ roles 对比图成功（JAMA + theme_nature）\n")
    },
    error = function(e) {
      cat("⚠ roles 对比图失败:", e$message, "\n")
    }
  )

  tryCatch(
    {
      p_comp_deg <- plot.NetworkComparison(
        comparison,
        type = "degree",
        palette = "D3"
      )
      cat("✓ degree 对比图成功（默认 theme_pub_stat）\n")
    },
    error = function(e) {
      cat("⚠ degree 对比图失败:", e$message, "\n")
    }
  )
} else {
  cat("⏭ ggplot2 未安装，跳过\n")
}

# ============================================================================
# TEST 10: 导出功能
# ============================================================================
cat("\n── TEST 10: 导出功能 ──\n")

edge_tbl <- na$get_edge_table()
cat("✓ 边表导出:", nrow(edge_tbl), "行\n")

adj_mat <- na$get_adjacency_matrix(sparse = FALSE)
cat("✓ 邻接矩阵导出:", nrow(adj_mat), "×", ncol(adj_mat), "\n")

# ============================================================================
# TEST 11: 无 palette_system 时的降级测试
# ============================================================================
cat("\n── TEST 11: 降级行为验证 ──\n")

# 验证 get_colors() 存在
cat("  get_colors 可用:", exists("get_colors", mode = "function"), "\n")
cat(
  "  scale_fill_pub_d 可用:",
  exists("scale_fill_pub_d", mode = "function"),
  "\n"
)
cat("  theme_pub_stat 可用:", exists("theme_pub_stat", mode = "function"), "\n")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("  全部测试完成\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
