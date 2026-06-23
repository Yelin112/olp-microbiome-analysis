# ============================================
# NetworkAnalyzer 测试和演示脚本
# ============================================

# 加载必要的包
if (!requireNamespace("R6", quietly = TRUE)) {
  install.packages("R6")
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(R6)
library(igraph)
library(ggplot2)

# 加载NetworkAnalyzer类
source("Scripts\\NetworkAnalyzer.R")

# ============================================
# 测试1: 从模拟丰度数据构建网络
# ============================================

cat("=== 测试1: 从丰度矩阵构建网络 ===\n")

# 模拟微生物丰度数据
set.seed(123)
n_features <- 50 # 50个OTU/ASV
n_samples <- 30 # 30个样本

# 创建丰度矩阵 (行=特征, 列=样本)
abundance <- matrix(
  rpois(n_features * n_samples, lambda = 20),
  nrow = n_features,
  ncol = n_samples
)

# 添加一些特征间的相关性结构
for (i in seq(1, n_features, by = 5)) {
  indices <- i:min(i + 4, n_features)
  # 同一组内的特征有相似的丰度模式
  base_pattern <- rnorm(n_samples, mean = 30, sd = 5)
  for (j in indices) {
    abundance[j, ] <- abundance[j, ] + base_pattern + rnorm(n_samples, sd = 5)
  }
}

# 命名
rownames(abundance) <- paste0("OTU_", 1:n_features)
colnames(abundance) <- paste0("Sample_", 1:n_samples)

# 创建元数据
metadata <- data.frame(
  Sample = colnames(abundance),
  Group = rep(c("Control", "Treatment"), each = n_samples / 2),
  row.names = colnames(abundance)
)

# 创建网络分析对象
net <- NetworkAnalyzer$new(
  data = abundance,
  metadata = metadata,
  data_type = "abundance",
  filter_threshold = 0.001
)

# 构建网络
net$build_network(
  method = "correlation",
  cor_method = "spearman",
  cor_threshold = 0.6,
  p_threshold = 0.05,
  use_abs = FALSE # 保留正负相关
)

# 检测模块
net$detect_modules(method = "louvain")

# 计算网络属性
net$calculate_properties(calculate_roles = TRUE)

# 查看网络信息
print(net)

# 查看节点表
cat("\n节点属性表 (前10行):\n")
nodes <- net$get_node_table()
print(head(nodes, 10))

# 查看节点角色分布
cat("\n节点角色分布:\n")
print(table(nodes$role))

# ============================================
# 测试2: 网络可视化
# ============================================

cat("\n=== 测试2: 网络可视化 ===\n")

# 基础igraph可视化
par(mfrow = c(1, 2))

# 使用Fruchterman-Reingold布局
net$plot_network(
  layout = "fr",
  color_by = "module",
  size_by = "degree"
)
title("Fruchterman-Reingold Layout")

# 使用Kamada-Kawai布局
net$plot_network(
  layout = "kk",
  color_by = "module",
  size_by = "degree"
)
title("Kamada-Kawai Layout")
par(mfrow = c(1, 1))

# ============================================
# 测试3: 分组网络比较
# ============================================

cat("\n=== 测试3: 分组网络比较 ===\n")

# 分别构建两组的网络
abundance_control <- abundance[, metadata$Group == "Control"]
abundance_treatment <- abundance[, metadata$Group == "Treatment"]

# 控制组网络
net_control <- NetworkAnalyzer$new(
  data = abundance_control,
  data_type = "abundance",
  filter_threshold = 0.001
)
net_control$build_network(
  method = "correlation",
  cor_method = "spearman",
  cor_threshold = 0.6,
  p_threshold = 0.05
)
net_control$detect_modules(method = "louvain")
net_control$calculate_properties(calculate_roles = TRUE)

cat("\n控制组网络:\n")
print(net_control)

# 处理组网络
net_treatment <- NetworkAnalyzer$new(
  data = abundance_treatment,
  data_type = "abundance",
  filter_threshold = 0.001
)
net_treatment$build_network(
  method = "correlation",
  cor_method = "spearman",
  cor_threshold = 0.6,
  p_threshold = 0.05
)
net_treatment$detect_modules(method = "louvain")
net_treatment$calculate_properties(calculate_roles = TRUE)

cat("\n处理组网络:\n")
print(net_treatment)

# 比较两个网络
cat("\n开始网络比较...\n")
comparison <- net_control$compare_with(net_treatment, permutations = 100)

cat("\n网络比较结果:\n")
print(comparison)

# 可视化比较结果
if (requireNamespace("ggplot2", quietly = TRUE)) {
  cat("\n生成比较可视化...\n")

  # 属性比较
  p1 <- plot(comparison, type = "properties")
  print(p1)

  # 角色比较
  if (!is.null(comparison$node_roles)) {
    p2 <- plot(comparison, type = "roles")
    print(p2)
  }

  # 度分布比较
  p3 <- plot(comparison, type = "degree")
  print(p3)
}

# ============================================
# 测试4: 子网络提取
# ============================================

cat("\n=== 测试4: 子网络提取 ===\n")

# 按模块提取
modules <- unique(V(net$network)$module)
if (length(modules) > 0) {
  subnet_m1 <- net$subset_network(
    module = modules[1],
    remove_isolates = TRUE
  )
  cat(sprintf("\n提取模块 %s 的子网络:\n", modules[1]))
  print(subnet_m1)
}

# 提取高度节点
subnet_hubs <- net$subset_network(
  min_degree = 5,
  remove_isolates = TRUE
)
cat("\n提取高度节点 (degree >= 5) 的子网络:\n")
print(subnet_hubs)

# ============================================
# 测试5: 特定节点分析 (模拟Rothia)
# ============================================

cat("\n=== 测试5: 特定节点分析 ===\n")

# 选择一个节点作为目标 (模拟 Rothia_mucilaginosa)
target_node <- rownames(abundance)[1]
cat(sprintf("\n分析目标节点: %s\n", target_node))

# 在控制组中的属性
nodes_control <- net_control$get_node_table()
if (target_node %in% nodes_control$node) {
  target_control <- nodes_control[nodes_control$node == target_node, ]
  cat("\n控制组中的属性:\n")
  print(target_control)
}

# 在处理组中的属性
nodes_treatment <- net_treatment$get_node_table()
if (target_node %in% nodes_treatment$node) {
  target_treatment <- nodes_treatment[nodes_treatment$node == target_node, ]
  cat("\n处理组中的属性:\n")
  print(target_treatment)
}

# 比较该节点在两组中的度
if (
  target_node %in% nodes_control$node && target_node %in% nodes_treatment$node
) {
  cat("\n节点度比较:\n")
  cat(sprintf("  控制组: %d\n", target_control$degree))
  cat(sprintf("  处理组: %d\n", target_treatment$degree))
  cat(sprintf("  差异: %d\n", target_treatment$degree - target_control$degree))
}

# 提取包含该节点的子网络
neighbors_idx <- neighbors(
  net$network,
  which(V(net$network)$name == target_node)
)
neighbors_names <- V(net$network)$name[neighbors_idx]

cat(sprintf("\n%s 的邻居节点 (%d个):\n", target_node, length(neighbors_names)))
print(neighbors_names)

# 创建ego网络
ego_nodes <- c(target_node, neighbors_names)
ego_net <- net$subset_network(nodes = ego_nodes, remove_isolates = FALSE)

cat(sprintf("\n%s 的ego网络:\n", target_node))
print(ego_net)

# 可视化ego网络
ego_net$plot_network(
  layout = "fr", # 使用Fruchterman-Reingold布局
  color_by = "module"
)
title(paste("Ego network of", target_node))

# ============================================
# 测试6: 导出结果
# ============================================

cat("\n=== 测试6: 导出结果 ===\n")

# 导出节点表
node_table <- net$get_node_table()
write.csv(node_table, "test_node_attributes.csv", row.names = FALSE)
cat("节点表已导出到: test_node_attributes.csv\n")

# 导出边表
edge_table <- net$get_edge_table()
write.csv(edge_table, "test_edge_list.csv", row.names = FALSE)
cat("边表已导出到: test_edge_list.csv\n")

# 保存网络文件
net$save_network("test_network.graphml", format = "graphml")
cat("网络已保存到: test_network.graphml\n")

# ============================================
# 测试7: Zi-Pi 图
# ============================================

cat("\n=== 测试7: 节点角色可视化 (Zi-Pi图) ===\n")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  nodes <- net$get_node_table()

  p <- ggplot(nodes, aes(x = Pi, y = Zi, color = role, size = degree)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = 2.5, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0.62, linetype = "dashed", color = "grey40") +
    scale_color_manual(
      values = c(
        "Peripherals" = "#4DAF4A",
        "Connectors" = "#377EB8",
        "Module hubs" = "#E41A1C",
        "Network hubs" = "#984EA3"
      )
    ) +
    labs(
      title = "Node Roles (Zi-Pi Plot)",
      subtitle = paste(nrow(nodes), "nodes"),
      x = "Participation coefficient (Pi)",
      y = "Within-module degree (Zi)",
      color = "Role",
      size = "Degree"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    )

  print(p)
}

# ============================================
# 测试8: 不同输入格式
# ============================================

cat("\n=== 测试8: 不同输入格式测试 ===\n")

# 从邻接矩阵
adj_matrix <- as.matrix(net$get_adjacency_matrix())
net_from_adj <- NetworkAnalyzer$new(
  data = adj_matrix,
  data_type = "adjacency"
)
cat("\n从邻接矩阵创建网络:\n")
print(net_from_adj)

# 从边列表
edge_list <- net$get_edge_table()
net_from_edge <- NetworkAnalyzer$new(
  data = edge_list,
  data_type = "edgelist"
)
cat("\n从边列表创建网络:\n")
print(net_from_edge)

# 从igraph对象
igraph_obj <- net$network
net_from_igraph <- NetworkAnalyzer$new(
  data = igraph_obj,
  data_type = "igraph"
)
cat("\n从igraph对象创建网络:\n")
print(net_from_igraph)

# ============================================
# 完成
# ============================================

cat("\n")
cat("==========================================\n")
cat("所有测试完成！\n")
cat("==========================================\n")
cat("\n生成的文件:\n")
cat("  - test_node_attributes.csv\n")
cat("  - test_edge_list.csv\n")
cat("  - test_network.graphml\n")
cat("\n")

# 清理
cat("是否要清理测试文件? (手动删除以下文件)\n")
cat("  rm test_*.csv test_*.graphml\n")
