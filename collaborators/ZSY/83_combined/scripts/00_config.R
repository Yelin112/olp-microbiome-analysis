# ================================================================
# 00_config.R —— 全局参数 / 路径 / 包加载
# 项目：83 细胞系两批测序合并分析（LM_83 项目里的"83"细胞系 + 新的 A83 批次）
# 两批为同一套实验设计（Ctrl/NF/CAF/shART3 各3重复）的独立生物学批次
# （不同测序深度/建库批次），用 DESeq2 design = ~Batch + Condition 校正批次效应、
# 合并两批统计功效（每组 3+3=6 重复）。
# ================================================================

# ── 0. 包检查与加载 ──────────────────────────────────────────

required_pkgs <- c(
  "DESeq2",
  "ClusterGVis",
  "org.Hs.eg.db",
  "clusterProfiler",
  "enrichplot",
  "limma",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "ggsci",
  "ggrepel",
  "ComplexHeatmap",
  "grid",
  "writexl",
  "readxl",
  "Mfuzz",
  "factoextra",
  "pheatmap",
  "RColorBrewer"
)

not_installed <- required_pkgs[
  !sapply(required_pkgs, requireNamespace, quietly = TRUE)
]
if (length(not_installed) > 0) {
  message("请先安装以下包：", paste(not_installed, collapse = ", "))
  cat(sprintf(
    'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")\nBiocManager::install(c("%s"))\n',
    paste(not_installed, collapse = '", "')
  ))
  stop("缺少必要包，请安装后重新运行", call. = FALSE)
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ClusterGVis)
  library(org.Hs.eg.db)
  library(limma)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggsci)
  library(ComplexHeatmap)
  library(grid)
  library(writexl)
})

set.seed(42)
options(warn = 1)

# ── 1. 路径 ──────────────────────────────────────────────────

# 定位仓库根目录（本地/服务器通用）：优先 R_TOOLKIT_ROOT 环境变量，否则用 git
PROJECT_ROOT <- Sys.getenv("R_TOOLKIT_ROOT", unset = "")
if (!nzchar(PROJECT_ROOT)) PROJECT_ROOT <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE), error = function(e) "")
if (length(PROJECT_ROOT) != 1 || !nzchar(PROJECT_ROOT)) stop("找不到仓库根目录：请设置 R_TOOLKIT_ROOT 环境变量，或在仓库目录内运行")
source(file.path(PROJECT_ROOT, "init.R"))

# 本仓库（olp-microbiome-analysis）根目录：脚本本身在这个仓库里，用 git 定位
# （不能用 PROJECT_ROOT——它现在指向独立的 r-pub-toolkit 仓库，没有 collaborators/ 目录）
REPO_ROOT <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE), error = function(e) "")
if (length(REPO_ROOT) != 1 || !nzchar(REPO_ROOT)) stop("找不到本仓库根目录：请在仓库目录内运行脚本")

PROJ_DIR <- file.path(REPO_ROOT, "collaborators/ZSY/83_combined")
OUT_DIR <- file.path(PROJ_DIR, "outputs")

# 两个源项目的清洗后数据（01 脚本从这里读取，本项目不重复存 data/ 原始表）
LM83_PROJ_DIR <- file.path(REPO_ROOT, "collaborators/ZSY/LM_83")
A83_PROJ_DIR <- file.path(REPO_ROOT, "collaborators/ZSY/A83")

source(file.path(PROJECT_ROOT, "utils/helpers.R"))
source(file.path(PROJECT_ROOT, "utils/palette_system.R"))
source(file.path(PROJECT_ROOT, "utils/theme_system.R"))

# ── 2. 分组设计 ──────────────────────────────────────────────

# 批次标签：LM83 = 老项目里的"83"细胞系（TPM+gene_biotype注释，~15-18M reads/样本）
#           A83  = 新批次（FPKM+coding_type注释，~30-39M reads/样本）
BATCH_LEVELS <- c("LM83", "A83")

# 老批次 Condition 名 -> 统一命名（对齐新批次的命名）
CONDITION_MAP <- c(
  ctrl = "Ctrl",
  NFCM = "NF",
  CAFCM = "CAF",
  shCAFCM = "shART3"
)

GROUP_SETS <- list(
  full = c("Ctrl", "NF", "CAF", "shART3"),
  no_NF = c("Ctrl", "CAF", "shART3")
)

ACTIVE_VARIANT <- "full"

GROUP_COLORS <- setNames(get_colors("NPG", 4), GROUP_SETS$full)
BATCH_COLORS <- setNames(get_colors("Lancet", 2), BATCH_LEVELS)

# ── 3. 分析参数（集中管理，便于调整）────────────────────────

P <- list(
  min_count = 10, # 低表达过滤：单个样本中的最低 count
  min_samples = 6, # 至少有 N 个样本通过 min_count（24样本，按老流程3/12的比例放大到6/24）
  padj_cutoff = 0.05, # Wald 两两对比的差异显著性阈值
  lfc_cutoff = 1.0, # DEG 计数 / Volcano 用的 log2FC 阈值
  lrt_padj_cutoff = 0.05, # LRT 显著性阈值，决定进入聚类的基因全集
  min_std = 0.1, # clusterData() 内部低方差基因过滤阈值
  enrich_topn = 12, # 每个簇富集展示的 top term 数
  # elbow 图确定后手动填入，未填(NA)时聚类脚本会报错，
  # 避免在未经人工确认 K 值的情况下静默产出结果
  k_clusters = list(kmeans = 5, mfuzz = 5)
)

# ── 4. 输出目录 ──────────────────────────────────────────────

dir.create(file.path(OUT_DIR, "00_data"), recursive = TRUE, showWarnings = FALSE)
for (variant in names(GROUP_SETS)) {
  dir.create(
    file.path(OUT_DIR, variant, "DE"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  dir.create(
    file.path(OUT_DIR, variant, "ClusterGVis", "kmeans", "enrichment"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  dir.create(
    file.path(OUT_DIR, variant, "ClusterGVis", "mfuzz", "enrichment"),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ── 5. 通用绘图/保存辅助函数（沿用 LM_83 / A83 项目模式）─────

save_gg <- function(p, base_path, width = 10, height = 8) {
  ggsave(paste0(base_path, ".pdf"), plot = p, width = width, height = height)
  ggsave(
    paste0(base_path, ".png"),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
  message("  已保存: ", base_path, ".pdf / .png")
}

save_ht_pdf <- function(draw_func, filepath, width = 14, height = 12) {
  pdf(filepath, width = width, height = height)
  tryCatch(
    draw_func(),
    error = function(e) message("  [PDF 保存警告] ", conditionMessage(e))
  )
  dev.off()
  message("  已保存: ", filepath)
}

save_ht_png <- function(
  draw_func,
  filepath,
  width = 4200,
  height = 3600,
  res = 300
) {
  png(filepath, width = width, height = height, res = res)
  tryCatch(
    draw_func(),
    error = function(e) message("  [PNG 保存警告] ", conditionMessage(e))
  )
  dev.off()
  message("  已保存: ", filepath)
}

# 将 clusterData() 输出对象（ck/cm）截取为仅含指定簇号的子对象
subset_cluster_obj <- function(obj, selected) {
  wide <- obj$wide.res[obj$wide.res$cluster %in% selected, ]
  long <- obj$long.res[obj$long.res$cluster %in% selected, ]
  cluster_counts <- table(wide$cluster)
  long$cluster_name <- paste0(
    "cluster ",
    long$cluster,
    " (",
    cluster_counts[as.character(long$cluster)],
    ")"
  )
  out <- obj
  out$wide.res <- wide
  out$long.res <- long
  out$cluster.list <- obj$cluster.list[paste0("C", selected)]
  out
}

# ClusterGVis::enrichCluster() 内部用未校正的原始 pvalue < pvalueCutoff 筛选/排序
# term（qvalueCutoff/pvalueCutoff 传给 enrichGO/enrichKEGG 时都设成1，即先取全量
# 结果，再手动按原始 pvalue 过滤，p.adjust/qvalue 列在返回前被丢弃，不受多重检验
# 校正）。此函数是它的替代实现：改用 BH 校正后的 p.adjust 过滤/排序/取 topn，
# 输出结构与 enrichCluster() 完全一致（group/Description/pvalue位置放p.adjust值/
# ratio[/geneID]），因此可以直接替换 enrichCluster() 传给 visCluster(annoTermData=)，
# 无需改动可视化代码（visCluster 按列位置而非列名解析，goSize="pval"会自动读到
# 这里放的 p.adjust 值）。（与 LM_83 / A83 项目 00_config.R 中的同名函数保持一致）
enrichCluster_padj <- function(
  object,
  type = c("BP", "MF", "CC", "KEGG"),
  OrgDb = NULL,
  idTrans = TRUE,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  readable = TRUE,
  organism = "hsa",
  padjCutoff = 0.05,
  topn = 12,
  addGene = FALSE,
  useInternalData = FALSE
) {
  type <- match.arg(type)
  enrich.data <- object$wide.res
  purrr::map_df(seq_len(length(unique(enrich.data$cluster))), function(x) {
    cl_id <- unique(enrich.data$cluster)[x]
    tmp <- dplyr::filter(enrich.data, cluster == cl_id)
    if (idTrans) {
      gene.ent <- clusterProfiler::bitr(
        tmp$gene, fromType = fromType, toType = toType, OrgDb = OrgDb
      )
      target.gene <- unlist(gene.ent[, toType])
    } else {
      target.gene <- tmp$gene
    }

    if (type %in% c("BP", "MF", "CC")) {
      ego <- clusterProfiler::enrichGO(
        gene = target.gene, keyType = toType, OrgDb = OrgDb, ont = type,
        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = readable
      )
    } else {
      ego <- clusterProfiler::enrichKEGG(
        gene = target.gene, keyType = "kegg", organism = organism, universe = NULL,
        pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1,
        use_internal_data = useInternalData
      )
      if (readable) ego <- clusterProfiler::setReadable(ego, OrgDb = OrgDb, keyType = toType)
    }

    enrich_res <- data.frame(ego)
    if (nrow(enrich_res) == 0) return(NULL)

    df <- dplyr::arrange(
      dplyr::mutate(
        dplyr::filter(enrich_res, p.adjust < padjCutoff),
        group = paste0("C", cl_id)
      ),
      p.adjust
    )
    if (nrow(df) == 0) return(NULL)

    top <- if (length(topn) == 1) topn else topn[x]
    df <- dplyr::slice_head(df, n = top)
    df1 <- purrr::map_df(seq_len(nrow(df)), function(i) {
      tmp1 <- df[i, ]
      size <- unlist(strsplit(as.character(tmp1$GeneRatio), split = "/"))
      tmp1$ratio <- (as.numeric(size[1]) / as.numeric(size[2])) * 100
      tmp1
    })
    cols <- if (addGene) {
      c("group", "Description", "p.adjust", "ratio", "geneID")
    } else {
      c("group", "Description", "p.adjust", "ratio")
    }
    stats::na.omit(dplyr::select(df1, dplyr::all_of(cols)))
  })
}

message("00_config.R 加载完成 | ACTIVE_VARIANT = ", ACTIVE_VARIANT)
