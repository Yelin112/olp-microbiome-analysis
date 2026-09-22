# ================================================================
# 00_config.R —— 全局参数 / 路径 / 包加载
# 项目：LM_83 CAFCM 上清诱导转录重编程 + shCAFCM 救援效应
# 细胞系：83 / LM  |  分组：ctrl / NFCM / CAFCM / shCAFCM（各3重复）
# ================================================================

# ── 0. 包检查与加载 ──────────────────────────────────────────

required_pkgs <- c(
  "DESeq2",
  "ClusterGVis",
  "org.Hs.eg.db",
  "clusterProfiler",
  "enrichplot",
  "data.table",
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

PROJ_DIR <- file.path(REPO_ROOT, "collaborators/ZSY/LM_83")
DATA_DIR <- file.path(PROJ_DIR, "data")
OUT_DIR <- file.path(PROJ_DIR, "outputs")

source(file.path(PROJECT_ROOT, "utils/helpers.R"))
source(file.path(PROJECT_ROOT, "utils/palette_system.R"))
source(file.path(PROJECT_ROOT, "utils/theme_system.R"))

# ── 2. 分组设计 ──────────────────────────────────────────────

CELL_LINES <- c("83", "LM")

# 有序分组：ctrl(基线) -> NFCM(正常成纤维细胞上清) -> CAFCM(癌相关成纤维细胞上清)
#           -> shCAFCM(敲低目的基因后的CAFCM，预期部分救援回ctrl/NFCM)
GROUP_SETS <- list(
  full = c("ctrl", "NFCM", "CAFCM", "shCAFCM"),
  no_NFCM = c("ctrl", "CAFCM", "shCAFCM")
)

# 当前激活的分组变体；后续如需去掉 NFCM 重新分析救援效应，
# 只需改为 "no_NFCM" 并重跑 02→04 脚本，01 的清洗结果无需重新生成
ACTIVE_VARIANT <- "full"

GROUP_COLORS <- setNames(get_colors("NPG", 4), GROUP_SETS$full)

# ── 3. 分析参数（集中管理，便于调整）────────────────────────

P <- list(
  min_count = 10, # 低表达过滤：单个样本中的最低 count
  min_samples = 3, # 至少有 N 个样本通过 min_count
  padj_cutoff = 0.05, # Wald 两两对比的差异显著性阈值
  lfc_cutoff = 1.0, # DEG 计数 / Volcano 用的 log2FC 阈值
  tpm_cutoff = 1, # 聚类前的最低表达量过滤（任一组均值 TPM，可调）
  # LRT 显著性阈值，决定进入聚类的基因全集；83信号较弱，放宽至0.10以获得更多聚类基因，
  # LM 保持标准的0.05（信号已经很强，不需要放宽）
  lrt_padj_cutoff = list("83" = 0.10, "LM" = 0.05),
  min_std = 0.1, # clusterData() 内部低方差基因过滤阈值
  enrich_topn = 12, # 每个簇富集展示的 top term 数
  # elbow 图确定后手动填入对应细胞系/方法的簇数，未填(NA)时聚类脚本会报错，
  # 避免在未经人工确认 K 值的情况下静默产出结果
  k_clusters = list(
    "83" = list(kmeans = 5, mfuzz = 5),
    "LM" = list(kmeans = 5, mfuzz = 5)
  )
)

# ── 4. 输出目录 ──────────────────────────────────────────────

dir.create(
  file.path(OUT_DIR, "00_data"),
  recursive = TRUE,
  showWarnings = FALSE
)
for (cl in CELL_LINES) {
  for (variant in names(GROUP_SETS)) {
    dir.create(
      file.path(OUT_DIR, cl, variant, "DE"),
      recursive = TRUE,
      showWarnings = FALSE
    )
    dir.create(
      file.path(OUT_DIR, cl, variant, "ClusterGVis", "kmeans", "enrichment"),
      recursive = TRUE,
      showWarnings = FALSE
    )
    dir.create(
      file.path(OUT_DIR, cl, variant, "ClusterGVis", "mfuzz", "enrichment"),
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
}

# ── 5. 通用绘图/保存辅助函数（沿用 WY/Semaglutide 项目模式）───

# 保存 ggplot 对象为 PDF + PNG
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

# 保存 ComplexHeatmap / visCluster 图为 PDF
save_ht_pdf <- function(draw_func, filepath, width = 14, height = 12) {
  pdf(filepath, width = width, height = height)
  tryCatch(
    draw_func(),
    error = function(e) message("  [PDF 保存警告] ", conditionMessage(e))
  )
  dev.off()
  message("  已保存: ", filepath)
}

# 将 clusterData() 输出对象（ck/cm）截取为仅含指定簇号的子对象，
# 供 enrichCluster()/visCluster() 使用。
# 注：ClusterGVis 0.99.9 的 visCluster(subgroupAnno=) 与 annoTermData 联用时
# 存在内部 bug ("align_to 与 text 需同名列表")，因此改用手动截取子对象的方式，
# 效果等价（只富集/只展示选中簇），且已验证可行。
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

# 保存 ComplexHeatmap / visCluster 图为 PNG
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

# ClusterGVis::enrichCluster() 内部用未校正的原始 pvalue < pvalueCutoff 筛选/排序
# term（qvalueCutoff/pvalueCutoff 传给 enrichGO/enrichKEGG 时都设成1，即先取全量
# 结果，再手动按原始 pvalue 过滤，p.adjust/qvalue 列在返回前被丢弃，不受多重检验
# 校正）。此函数是它的替代实现：改用 BH 校正后的 p.adjust 过滤/排序/取 topn，
# 输出结构与 enrichCluster() 完全一致（group/Description/pvalue位置放p.adjust值/
# ratio[/geneID]），因此可以直接替换 enrichCluster() 传给 visCluster(annoTermData=)，
# 无需改动可视化代码（visCluster 按列位置而非列名解析，goSize="pval"会自动读到
# 这里放的 p.adjust 值）。
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
