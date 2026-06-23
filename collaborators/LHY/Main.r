## from CRAN
CRAN.packages <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies = TRUE)
  }
}
CRAN.packages(c(
  "ROCR",
  "quantmod",
  "xts",
  "RcppArmadillo",
  "Rcpp",
  "MASS",
  "parallel",
  "devtools",
  "methods",
  "igraph",
  "gridExtra",
  "grid",
  "ggplot2",
  "gplots",
  "reticulate"
))

## from Bioconductor
Bioconductor.packages <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    BiocManager::install(new.pkg, dependencies = TRUE)
  }
}
Bioconductor.packages(c(
  "DESeq2",
  "limma",
  "ConsensusClusterPlus",
  "pathview",
  "KEGGgraph",
  "KEGGREST",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "org.Rn.eg.db",
  "org.Dm.eg.db",
  "reactome.db"
))

library(devtools)
install_github("https://github.com/CAMO-R/Rpackage")


library(CAMO)

data("hb") #load human burn data
data("hb.group")


summaryDE = indDE(
  data = hb,
  group = as.factor(hb.group),
  data.type = "microarray",
  case.label = "2",
  ctrl.label = "1"
) #differential anlysis
hb_pData = summaryDE[, c(3, 1)]
hb_MCMCout = bayes(hb_pData, seed = 12345) #Bayesian differential analysis

data("hs") #load human sepsis data
data("hs.group")
summaryDE = indDE(
  data = hs,
  group = as.factor(hs.group),
  data.type = "microarray",
  case.label = "2",
  ctrl.label = "1"
) #differential anlysis
hs_pData = summaryDE[, c(3, 1)]
hs_MCMCout = bayes(hs_pData, seed = 12345) #Bayesian differential analysis

data("mb") #load mouse burn data
data("mb.group")
summaryDE = indDE(
  data = mb,
  group = as.factor(mb.group),
  data.type = "microarray",
  case.label = "2",
  ctrl.label = "1"
) #differential anlysis
mb_pData = summaryDE[, c(3, 1)]
mb_MCMCout = bayes(mb_pData, seed = 12345) #Bayesian differential analysis

data("ms") #load mouse sepsis data
data("ms.group")
summaryDE = indDE(
  data = ms,
  group = as.factor(ms.group),
  data.type = "microarray",
  case.label = "2",
  ctrl.label = "1"
) #differential anlysis
ms_pData = summaryDE[, c(3, 1)]
ms_MCMCout = bayes(ms_pData, seed = 12345) #Bayesian differential analysis


mcmc.list = list(hb_MCMCout, hs_MCMCout, mb_MCMCout, ms_MCMCout)
species = c(rep("human", 2), rep("mouse", 2)) #specify species for each MCMC matrix
data(hm_orth) ##load orthologs file, retrieved from:https://fgr.hms.harvard.edu/diopt
mcmc.merge.list <- merge(
  mcmc.list,
  species = species,
  ortholog.db = hm_orth,
  reference = 1
)
save(mcmc.merge.list, file = "mcmc.merge.list.RData")


dataset.names = c("hb", "hs", "mb", "ms") #specify study names for each merged MCMC matrix
set.seed(12345)
ACS_ADS_global = multi_ACS_ADS_global(
  mcmc.merge.list,
  dataset.names,
  measure = "Fmeasure",
  B = 1000
)
save(ACS_ADS_global, file = "ACS_ADS_global.RData")


#Select pathways by meta-enrichement-analysis
data(human.pathway.list) #load pathway database, this includes KEGG and Reactome database
select.pathway = pathSelect(
  mcmc.merge.list,
  pathway.list,
  pathwaysize.lower.cut = 5,
  pathwaysize.upper.cut = 200,
  overlapsize.cut = 5,
  med.de.cut = 3,
  min.de.cut = 0,
  qfisher.cut = 0.05
)
(K = length(select.pathway))
select.pathway.list = pathway.list[select.pathway]
# save(select.pathway.list,file=paste(WD,"/select_",K,"_pathways.RData",sep=""))

set.seed(12345)
ACS_ADS_pathway = multi_ACS_ADS_pathway(
  mcmc.merge.list,
  dataset.names,
  select.pathway.list,
  measure = "Fmeasure",
  B = 1000,
  parallel = F
)
save(ACS_ADS_pathway, file = "ACS_ADS_pathway.RData")
