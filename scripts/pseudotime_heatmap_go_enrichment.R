library(ClusterGVis)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(monocle)
library(WGCNA)
library(hdWGCNA)
library(future)
library(tidyverse,ggplot2)
library(patchwork)
library(reshape2)
library(ggdark)
library(cluster)
library(png)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(cowplot)
library(patchwork)
library(scDblFinder)
library(harmony)
library(viridis)
library(RColorBrewer)
library(scCustomize)
library(SoupX)
library(DropletUtils)
#library(DoubletFinder)
# enable parallel processing for network analysis (optional)
enableWGCNAThreads(nThreads = 8)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

setwd("C:/Users/zhangdeng/Desktop/single_cell_MX")
HSMM <- readRDS("HSMM_myo_final.rds")

library('org.MvBL.eg.db')
library(clusterProfiler)
library(ClusterGVis)

raw <- read.csv(file.path("umap基因.csv"), header = FALSE,
                stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
raw <- raw[!is.na(raw[, 1]) & raw[, 1] != "", , drop = FALSE]
raw[, 1] <- trimws(raw[, 1])
raw[, 2] <- trimws(raw[, 2])
colnames(raw) <- c("gene", "Symbol")
cat("Markers in CSV:", nrow(raw), "\n")
allgenes <- rownames(fData(HSMM))
raw <- raw[raw$gene %in% allgenes, ]
cat("Markers in CDS:", nrow(raw), "\n")

gene_list <- raw$gene
names(gene_list) <- raw$Symbol

mark_labels <- raw$Symbol

diff <- read.csv("diff_test_results.csv",
                 stringsAsFactors=FALSE)
sig <- diff$gene_short_name[!is.na(diff$qval) & diff$qval < 0.01]
sig <- intersect(sig, rownames(HSMM))
df <- plot_pseudotime_heatmap2(HSMM[sig,], num_clusters = 3,cores = 1)

genus    <- "Erythranthe"   # 属名
species  <- "verbenacea"    # 种名
pkg_name <- paste0("org.", toupper(substr(genus, 1, 1)), species, ".eg.db")
install.packages(file.path('orgdb', pkg_name), repos = NULL, type = "source")
library(pkg_name, character.only = TRUE)

genes_cluster1 <- df$wide.res %>%
  filter(cluster == 1) %>%
  pull(gene)
genes_cluster2 <- df$wide.res %>%
  filter(cluster == 2) %>%
  pull(gene)
genes_cluster3 <- df$wide.res %>%
  filter(cluster == 3) %>%
  pull(gene)

ego1 <- enrichGO(
  gene         =sub("\\(.*\\)", "", genes_cluster1),       # 一组 AT 编号
  OrgDb        = org.Everbenacea.eg.db,
  keyType      = "GID",
  ont          = "BP",        # BP: Biological Process
  pAdjustMethod= "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

ego2 <- enrichGO(
  gene         = sub("\\(.*\\)", "", genes_cluster2),       # 一组 AT 编号
  OrgDb        = org.Everbenacea.eg.db,
  keyType      = "GID",
  ont          = "BP",        # BP: Biological Process
  pAdjustMethod= "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

ego3 <- enrichGO(
  gene         = sub("\\(.*\\)", "", genes_cluster3),       # 一组 AT 编号
  OrgDb        = org.Everbenacea.eg.db,
  keyType      = "GID",
  ont          = "BP",        # BP: Biological Process
  pAdjustMethod= "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

ego1 <- gofilter(ego1, level = 5)
ego2 <- gofilter(ego2, level = 5)
ego3 <- gofilter(ego3, level = 5)
# 假设你的 enrichGO 结果保存在 ego 对象中


ego1_df <- as.data.frame(ego1)[1:5,]
ego1_df$group <- "C1"

ego2_df <- as.data.frame(ego2)[c(10,11,12,13,15),]
ego2_df$group <- "C2"

ego3_df <- as.data.frame(ego3)[1:5,]
ego3_df$group <- "C3"

all_ego_df <- rbind(ego1_df, ego2_df, ego3_df)

# 提取并重命吊你要的列
formatted_ego <- all_ego_df[, c( "group", "Description", "pvalue", "Count")]

# 重命吊列（ID → GO term, Count → ratio）
colnames(formatted_ego) <- c( "group", "Description", "pvalue", "ratio")




mp <- setNames(raw$Symbol, raw$gene)
df$wide.res$gene <- ifelse(df$wide.res$gene %in% raw$gene,
                           mp[df$wide.res$gene], df$wide.res$gene)
df$long.res$gene <- ifelse(df$long.res$gene %in% raw$gene,
                           mp[df$long.res$gene], df$long.res$gene)

mark_labels <- raw$Symbol


pdf(file = "gh2.pdf",height = 8,width = 16,onefile = F)
visCluster(object = df,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = mark_labels,
           markGenes.side = "left",
           annoTerm.data = formatted_ego,
           go.col = rep(jjAnno::useMyCol("calm",n = 3),each = 5),
           add.bar = T,
           line.side = "left")
dev.off()
