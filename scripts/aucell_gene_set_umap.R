library(AUCell)

gene_set <- c(
  'MV7G041300',
  'MV1G093900'
)

# 提取表达矩阵
expr_mat <- GetAssayData(combined, layer = "data")
gene_set <- intersect(gene_set, rownames(expr_mat))
gene_sets <- list("CRC_YAB_gene_set" = gene_set)
rankings <- AUCell_buildRankings(expr_mat, plotStats = FALSE)
aucell_result <- AUCell_calcAUC(gene_sets, rankings)
auc_scores <- getAUC(aucell_result)
combined$AUCell_score <- as.numeric(auc_scores[1,])
assignments <- AUCell_exploreThresholds(aucell_result, plotHist = FALSE)

FeaturePlot(combined, features = "AUCell_score",order=T, reduction = "umap") +
  scale_colour_viridis_c(name = "AUCell Score") +
  ggtitle("UMAP colored by AUCell Score\n(ARF3,ARF4)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))+
  geom_polygon(
    data = outline,
    aes(
      UMAP_1,
      UMAP_2
    ),
    fill = NA,
    colour = "red",
    linewidth = 1.2
  )
