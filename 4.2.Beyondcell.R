
#install_github("ltorgo/DMwR2",ref="master")
#devtools::install_github('cran/DMwR')
#devtools::install_local("beyondcell-master.zip")

library(Seurat)
library(beyondcell)
library(qs)
sc = CRC2
DefaultAssay(sc) = "RNA"


#（PSc）：Captures the transcriptional changes induced by a drug
#（SSc）： Captures the sensitivity to a given drug.
# Generate geneset object with one of the ready to use signature collections.
gs <- GetCollection(PSc)
# You can deactivate the functional pathways option if you are not interested in evaluating them.
nopath <- GetCollection(PSc, include.pathways = FALSE)

# Compute BCS for the PSc. This might take a few minutes depending on the size of your dataset.
bc <- bcScore(sc, gs, expr.thres = 0.1) 

qsave(bc,file = 'bc.qs')
bc =qread('bc.qs')

# Run the UMAP reduction. 
bc <- bcUMAP(bc, k.neighbors = 4, res = 0.2)
# Run the bcUMAP function again, specifying the number of principal components you want to use.
bc <- bcUMAP(bc, pc = 10, k.neighbors = 4, res = 0.1)

# Visualize whether the cells are clustered based on the number of genes detected per each cell.
bcClusters(bc, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE, pt.size = 0.8)

bcClusters(bc, UMAP = "beyondcell", idents = "celltype", factor.col = TRUE, pt.size = 0.8)

bc <- bcRegressOut(bc, vars.to.regress = "nFeature_RNA")


# Recompute the UMAP.
bc <- bcUMAP(bc, pc = 10, k.neighbors = 20, res = 0.2)
# Visualize the UMAP.
bcClusters(bc, UMAP = "beyondcell", idents = "bc_clusters_res.0.2", pt.size = 0.6)+ggtitle(label = 'Beyondcell classify all epithelial cells into different TCs')+theme_bw()
# Visualize the therapeutic clusters.
bcClusters(bc, UMAP = "Seurat", idents = "bc_clusters_res.0.2", pt.size = 0.6)+ggtitle(label = 'Mapping TCs to Seurat dimensionality reduction results')+theme_bw()
bcClusters(bc, UMAP = "Seurat", idents = "celltype", pt.size = 0.6)+ggtitle(label = 'Celltype')+theme_bw()+
  scale_color_manual(values = c("E0-NDRG1+cell" =colors2[1], "E1-CLCA1+cell" = colors2[2],'E2-BIRC5+cell' = colors2[3],'E3-CYP2W1+cell' = colors2[4],'E4-CLCA4+cell' = colors2[5],'E5-NOTUM+cell' =colors2[6]))

# Obtain condition-based statistics.
#bc@meta.data$sample_type
bc <- bcRanks(bc, idents = "celltype")

# Obtain unextended therapeutic cluster-based statistics.
bc <- bcRanks(bc, idents = "bc_clusters_res.0.2")


# Explore the statistics table.
head(bc@ranks$celltype)
bc_rank =  as.data.frame(bc@ranks[["bc_clusters_res.0.2"]])
View(bc_rank)
write.csv(bc_rank,file = 'beyondcell result.csv')
#可以使用该函数获取汇总表。此表包括 SP、均值、中位数、标准差、方差、最小值、最大值、NaN 比例等指标 
#以及每个签名的残差平均值。此外，还计算了签名等级 考虑 SP 和残差均值。
#此表旨在帮助您 特定细胞/斑点组的候选药物的优先级。
#The ranking returned by orders the drug signatures from most to least sensitive for a specific group of cells/spots.
#This kind of rank might be useful to inspect intratumoural heterogeneity (ITH). 

#The scaled BCS ranges from 0 to 1 and measures the cell perturbation susceptibility 
#(when using the PSC) or the predicted sensitivity to a given drug (when using the SSC).
#缩放后的 BCS 范围从 0 到 1，用于衡量细胞扰动的易感性（使用 PSC 时）或对特定药物的预测敏感性（使用 SSC 时）。

#The BCS can also be used to evaluate the cells’ functional status if functional signatures are applied.
#Furthermore, a switch point (SP) is calculated for each analysed signature, by determining the value in the 0 to 1 
#scale where cells switch from a down-regulated status to an up-regulated one. 
#此外，为每个分析的签名计算一个转换点（SP），通过确定在 0 到 1 的量表上细胞从下调状态转换为上调状态的值

#Thus, the most therapeutically-homogeneous tumours would be those in which each and every one of their cells responds with the 
#same directionality to a certain drug, either towards sensitivity (SP = 0) or resistance (SP = 1), 
#因此，最具治疗同质性的肿瘤将是那些其每个细胞对某一药物都以相同方向响应的肿瘤，无论是朝向敏感性（SP = 0）还是抗性（SP = 1）。
#while a heterogeneous response would be represented by intermediate SPs
#而异质性响应将由中间的 SP 值表示。

#敏感性（SP = 0）或耐药性（SP = 1）反应

bcClusters(bc, UMAP = "Seurat", idents = "celltype", pt.size = 0.8)

FindDrugs(bc, "PNU-74654")
FindDrugs(bc, "TALNIFLUMATE")
FindDrugs(bc, "Gefitinib")
FindDrugs(bc, "Erlotinib")
FindDrugs(bc, "Lapatinib")
FindDrugs(bc, "XAV-939")
FindDrugs(bc, "IWP-2")
FindDrugs(bc, "RUXOLITINIB")
FindDrugs(bc, "TOFACITINIB")
FindDrugs(bc, "Vemurafenib")
FindDrugs(bc, "Dabrafenib")
FindDrugs(bc, "Trametinib")
FindDrugs(bc, "Selumetinib") 
FindDrugs(bc, "Bortezomib")
FindDrugs(bc, "BX-795")



bcSignatures(bc, UMAP = "Seurat", signatures = list(values = "sig-6957"), pt.size = 0.8)+ theme_bw() ##change the values FindDrugs()

bcSignatures(bc, UMAP = "Seurat", signatures = list(values = "sig-6957"), pt.size = 0.8)+ 
 # scale_color_viridis_c()+ 
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 10), 
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 12, face = "bold")) + 
  guides(color = guide_legend(title.position = "top", 
                              title.hjust = 0.5, 
                              override.aes = list(size = 5))) + 
  labs(x = "Dimension 1", 
       y = "Dimension 2", 
       color = "Beyondcell - Rank of drug sensitivity from low to high") + 
  theme(panel.grid.major = element_line(size = 0.5), 
        panel.grid.minor = element_line(size = 0.2))


bc4Squares(bc,idents = "bc_clusters_res.0.2", lvl = '0',pt.size = 0.1)


# General view.
bcHistogram(bc, signatures = "sig-6957", idents = NULL)
# Condition-based histograms.
bcHistogram(bc, signatures = "sig-6957", idents = "celltype")
bcHistogram(bc, signatures = "sig-7134", idents = "celltype")
bcHistogram(bc, signatures = "sig-6611", idents = "celltype")
bcHistogram(bc, signatures = "sig-20603", idents = "celltype")
bcHistogram(bc, signatures = "sig-6957", idents = "celltype")



bc_rank0 = bc_rank[,1:16]
tyrosine_kinase_receptor_rows <- which(grepl(c("tyrosine kinase receptor"), bc_rank0$MoAs, ignore.case = TRUE))
bc_rank_0 <- bc_rank0[tyrosine_kinase_receptor_rows, ]
bc_rank_0 = bc_rank_0[order(bc_rank_0$switch.point.0),]
write.csv(bc_rank_0,file = 'bc_rank_0.csv')


BETA_CATENIN_rows <- which(grepl(c("BETA-CATENIN"), bc_rank0$MoAs, ignore.case = TRUE))
bc_rank_BETA_CATENIN_rows <- bc_rank[BETA_CATENIN_rows, ]
bc_rank_BETA_CATENIN_rows = bc_rank_BETA_CATENIN_rows[order(bc_rank_BETA_CATENIN_rows$switch.point.0),]
write.csv(bc_rank_BETA_CATENIN_rows,file = 'bc_rank_BETA_CATENIN_rows.csv')


INTERLEUKIN_rows <- which(grepl(c("INTERLEUKIN"), bc_rank0$MoAs, ignore.case = TRUE))
bc_rank_INTERLEUKIN_rows <- bc_rank[INTERLEUKIN_rows, ]
bc_rank_INTERLEUKIN_rows = bc_rank_INTERLEUKIN_rows[order(bc_rank_INTERLEUKIN_rows$switch.point.0),]
write.csv(bc_rank_INTERLEUKIN_rows,file = 'bc_rank_INTERLEUKIN_rows.csv')

MAPK_INHIBITOR_rows <- which(grepl(c("MAPK INHIBITOR"), bc_rank0$MoAs, ignore.case = TRUE))
bc_rank_MAPK_INHIBITOR_rows <- bc_rank[MAPK_INHIBITOR_rows, ]
bc_rank_MAPK_INHIBITOR_rows = bc_rank_MAPK_INHIBITOR_rows[order(bc_rank_MAPK_INHIBITOR_rows$switch.point.0),]
write.csv(bc_rank_MAPK_INHIBITOR_rows,file = 'bc_rank_MAPK_INHIBITOR_rows.csv')

MAPK_INHIBITOR_rows <- which(grepl(c("MAPK INHIBITOR"), bc_rank0$MoAs, ignore.case = TRUE))
bc_rank_MAPK_INHIBITOR_rows <- bc_rank[MAPK_INHIBITOR_rows, ]
bc_rank_MAPK_INHIBITOR_rows = bc_rank_MAPK_INHIBITOR_rows[order(bc_rank_MAPK_INHIBITOR_rows$switch.point.0),]
write.csv(bc_rank_MAPK_INHIBITOR_rows,file = 'bc_rank_MAPK_INHIBITOR_rows.csv')




bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-7039')), pt.size = 0.8)+ theme_bw()
bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-7222')), pt.size = 0.8)+ theme_bw()
bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-24')), pt.size = 0.8)+ theme_bw()
bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-1386')), pt.size = 0.8)+ theme_bw()
bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-19256')), pt.size = 0.8)+ theme_bw()
bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-17067')), pt.size = 0.8)+ theme_bw()
bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-6957')), pt.size = 0.8)+ theme_bw()
bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-2431')), pt.size = 0.8)+ theme_bw()
bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-570')), pt.size = 0.8)+ theme_bw()
bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-16365')), pt.size = 0.8)+ theme_bw()
bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-11767')), pt.size = 0.8)+ theme_bw()
bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-1376')), pt.size = 0.8)+ theme_bw()
bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-3552')), pt.size = 0.8)+ theme_bw()
bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-16073')), pt.size = 0.8)+ theme_bw()
bcSignatures(bc, UMAP = "Seurat", signatures = list(values = c('sig-3379')), pt.size = 0.8)+ theme_bw()














