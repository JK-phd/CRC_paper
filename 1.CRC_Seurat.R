library(Seurat)
library(dplyr)
library(AUCell)
library(magrittr)
library(data.table)
library(Matrix)
library(devtools)
library(RcppArmadillo)
library(Rcpp)
library(scales)
library(gplots)
library(ggplot2)
library(cowplot)
library(tibble)
library(data.table)
library(harmony)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(ggrepel)
library(rjags)
library(infercnv)
library(CellChat)
library(ComplexHeatmap)
#library(copykat)
library(scRNAtoolVis)
library(RColorBrewer)
library(MetBrewer)
library(tidyverse)
mycolor <- c(
  "#84B8D7", "#4391A9", "#8DC594", "#72BF5A", "#B49D99", "#F6955A", "#D9CE99", 
  "#FDAE53", "#408DBF", "#D8AC60", "#969D62", "#92CF72", "#659E47", "#F47A79", 
  "#E93A3B", "#B294C7", "#815AA8", "#EBD57C", "#C48244", "#B15928", "#E94330", 
  "#E4986B", "#F06C45", "#FE9E37", "#984ea3", "#377eb8", "#4daf4a", "#e41a1c", 
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#3d7bb0", 
  "#a2c8db", "#5aa554", "#afd195", "#d83f36", "#e99997", "#f2bc7c", "#e68740", 
  "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87", 
  "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", "#8C549C", "#585658", 
  "#9FA3A8", "#E0D4CA", "#a6cee3", "#fdbf6f", "#b2df8a", "#ffff99", "#33a02c", 
  "#1f78b4", "#ff7f00", "#cab2d6", "#6a3d9a", "#fb9a99", "#e31a1c", "#7db0e9", 
  "#3d7bb0", "#57a65b", "#7fcd87", "#b9364e", "#c56585", "#7e4878", "#513c65"
)
colors22 <- c(
  "C0-CD4_T cell" = "#84B8D7", 
  "C1-Follicular B cell" = "#4391A9", 
  "C2-NK" = "#8DC594", 
  "C3-Stem cell" = "#72BF5A", 
  "C4-Plasma cell" = "#B49D99", 
  "C5-CD14_Myeloid cell" = "#F6955A", 
  "C6-NK" = "#D9CE99", 
  "C7-Ciliated cell" = "#FDAE53", 
  "C8-Transit Amplifying Cell" = "#408DBF", 
  "C9-TAM" = "#D8AC60", 
  "C10-Goblet cell" = "#969D62", 
  "C11-Germinal center B cell" = "#92CF72", 
  "C12-Progenitor cell" = "#659E47", 
  "C13-Fibroblast" = "#F47A79", 
  "C14-Transit amplifying cell" = "#E93A3B", 
  "C15-Brush cell" = "#B294C7", 
  "C16-Paneth cell" = "#815AA8", 
  "C17-Progenitor cell" = "#EBD57C", 
  "C18-Enteroendocrine cell" = "#C48244", 
  "C19-Mast cell" = "#B15928", 
  "C20-Endothelial cell" = "#E94330", 
  "C21-Stem cell" = "#E4986B"
)

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
 
if(F){
setwd("path/to/CRC ScRNA/data/") #rawdata
fs=list.files('./','features.tsv.gz')
fs
samples1=gsub('_features.tsv.gz','',fs)
samples1
samples2=samples1

library(stringr)
dir1 <- 'data/Patient0_carcinoma/'
dir2 <- 'data/Patient1_adenoma/'
dir3 <- 'data/Patient1_carcinoma/'
dir4 <- 'data/Patient1_normal/'
dir5 <- 'data/Patient2_adenoma/'
dir6 <- 'data/Patient2_carcinoma/'
dir7 <- 'data/Patient2_normal/'
dir8 <- 'data/Patient2_para-cancer/'
dir9 <- 'data/Patient3_adenoma_1/'
dir10 <- 'data/Patient3_adenoma_2/'
dir11 <- 'data/Patient3_blood/'
dir12 <- 'data/Patient3_carcinoma/'
dir13 <- 'data/Patient3_normal/'
# 
dir.create(dir1, showWarnings = FALSE, recursive = TRUE)
dir.create(dir2, showWarnings = FALSE, recursive = TRUE)
dir.create(dir3, showWarnings = FALSE, recursive = TRUE)
dir.create(dir4, showWarnings = FALSE, recursive = TRUE)
dir.create(dir5, showWarnings = FALSE, recursive = TRUE)
dir.create(dir6, showWarnings = FALSE, recursive = TRUE)
dir.create(dir7, showWarnings = FALSE, recursive = TRUE)
dir.create(dir8, showWarnings = FALSE, recursive = TRUE)
dir.create(dir9, showWarnings = FALSE, recursive = TRUE)
dir.create(dir10, showWarnings = FALSE, recursive = TRUE)
dir.create(dir11, showWarnings = FALSE, recursive = TRUE)
dir.create(dir12, showWarnings = FALSE, recursive = TRUE)
dir.create(dir13, showWarnings = FALSE, recursive = TRUE)

# 
for (i in 1:length(samples2)) {
  patient <- str_split(samples2[i], "_", simplify = TRUE)[2]
  data_type <- str_split(samples2[i], "_", simplify = TRUE)[3]
  
  if (patient == "Patient0" && data_type == "carcinoma") {
    dir_path <- dir1
  } else if (patient == "Patient1" && data_type == "adenoma") {
    dir_path <- dir2
  } else if (patient == "Patient1" && data_type == "carcinoma") {
    dir_path <- dir3
  } else if (patient == "Patient1" && data_type == "normal") {
    dir_path <- dir4
  } else if (patient == "Patient2" && data_type == "adenoma") {
    dir_path <- dir5
  } else if (patient == "Patient2" && data_type == "carcinoma") {
    dir_path <- dir6
  } else if (patient == "Patient2" && data_type == "normal") {
    dir_path <- dir7
  } else if (patient == "Patient2" && data_type == "para-cancer") {
    dir_path <- dir8
  } else if (patient == "Patient3" && data_type == "adenoma_1") {
    dir_path <- dir9
  } else if (patient == "Patient3" && data_type == "adenoma_2") {
    dir_path <- dir10
  } else if (patient == "Patient3" && data_type == "blood") {
    dir_path <- dir11
  } else if (patient == "Patient3" && data_type == "carcinoma") {
    dir_path <- dir12
  } else if (patient == "Patient3" && data_type == "normal") {
    dir_path <- dir13
  }
  
  file.copy(from = paste0(samples1[i], '_features.tsv.gz'),
            to = file.path(dir_path, 'features.tsv.gz'))
  file.copy(from = paste0(samples1[i], '_matrix.mtx.gz'),
            to = file.path(dir_path, 'matrix.mtx.gz'))
  file.copy(from = paste0(samples1[i], '_barcodes.tsv.gz'),
            to = file.path(dir_path, 'barcodes.tsv.gz'))
}

}

setwd('data/')
# 
dir <- list.dirs('./')[-1] 
dir

CRC <- CreateSeuratObject(
  counts = Read10X(data.dir = dir),min.cells = 5,
  min.features = 300 
  
)

print(CRC)
head(CRC)
table(CRC$orig.ident)

sample_map <- c(
  "1" = "P0_carcinoma",
  "2" = "P1_adenoma",
  "3" = "P1_carcinoma",
  "4" = "P1_normal",
  "5" = "P2_adenoma",
  "6" = "P2_carcinoma",
  "7" = "P2_normal",
  "8" = "P2_paracancer",
  "9" = "P3_adenoma1",
  "10" = "P3_adenoma2",
  "11" = "P3_blood",
  "12" = "P3_carcinoma",
  "13" = "P3_normal"
)

CRC$sample <- sample_map[as.character(CRC$orig.ident)]
table(CRC$sample)

dim(CRC)  
CRC[["percent.mt"]] <- PercentageFeatureSet(CRC, pattern = "^MT-.")
head(CRC@meta.data,5)
summary(CRC@meta.data$nCount_RNA)
CRC = subset(CRC,subset =nFeature_RNA >200 &
                nFeature_RNA <8000 & percent.mt < 10)
table(CRC$sample)
CRC = CRC[,!CRC$sample=='P3_blood']

dim(CRC) #[1] 23041 31542
CRC$sample_type <- ifelse(CRC$sample %in% c('P1_normal', 'P2_normal', 'P3_normal'), 'normal',
                          ifelse(CRC$sample %in% c('P1_adenoma', 'P2_adenoma', 'P3_adenoma1', 'P3_adenoma2'), 'adenoma',
                                 ifelse(CRC$sample %in% c('P0_carcinoma', 'P1_carcinoma', 'P2_carcinoma', 'P3_carcinoma'), 'carcinoma',
                                        ifelse(CRC$sample %in% c('P2_paracancer'), 'paracancer', NA))))
table(CRC$sample_type)
CRC<- NormalizeData(CRC)
CRC = FindVariableFeatures(CRC,selection.method = "vst", nfeatures = 2000, verbose = F)
CRC <- ScaleData(CRC)
CRC <- RunPCA(CRC, verbose = FALSE)

CRC <- harmony::RunHarmony(CRC,"sample", plot_convergence = TRUE)

ElbowPlot(CRC, ndims = 30,reduction = "harmony")

CRC <- Seurat::RunUMAP(CRC,reduction = "harmony", dims = 1:20) 
CRC <- Seurat::FindNeighbors(CRC,reduction = "harmony", dims = 1:20) 

CRC <- Seurat::FindClusters(CRC,resolution =0.6,graph.name="RNA_snn")

DimPlot(CRC,group.by = "seurat_clusters",label=T,cols = mycolor)
DimPlot(CRC,group.by = "sample",label=F,cols = mycolor)
DimPlot(CRC,group.by = "sample_type",label=F,cols = colors5)

#sceasy::convertFormat(CRC, from="seurat", to="anndata", outFile='CRC.h5ad')
library(qs)
qsave(CRC,file = 'CRC_ScRNA.qs')

####___________________________________
CRC = qread('CRC_ScRNA.qs')
library(COSG)
library(dplyr)
Indets(CRC)
COSG_markers <- cosg(CRC,groups='all',assay='RNA',slot='data',mu=1,n_genes_user=20)

CRC_markergene = as.data.frame(COSG_markers$names)
view(CRC_markergene)

markergene = head(CRC_markergene,2)
###Fig 1B
DotPlot(object = CRC, features =  c('PTPRC','EPCAM',markergene),scale = T,group.by = "seurat_clusters") + 
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
  
DotPlot(object = CRC, features =  c("^IG-"),scale = T,group.by = "celltype") + 
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

  
CRC_seurat_cluster2 = c('0' = 'CD4_T cell',
                         '1' = 'Follicular B cell',
                         '2' = 'NK',
                         '3'= 'Steam cell',
                         '4' = 'Plasma cell',
                         '5' = 'CD14_Myeloid cell',
                         '6' = 'NK',
                         '7' = 'Ciliated cell',
                         '8' = 'Transit Amplifying Cell',
                         '9' = 'TAM',
                         '10' = 'Goblet cell',
                        '11' = 'Germinal center B cell',
                        '12' = 'Progenitor cell',
                        '13' = 'Fibroblast',
                        '14' = 'Transit amplifying cell',
                        '15' = 'Brush cell',
                        '16' = 'Paneth cell',
                        '17' = 'Progenitor cell',
                        '18' = 'Enteroendocrine cell',
                        '19' = 'Mast cell',
                        '20' = 'Endothelial cell',
                        '21' = 'Stem cell')

CRC[['seurat_clusters2']] = unname(CRC_seurat_cluster2[CRC@meta.data$seurat_clusters])

CRC$new_labers <- paste0('C', CRC$seurat_clusters, '-', CRC$seurat_clusters2)
a <- names(table(CRC$new_labers))
b <- gsub("-.*$", "", a)
c <- a[order(as.numeric(gsub("C", "", b)))]
CRC$celltype <- factor(CRC$new_labers, levels=c)

### Fig 1A
DimPlot(CRC, group.by="celltype", label=T, pt.size = 0.5,label.size = 3,label.box = F,reduction='umap',repel =T,cols = c(mycolor),raster=FALSE)+guides(color=guide_legend(ncol = 1,override.aes = list(size = 3)))
DimPlot(CRC, group.by="seurat_clusters", label=T, split.by = 'sample_type',pt.size = 0.5,label.size = 3,label.box = F,reduction='umap',repel =T,cols = c(mycolor),raster=FALSE)+guides(color=guide_legend(ncol = 1,override.aes = list(size = 2)))+NoLegend()

########
library(SCPA)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(magrittr)
library(qs)
library(msigdbr)
library(ggrepel)
library(ggplot2)
pathways <- c("hallmark", "reactome")
hkr_sets <- msigdbr("Homo sapiens") %>%
  filter(grepl(paste(pathways, collapse = "|"), gs_name, ignore.case = T)) %>%
  format_pathways()
me_pathways <- "combined_metabolic_pathways.csv"
table(CRC$sample_type)

###Fig1 C-F_______________
### G1 vs G2, if FC<0, higher in G2, FC>0, higher in G1;

C = compare_seurat(CRC,group1 = "sample_type", 
                           group1_population = c("adenoma", "normal"),
                           pathways = me_pathways,downsample = 5000,min_genes = 10)  

D = compare_seurat(CRC,group1 = "sample_type", 
                           group1_population = c("adenoma", "normal"), 
                           pathways = hkr_sets,downsample = 5000,min_genes = 10)  

E = compare_seurat(CRC,
                     group1 = "sample_type", 
                     group1_population = c("carcinoma", "adenoma"),
                     pathways = me_pathways,downsample = 3000,min_genes = 10)
F = compare_seurat(CRC,
                     group1 = "sample_type", 
                     group1_population = c("carcinoma", "adenoma"),
                     pathways = hkr_sets,downsample = 3000,min_genes = 10)

ggplot(C[c(1:10),], aes(x = qval, y = reorder(Pathway, qval), color = FC > 0)) +  ### C D E F
  geom_point(size = 6) + 
  geom_segment(aes(x = 0, xend = qval, y = Pathway, yend = Pathway), 
               color = ifelse(C[c(1:10),"FC"] > 0, "#bb78f3", "#c5e9c7"),  ### C D E F
               size = 1) + 
  scale_color_manual(values = c("FALSE" = "#c5e9c7", "TRUE" = "#bb78f3"), 
                     labels = c("normal↑", "adenoma↑")) + 
  labs(x = "Q-value", y = "Pathway") + 
  ggtitle("Metabolic pathway perturbation analysis") + 
  theme_classic() +  
  guides(color=guide_legend(ncol = 1,override.aes = list(size = 5),title = ''))+
  scale_y_discrete(labels = function(x) str_wrap(tolower(gsub("-", " ", gsub("_", " ", x))), width = 30))+
  labs(x = "SCPA qvals", y = "Pathway")+
  theme(axis.text.x = element_text(face = 'italic'),
        axis.text.y = element_text(face = 'italic')) 

#compare_seurat(CRC, 
#                     group1 = "sample_type", 
#                     group1_population = c("carcinoma", "normal"),
#                     pathways = me_pathways,downsample = 3000,min_genes = 10)


###########pathway scores
library(ggplot2)
library(GSEABase)
library(AUCell)
library(ggpubr)

#REACTOME_CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL

####Fig 2A-D

vascular_pathways <- msigdbr("Homo sapiens") %>%
  filter(grepl("REACTOME_CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL", gs_name, ignore.case = T))


geneSets <- vascular_pathways$gene_symbol
matrix <- Seurat::GetAssayData(object=CRC, slot="counts", assay="RNA")
cells_rankings <- AUCell_buildRankings(matrix, nCores=6, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

signature_exp <- data.frame(t(getAUC(cells_AUC)))
CRC <- Seurat::AddMetaData(CRC, as.data.frame(signature_exp))
CRC$sample_type
CRC$geneset1 = CRC$geneSet
median(CRC$geneset1) #
quantile(CRC$geneset1, 0.75) #

a1 = CRC@meta.data %>% 
  group_by(celltype) %>% 
  summarise(m = median(geneset1)) %>% 
  arrange(desc(m)) %>% 
  pull(celltype)

CRC$G1 = factor(CRC$celltype,levels = a1)

p.vascular <- ggviolin(CRC@meta.data, x = "G1", y = "geneset1",
                       color = "G1",
                       add.params = list(color = "G1"),
                       size = 0.2,
                       add = c("boxplot"),
                       palette =colors22) + 
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 60, vjust = 0.5, hjust = 0.5, size = 10, color = c('red','red','red','red',
                                                                                                                 'red','red','red','red',
                                                                                                                 'black','black','black','black',
                                                                                                                 'black','black','black','black',
                                                                                                                 'black','black','black','black',
                                                                                                                 'black','black'), face = "bold")) + 
  #stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  labs(x = '')+ 
  #geom_hline(aes(yintercept=0.05752154), colour="black", linetype="dashed")+
  geom_hline(aes(yintercept=0.04841445), colour="black", linetype="dashed")+ #yintercept=median(CRC$geneset1)
  NoLegend()+
  ggtitle(label = 'Cell surface interactions at the vascular wall - AUCell score')+
  theme(plot.title = element_text(hjust = 0))+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) + 
  theme(panel.grid.major = element_line(size = 0.5, color = "gray90")) + 
  theme(panel.grid.minor = element_line(size = 0.25, color = "gray90"))
 

p.vascular
###_____________________________________________
###Fig 2E
CRC$seurat_clusters <- paste0('C', CRC$RNA_snn_res.0.6)

highlighted_genes <- c("IGKC", "IGHA2",'IGHG4','IGHG2','IGHA1','IGHG1','IGHG3','IGHM','IGLC2','IGLC3')
features = c(grep("^IGK", rownames(CRC), value = TRUE), 
             grep("^IGH", rownames(CRC), value = TRUE), 
             grep("^IGL", rownames(CRC), value = TRUE),
             grep("^IGG", rownames(CRC), value = TRUE),'CD40LG')
gene_colors <- ifelse(features %in% c(highlighted_genes), "red", "black")
gene_sizes <- ifelse(features %in% c(highlighted_genes), 10, 6)

DotPlot(object = CRC, 
        features = c(grep("^IGK", rownames(CRC), value = TRUE), 
                     grep("^IGH", rownames(CRC), value = TRUE), 
                     grep("^IGL", rownames(CRC), value = TRUE),
                     grep("^IGG", rownames(CRC), value = TRUE),'CD40LG'),scale = T, 
        group.by = "seurat_clusters") +scale_colour_gradientn(colors = brewer.pal(9, "YlGnBu")) + 
        theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust = 1,color = c('black','black','black','black','black','black','black',
                                                                              'black','black','black','black',
                                                                              'black','black','black','black',
                                                                              'black','red','black','black','black',
                                                                              'black','black')),
        axis.text.y = element_text(colour = gene_colors,size = c(gene_sizes))) + 
        FontSize(x.text = 8, y.text =6)+coord_flip()


##______________________________________________

## Fig 2H
Pla = CRC[,CRC$celltype=='C4-Plasma cell']
highlighted_genes2 <- c('IGHG3')
features2 =  c("IGKC", "IGHA2",'IGHG4','IGHG2','IGHA1','IGHG1','IGHG3','IGHM','IGLC2','IGLC3')
gene_colors2 <- ifelse(features2 %in% c(highlighted_genes2), "red", "black")
gene_sizes2 <- ifelse(features2 %in% c(highlighted_genes2), 12, 9)

DotPlot(object = Pla, features =  c("IGKC", "IGHA2",'IGHG4','IGHG2','IGHA1','IGHG1','IGHG3','IGHM','IGLC2','IGLC3'),scale = T,
        group.by = "sample_type",dot.scale = 10)+#coord_flip() + ##celltype_l3. ###seurat_clusters
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+theme(axis.text.x = element_text(angle = 45, hjust = 1,color = c('black','black','black','black','black','black','red','black','black','black')),axis.text.y = element_text(colour = gene_colors2,size = c(gene_sizes2)))+
  FontSize(x.text = 10)+ggtitle('C4-Plasma cells - Immunoglobulin-gene expression')+
  theme(plot.title = element_text(hjust = 0.3, size = 9, face = "bold"))
  

dim(Pla)
####
###_________________________________________________________________________
###Neutrophil
table(CRC$celltype)
CRC3 = CRC[,CRC$celltype%in%c('C5-CD14_Myeloid cell','C9-TAM')]

table(CRC3$celltype)
##return to this step:remove HMMR+ cell

CRC3<- NormalizeData(CRC3)
CRC3 = FindVariableFeatures(CRC3,selection.method = "vst", nfeatures = 2000, verbose = F)

CRC3 <- ScaleData(CRC3,verbose = FALSE)
CRC3 <- RunPCA(CRC3,features = VariableFeatures(object = CRC3),verbose = FALSE)
CRC3 <- harmony::RunHarmony(CRC3,"sample", plot_convergence = TRUE)

ElbowPlot(CRC3, ndims = 50,reduction = "harmony")

CRC3 <- Seurat::RunUMAP(CRC3,reduction = "harmony", dims = 1:20) 
CRC3 <- Seurat::FindNeighbors(CRC3,reduction = "harmony", dims = 1:20) 

CRC3 <- Seurat::FindClusters(CRC3,resolution =0.3,graph.name="RNA_snn")

DimPlot(CRC3,group.by = "seurat_clusters",label=T,repel = T,cols = my36colors,pt.size = 0.8)
DimPlot(CRC3,label=F,cols = my36colors,pt.size = 0.8,split.by = 'sample_type')
  DimPlot(CRC3,group.by = "sample",label=F,cols = colors2,pt.size = 0.8)

table(CRC3$seurat_clusters)  
library(COSG)
library(dplyr)
COSG3_markers <- cosg(
  CRC3,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=20)

CRC3_markergene = as.data.frame(COSG3_markers$names)
view(CRC3_markergene)
comma_separated_string3 <- apply(t(CRC3_markergene), 1, function(x) paste(x, collapse = ","))
comma_separated_string3 = as.matrix(comma_separated_string3)


## NANs：SELL、PTSG2、CXCR2、CXCR1、FCGR3B、MME、S100A8、S100A9、S100A12
##TAN：OLR1（LOX-1）、VEGFA、CD83、ICAM1和CXCR4高表达，
#CXCR1、CXCR2、PTGS2、SELL（CD62L）、CSF3R和FCGR3B（CD16B）低表达
markergene3 = head(CRC3_markergene,3)
DotPlot(object = CRC3, features =  c('PTPRC','EPCAM',unique(markergene3)),scale = T,group.by = "seurat_clusters") +
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

DotPlot(object = CRC3, features =  c('VEGFA','CCL20','IL1RN','FCGR3B','CXCR2','G0S2','FCGR2A','PTGS2','IL1B'),scale = T,group.by = "celltype") +
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

DotPlot(object = CRC, features =  c('VEGFA','CCL20','IL1RN','FCGR3B','G0S2','CSF3R','IL1B'),scale = T,group.by = "celltype") + 
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

CRC3_seurat_cluster2 = c('0' = 'classical Monocyte',
                        '1' = 'TAM',
                        '2' = 'Macrophage',
                        '3'= 'VEGFA_TAN',
                        '4' = 'CD1C_DC',
                        '5' = 'LAMP3_DC',
                        '6' = 'non-classical Monocyte',
                        '7' = 'IFN_TAM',
                        '8' = 'Neutrophil(NAN)'
                        )

CRC3[['seurat_clusters2']] = unname(CRC3_seurat_cluster2[CRC3@meta.data$seurat_clusters])


CRC3$new_labers <- paste0('C', CRC3$seurat_clusters, '-', CRC3$seurat_clusters2)
a <- names(table(CRC3$new_labers))
b <- gsub("-.*$", "", a)
c <- a[order(as.numeric(gsub("C", "", b)))]
CRC3$celltype <- factor(CRC3$new_labers, levels=c)


DimPlot(CRC3,group.by = "celltype",label=F,repel = T,label.size = 3,cols = my36colors,pt.size = 0.8)+
  theme(legend.position = 'right',text = element_text(size = 12))


neutrophil_degranulation_pathways <- msigdbr("Homo sapiens") %>%
  filter(grepl("REACTOME_NEUTROPHIL_DEGRANULATION", gs_name, ignore.case = T)) 

if(T){
  geneSets3 <- neutrophil_degranulation_pathways$gene_symbol
   matrix <- Seurat::GetAssayData(object=CRC3, slot="scale.data", assay="RNA")
  cells_rankings <- AUCell_buildRankings(matrix, nCores=6, plotStats=TRUE)
  cells_AUC <- AUCell_calcAUC(geneSets3, cells_rankings)
  
  signature_exp <- data.frame(t(getAUC(cells_AUC)))
  CRC3 <- Seurat::AddMetaData(CRC3, as.data.frame(signature_exp))
  CRC3$geneSet
  CRC3$geneset3 = CRC3$geneSet
  median(CRC3$geneset3)
  quantile(CRC3$geneset3, 0.75)
}

my_comparisons <- list(c('C8-Neutrophil(NAN)','C3-VEGFA_TAN'))
neu3 = CRC3@meta.data %>% 
  group_by(celltype) %>% 
  summarise(m = median(geneset3)) %>% 
  arrange(desc(m)) %>% 
  pull(celltype)

CRC3$neutrophil = factor(CRC3$celltype,levels = neu3)
p.neutrophil <- ggviolin(CRC3@meta.data, x = "neutrophil", y = "geneset3",
                                       color = "neutrophil",
                                       add.params = list(color = "neutrophil"),
                                       size = 0.2,
                                       add = c("boxplot"),
                                       palette =my36colors) + 
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 60,hjust = 0.5,vjust = 0.5,size = 10, face = "bold")) + 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  labs(x = '')+
  geom_hline(aes(yintercept=0.06585859), colour="black", linetype="dashed")+
  NoLegend()+
  ggtitle(label = 'Neutrophil degranulation- All myeloid cells - AUCell score')+
  theme(plot.title = element_text(hjust = 0))+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) + 
  theme(panel.grid.major = element_line(size = 0.5, color = "gray90")) + 
  theme(panel.grid.minor = element_line(size = 0.25, color = "gray90"))

p.neutrophil


####


###_________________________________________________________________________
###Epi
table(CRC$celltype)
table(CRC2$sample_type)

DotPlot(object = CRC, features =  c('NOTUM','MUC2','FCGBP','REG1A','CLCA1','CLCA4','GUCA2A','SOX17'),scale = T,group.by = "celltype") + ##celltype_l3. ###seurat_clusters
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

CRC2 = CRC[,CRC$celltype%in%c('C3-Steam cell','C7-Ciliated cell','C8-Transit Amplifying Cell','C10-Goblet cell','C12-Progenitor cell','C14-Transit amplifying cell','C15-Brush cell',
                              'C16-Paneth cell','C17-Progenitor cell','C21-Steam cell')]
##return this step to clean inummone cell
#CRC2 = CRC2[,!CRC2$seurat_clusters%in%c('3')]
dim(CRC2)
CRC2<- NormalizeData(CRC2)
CRC2 = FindVariableFeatures(CRC2,selection.method = "vst", nfeatures = 2000, verbose = F)

CRC2 <- ScaleData(CRC2,verbose = FALSE)
CRC2 <- RunPCA(CRC2,features = VariableFeatures(object = CRC2),verbose = FALSE)
CRC2 <- harmony::RunHarmony(CRC2,"sample", plot_convergence = TRUE)

ElbowPlot(CRC2, ndims = 50,reduction = "harmony")

CRC2 <- Seurat::RunUMAP(CRC2,reduction = "harmony", dims = 1:4) 
CRC2 <- Seurat::FindNeighbors(CRC2,reduction = "harmony", dims = 1:4) 

CRC2 <- Seurat::FindClusters(CRC2,resolution =0.1,graph.name="RNA_snn")

DimPlot(CRC2,group.by = "seurat_clusters",label=T,cols = mycolor,pt.size = 0.8)+
DimPlot(CRC2,group.by = "sample_type",label=F,cols = colors3,pt.size = 0.8)
DimPlot(CRC2,group.by = "sample",label=F,cols = colors2,pt.size = 0.8)


CRC2_seurat_cluster2 = c('0' = 'NDRG1+cell',
                        '1' = 'CLCA1+cell',
                        '2' = 'BIRC5+cell',
                        '3'= 'CYP2W1+cell',
                        '4' = 'CLCA4+cell',
                        '5' = 'NOTUM+cell')

CRC2[['seurat_clusters2']] = unname(CRC2_seurat_cluster2[CRC2@meta.data$seurat_clusters])
CRC2$new_labers <- paste0('E', CRC2$seurat_clusters, '-', CRC2$seurat_clusters2)
a <- names(table(CRC2$new_labers))
b <- gsub("-.*$", "", a)
c <- a[order(as.numeric(gsub("E", "", b)))]
CRC2$celltype <- factor(CRC2$new_labers, levels=c)


DimPlot(CRC2, group.by="celltype", label=F, pt.size = 0.5,label.size = 3,label.box = F,reduction='umap',repel =T,cols = c(colors2),raster=FALSE)+#theme(legend.text = element_text(size = 10,face = 'bold'))+
DimPlot(CRC2,group.by = "sample_type",label=F,cols = colors3,pt.size = 0.5)

DimPlot(CRC2, group.by="seurat_clusters", label=T, split.by = 'sample_type',pt.size = 0.5,label.size = 3,label.box = F,reduction='umap',repel =T,cols = c(colors2),raster=FALSE)+guides(color=guide_legend(ncol = 1,override.aes = list(size = 2)))+NoLegend()


COSG2_markers <- cosg(
  CRC2,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=20)

CRC2_markergene = as.data.frame(COSG2_markers$names)
view(CRC2_markergene)
comma_separated_string_CRC2 <- apply(t(CRC2_markergene), 1, function(x) paste(x, collapse = ","))
comma_separated_string_CRC2 = as.matrix(comma_separated_string_CRC2)


#markergene2 = head(CRC2_markergene,5)
markergene2
markergene2.2 = list('0'= c('ERO1A','NDRG1','VEGFA'),
                   '1' = c('CLCA1','SPINK4','REG1A'),
                   '2' = c('BIRC5','MKI67','PCLAF'),
                   '3' = c('CYP2W1','GABRA4','HABP2'),
                   '4' = c('CLCA4','CA4','GUCA2A'),
                   '5' = c('NOTUM','GPR155','MAOB'))

DotPlot(object = CRC2, features =  c('PTPRC','EPCAM','LGR5',head(markergene2.2)),scale = T,group.by = "seurat_clusters") + ##celltype_l3. ###seurat_clusters
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))



library(qs)
qsave(CRC2,file = 'CRC2.qs')









