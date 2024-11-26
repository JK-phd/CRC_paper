
library(Seurat)
library(monocle)
library(qs)
source('preR.R')
CRC2 =qread(file = 'CRC2.qs')
table(CRC2$celltype,CRC2$sample_type)


expr_matrix= as(as.matrix(CRC2@assays$RNA@counts), 'sparseMatrix')
dim(expr_matrix)
sample_sheet<-CRC2@meta.data#
sample_sheet$celltype = CRC2$celltype

colnames(sample_sheet)
gene_annotation=data.frame(gene_short_name=rownames(CRC2))
rownames(gene_annotation)=rownames(CRC2)
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd,
                      expressionFamily=negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds) 
colnames(cds)
rownames(cds)

if(T){  
  diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = '~celltype',cores = 6) 
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.00000000000000000000000000000000000001))
  qsave(ordering_genes,file = 'ordering_genes.qs')
}


length(ordering_genes)
if(T){
  RPL_genes <- grep(pattern = c("^RPL"), x = ordering_genes, value = TRUE)
  RPS_genes <- grep(pattern = c("^RPS"), x = ordering_genes, value = TRUE)
  IGH_genes <- grep(pattern = c("^IGH"), x = ordering_genes, value = TRUE)
  IGL_genes <- grep(pattern = c("^IGL"), x = ordering_genes, value = TRUE)
  IGK_genes <- grep(pattern = c("^IGK"), x = ordering_genes, value = TRUE)
  MT_genes <- grep(pattern = c("^MT-"), x = ordering_genes, value = TRUE)
  rm_genes <- c(RPL_genes,RPS_genes,
                MT_genes,
                IGK_genes,IGH_genes,IGL_genes
    )
}
ordering_genes <- ordering_genes[!(ordering_genes %in% rm_genes)]


length(ordering_genes)

cds <- setOrderingFilter(cds, ordering_genes)
#plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2, 
                       reduction_method="DDRTree")


#cds = orderCells(cds,root_state = 6)
cds = orderCells(cds)

plot_cell_trajectory(cds)

dim(cds)
cell_type_color <- unique((pData(cds))$celltype)
cell_type_color
cols = c(colors2)[c(1,3,2,5,6,4)]

names(cols) <- cell_type_color
dim(cds)
#cds$CytoTRACE2_Potency = CRC2$CytoTRACE2_Potency 

print(head(pData(cds)))
m1 =  plot_cell_trajectory(cds, color_by = 'celltype',cell_size = 0.5)+
  facet_wrap('~celltype')+
  scale_color_manual(values = c(cols))+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+theme(legend.position = 'right')+guides(color=guide_legend(ncol = 1,override.aes = list(size = 2)))+NoLegend()
m1

#res = qread('CRC2_cytotrace.qs')
#CRC2$CytoTRACE = res$CytoTRACE
#cds$CytoTRACE = CRC2$CytoTRACE
#m2 =  plot_cell_trajectory(cds, color_by = 'CytoTRACE',cell_size = 1)+
#    scale_colour_viridis_c()+
#  theme(axis.text.x = element_text(size = 10),
#        axis.text.y = element_text(size = 10),
#        axis.title.x = element_text(size = 15),
#        axis.title.y = element_text(size = 15))+theme(legend.position = 'right')+guides(color=guide_legend(ncol =1,override.aes = list(size = 3)))#
#
#m2


m3 =  plot_cell_trajectory(cds, color_by = 'State',cell_size = 1)+
    #facet_wrap('~State')+
  scale_color_manual(values = c(colors6))+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+theme(legend.position = 'right')+guides(color=guide_legend(ncol =1,override.aes = list(size = 3)))

m3

cds$Lineage1

m4 = plot_cell_trajectory(cds, color_by = 'sample_type',cell_size = 1)+
  #  facet_wrap('~copykat.tumor.pred')+
    scale_color_manual(values = c(colors3))+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+theme(legend.position = 'right')+guides(color=guide_legend(ncol = 1,override.aes = list(size = 3)))
m4


m5 = plot_cell_trajectory(cds, color_by = 'Pseudotime',cell_size = 1,show_backbone = T)+
  scale_colour_viridis()+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+theme(legend.position = 'right')

m5

m6 = plot_cell_trajectory(cds, color_by = 'sample',cell_size = 1,show_backbone = T)+
  scale_color_manual(values = c(colors6[5:20]))+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+theme(legend.position = 'right')+theme(legend.position = 'right')+guides(color=guide_legend(ncol =1,override.aes = list(size = 2)))

m6

CRC2$State = cds$State
Idents(CRC2) = CRC2$State

colors7 = c('#e4c2ce','#c0a2ce','#a987c3','#8c6fa3','#884a93','#b96f94','#9e96c8')








