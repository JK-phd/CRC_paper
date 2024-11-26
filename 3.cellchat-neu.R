#____________________________________________________________________________________________________________________________________________________________________#
  # CRC   CellChat 
  #devtools::install_github("jinworks/CellChat")
  #devtools::install_github('immunogenomics/presto')
  library(CellChat)
  library(ggplot2)
  library(ggalluvial)
  library(svglite)
  library(Seurat)
  library(circlize)
  library(mindr)
  
  cellchat <- createCellChat(CRC@assays$RNA@data, meta = CRC@meta.data, group.by = "new_celltype")
  levels(cellchat@idents)
  groupSize <- as.numeric(table(cellchat@idents)) 
  
  
  CellChatDB <- CellChatDB.human 
  # showDatabaseCategory(CellChatDB)
  
  dplyr::glimpse(CellChatDB$interaction)
  CellChatDB.use = CellChatDB
  # CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling",'ECM-Receptor','Cell-Cell Contact')) # use Secreted Signaling for cell-cell communication analysis
  #CellChatDB.use2 <- subsetDB(CellChatDB, search = c('KEGG')) # use Secreted Signaling for cell-cell communication analysis
  
  cellchat@DB <- CellChatDB.use # set the used database in the object
  
  unique(CellChatDB$interaction$interaction_name)
 
   gc()
  cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  future::plan("multisession", workers = 4)  # do parallel  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)  
  
  
  cellchat <- computeCommunProb(cellchat,type = "triMean") 
  cellchat <- filterCommunication(cellchat, min.cells = 50)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  library(qs)
  qsave(cellchat,file = 'CRC_neu+allcelltype_cellchat.qs')
  
  levels(cellchat@idents)
  
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                   weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                   weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  mat <- cellchat@net$weight
  par(mfrow = c(4,6), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  
  
  cellchat@netP$pathways
  pathways.show <- c("MIF")  
  
  levels(cellchat@idents)   
  
  vertex.receiver = c(4,1,9) 
  
  netVisual_aggregate(cellchat, signaling = pathways.show,  
                      vertex.receiver = vertex.receiver,layout = "hierarchy")
  
  netVisual_bubble(cellchat, sources.use = c(1:31), 
                   targets.use = c(4,9), remove.isolate = FALSE)+
    theme(axis.text.x.bottom = element_text(angle = 90,size = 8,face = 'bold',color = c('black','black','black','black','black',
                                                                                                               'black','black','black','black','black',
                                                                                                               'black','black','black','black','black',
                                                                                                               'black','black','black','black','black',
                                                                                                               'black','black','black','black','black',
                                                                                                               'black','black','black','black','black')),
          axis.text.y.left = element_text(size = 8,
                                          face = "bold",color = c('black','black','black','black','black',
                                                                  'black','black','black','black','black',
                                                                  'black','black','black','black','black',
                                                                  'black','black','black','black','black',
                                                                  'black','black','black','black','black',
                                                                  'black','black','black','black','black',
                                                                  'black','black','black','black','black',
                                                                  'black','black','black','black','black',
                                                                  'black','black','black','black','red',   #red
                                                                  'black','black','black','black','black',
                                                                  'black','black','black','black','black',
                                                                  'black','black','black','black','black',
                                                                  'black','black','black','black','black',
                                                                  'black','black','black','black','black')))+
    coord_flip()+
    labs(x = '')
  
  
  
  
  
  
  
  
  
  
  
  
  
