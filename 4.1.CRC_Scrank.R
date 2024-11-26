

##Scrank
#devtools::install_github("rikenbit/rTensor")
#devtools::install_github("ZJUFanLab/scRank")

library(rTensor)
library(scRank)

#' utile_database
#' 
#' 1. database containing the drug information of inhibitor with drug name, drug target, and interaction type sourced from DGIdb updated on December 30, 2021.
#' 2. database containing the gene information of \code{human} and \code{mouse} updated on June 30, 2021.
#' 3. database containing the transcription factors information of \code{human} and \code{mouse}.
#' @source \url{dgidb.org/drugs} \url{https://www.ncbi.nlm.nih.gov/gene} \url{http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/}
"utile_database" 
library(qs)

CRC2_scrank <- CreateScRank(input = CRC2,
                            #meta = CRC2@meta.data,
                            species = 'human',
                            #cell_type = CRC2$celltype,
                            target = 'CTNNB1'
                            )

CRC2_scrank <- scRank::Constr_net(CRC2_scrank)

library(qs)
library(ggsci)
qsave(CRC2_scrank,file = 'CRC2_scrank.qs')
CRC2_scrank = qread(file = 'CRC2_scrank.qs')

CRC2_scrank <- scRank::rank_celltype(CRC2_scrank)

#CRC2_scrank@cell_type_rank
table(CRC2_scrank@net)
CRC2_scrank = scRank::init_mod(CRC2_scrank)

plot_dim(CRC2_scrank)

plot_dim(CRC2_scrank) + 
  scale_color_viridis_c()+ 
  theme_bw() + 
  theme(legend.position = "top", 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10), 
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 10, face = "bold")) + 
  guides(color = guide_legend(title.position = "top", 
                              title.hjust = 0.5, 
                              override.aes = list(size = 5))) + 
  labs(x = "Dimension 1", 
       y = "Dimension 2", 
       color = "scRank - Ranking of drug sensitivity targeted at CTNNB1 from low to high") + 
  theme(panel.grid.major = element_line(size = 0.5), 
        panel.grid.minor = element_line(size = 0.2))





