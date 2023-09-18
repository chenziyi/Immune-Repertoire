library(Platypus)
library(magrittr)

bcr.dir <- "/data/zychen/scRNASeq/Cellranger/BCR"
tcr.dir <- "/data/zychen/scRNASeq/Cellranger/TCR"
expres.dir <- "/data/zychen/scRNASeq/Cellranger/5mRNA/华大"

#bcr.dir <- "/home/chenziyi/remote1/zychen/scRNASeq/Cellranger/BCR"
#tcr.dir <- "/home/chenziyi/remote1/zychen/scRNASeq/Cellranger/TCR"
#expres.dir <- "/home/chenziyi/remote1/zychen/scRNASeq/Cellranger/华大"
mixcr <- "/data/zychen/software/anaconda3/bin/mixcr"

bcr.out.directory.list <- list.dirs(bcr.dir, recursive =F)  ### Set directory to the outs folder of cellranger vdj
tcr.out.directory.list <- list.dirs(tcr.dir, recursive =F)  ### Set directory to the outs folder of cellranger vdj
GEX.out.directory.list <- list.dirs(expres.dir, recursive =F) ### Set directory to the /outs folder of cellranger count

samplenames <- c("01", "02")

#Running the VDJ_GEX_matrix function
samplename = "01"
result.dir <- bcr.dir

for(sample in samplenames){
    bcr.out.directory.list.sub <- bcr.out.directory.list[grep(samplename, bcr.out.directory.list)]
    tcr.out.directory.list.sub <- tcr.out.directory.list[grep(samplename, tcr.out.directory.list)]
    GEX.out.directory.list.sub  <- GEX.out.directory.list[grep(samplename, GEX.out.directory.list)
    
    
    vdj_b_cell <- VDJ_analyze(VDJ.out.directory = paste(bcr.out.directory.list.sub, "/outs", sep=""),
                       filter.1HC.1LC=F
                      )
    #_________________________________
    vdj_amino_clonotype <- VDJ_clonotype(VDJ =vdj_b_cell[[1]], clone.strategy="cdr3.aa")
    vdj_germline_clonotype <- VDJ_clonotype(VDJ =vdj_b_cell,
                       clone.strategy="hvj.lvj", global.clonotype = TRUE
                      )
                      
    single_cell <- VDJ_per_clone(clonotype.list=, VDJ.Out.directory= bcr.out.directory)
    #single_cell_tcells <- VDJ_per_clone(clonotype.list=, VDJ.Out.directory= bcr.out.directory)
    
    vdj_region <- call_MIXCR(VDJ.per.clone=single_cell, mixcr.directory="", species="hsa")             
         
    vdj_germline <- VDJ_extract_germline(VDJ.per.clone=single_cell, mixcr.directory="", extract.VDJRegion=T,species="hsa")       
               
    vdj_clonotype_lineage <- VDJ_clonal_lineages(call_MIXCR.output =vdj_region,
                                           VDJ_extract_germline.output=vdj_germline,
                                           as.nucleotiDe=F,
                                           with.germline=T
                                           )                   
                              
    tree <- VDJ_tree(clonal.lineages=, with.germline=T, min.sequences=5, max.sequences=30, unique.sequences=T)                        
                              
    network_clones <- list()
    network_clones[1] <- vdj_b_cell[[1]][1:60, ]
    network_clones[2] <- vdj_b_cell[[3]][1:60, ]
              
    vdj_igraph_10 <- VDJ_network(network_clones[1:2], per.sample=F, distance.cutoff=10, connected=F) 
    vdj_igraph_14 <- VDJ_network(network_clones[1:2], per.sample=F, distance.cutoff=14, connected=F)               
    library(igraph)                             
    igraph::plot.igraph(vdj_igraph_10[[4]], vertex.label=na, vertex.size=7+(0.06 * vdj_igraph_10[[2]]$frequency), vertex.color=factor(vdj_igraph_10[[2]]$human))   
