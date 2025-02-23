suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(pbapply))

#-----FindMarkers-----#
allen_brain =readRDS('xxx/PaaSc_Data/MouseBrain/allen_brain.rds')
Idents(allen_brain) <- "subclass"
CellType=unique(allen_brain@meta.data$subclass)

allen_brain_markers  = pblapply(CellType,function(x){
                  print(x)  
                  markers =  FindMarkers(allen_brain,ident.1 = x)
                  markers =  rownames(markers[markers$p_val_adj < 0.05, ])
                  
    return(markers[1:200])
})

names(allen_brain_markers ) = CellType
#saveRDS(allen_brain_markers ,file='xxx/PaaSc_Data/MouseBrain/allen_brain_Top200_cellmarker.RDS')


#-----Data load & process-----#
suppressPackageStartupMessages(library(anndata))
suppressPackageStartupMessages(library(reticulate))
use_condaenv(condaenv = "/sibcb1/bioinformatics/liaoxiqi/miniconda3/envs/R4.1/bin/python", required = TRUE)
ad <- import("anndata")

adata <- ad$read_h5ad("/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/MouseBrain/adata_benchmark.h5ad")
colnames(adata$X) <- adata$var$features
obj <- Seurat::CreateSeuratObject(counts =  t(as.matrix(adata$X)), project = "Visium_Mousebrain",
                                         meta.data =  adata$obs,
                                         min.cells = 0, min.features = 0)
obj = computeMCA(obj,nmcs = 20)
#saveRDS(obj,file = '/sibcb1/bioinformatics/liaoxiqi/project/PaaSc/MouseBrain/Visium_Mousebrain.RDS')

#-----Compute Score-----#

#obj= readRDS('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/MouseBrain/Visium_Mousebrain.RDS')
#allen_brain_markers = readRDS('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/MouseBrain/allen_brain_Top200_cellmarker.RDS')
keggs = readRDS('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/MouseBrain/mouse_kegg_geneset.RDS')
gene_rate = getGeneRate(background.geneset = keggs, pathway.geneset = allen_brain_markers, mode ='single')
regression_data = doRegression(obj, gene.rate = gene_rate)
score_data = computeScore(obj, regression.data = regression_data, pvalue = 0.05, weight = FALSE, normalize = "z-score")
write.table(score_data,'/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/MouseBrain/Visium_Mousebrain_Structure_related_Top200markers_PaaSc_score.tsv',sep='\t',col.names=T,row.names=T,quote=F)

#-----Compute Spatial MI-----#
spatial=as.data.frame(adata$obsm)
x <- spatial$spatial.1
y <- spatial$spatial.2
spatialMI <- calculateSpatialMI(x = x, y = y, z = score_data)
spatialMI <- data.frame(spatialMI, qvalue = p.adjust(spatialMI$pvalue))
spatialMI[order(spatialMI$MI,decreasing = T),]
write.table(spatialMI[order(spatialMI$MI,decreasing = T),],'/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/MouseBrain/spatialMI.txt',sep = '\t',col.names = T,row.names = T,quote=F)
