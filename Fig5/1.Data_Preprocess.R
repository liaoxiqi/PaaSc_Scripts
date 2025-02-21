library(Seurat)
library(PaaSc)
library(reticulate)
use_condaenv(condaenv = "xxx/miniconda3/envs/R4.1/bin/python", required = TRUE)
ad <- import("anndata")


#-----H5ad to Seurat-----#
adata <- ad$read_h5ad("xxx/PaaSc_Data/TMS/tabula-muris-senis-facs-official-raw-obj.h5ad") 
colnames(adata$X) <- adata$var$features

celltype = c('proerythroblast','granulocyte monocyte progenitor cell','classical monocyte','dendritic cell',"CD4-positive, alpha-beta T cell","CD8-positive, alpha-beta T cell",'regulatory T cell','medium spiny neuron','interneuron','neuron','oligodendrocyte','oligodendrocyte precursor cell','Bergmann glial cell','astrocyte','hepatocyte','chondrocyte','bladder cell','ventricular myocyte','pancreatic B cell')

cells = rownames(adata$obs[adata$obs$cell_ontology_class %in% celltype, ])

obj <- Seurat::CreateSeuratObject(counts = t(as.matrix(adata$X))[,cells], project = "TMS facs",
                                         meta.data =  adata$obs[cells,],
                                         min.cells = 0, min.features = 0)
obj = computeMCA(obj,nmcs = 20)
saveRDS(obj,file = 'xxx/PaaSc_Data/TMS/TMS_facs_selected.RDS')

#-----Traits-----#

gs_files = list.files('xxx/gs_file/magma_10kb_top1000_zscore.75_traits.batch/',full.names = T)
gene_list_all <- lapply(gs_files, function(file) {
                          genesets <- read.table(file, header = TRUE)
                          gene_list <- setNames(
                                                lapply(genesets$GENESET, function(x) {
                                                  genes <- strsplit(x, ",")[[1]] 
                                                  genes <- sapply(genes, function(pair) strsplit(pair, ":")[[1]][1])
                                                  genes <- as.vector(genes)
                                                  return(genes)
                                                }),
                                                genesets$TRAIT
                          )
                      return(gene_list)
                    })

# Human
human_traits <- do.call(c, gene_list_all)

# Mouse
gene_map <- read.table('mouse_human_homologs.txt',header = T)                  
convert_to_mouse <- function(geneset, gene_map) {
        mouse_genes <- gene_map$MOUSE_GENE_SYM[match(geneset, gene_map$HUMAN_GENE_SYM)]
        mouse_genes <- mouse_genes[!is.na(mouse_genes)]
  return(mouse_genes)
}
mouse_traits <- lapply(human_traits, convert_to_mouse, gene_map = gene_map)
saveRDS(mouse_traits,file='mouse_traits_geneset.RDS')  

