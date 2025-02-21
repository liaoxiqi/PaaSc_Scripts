library(Seurat)
library(PaaSc)
library(CelliD)
library(GSDensity)


mouse_gene_sets = readRDS('xxx/PaaSc_Data/traits/mouse_traits_geneset.RDS')
keggs = readRDS('xxx/PaaSc_Data/traits/mouse_kegg_geneset.RDS')
obj = readRDS('xxx/PaaSc_Data/TMS/TMS_facs_selected.RDS')

###PaaSc
gene_rate = getGeneRate(background.geneset = keggs, pathway.geneset = mouse_gene_sets,mode ='single')
regression_data = doRegression(obj, gene.rate = gene_rate)
score_data = computeScore(obj, regression.data = regression_data, pvalue = 0.05, weight = FALSE, normalize = "z-score")

write.table(score_data,'xxx/PaaSc_Data/TMS/TMS_facs_selected_alltraits_PaaSc_score.tsv',sep='\t',col.names=T,row.names=T,quote=F)


###GSDensity
ce <- compute.mca(object = obj)
cells <- colnames(obj)
el <- compute.nn.edges(coembed = ce, nn.use = 300)
GSD <- run.rwr(el = el, gene_set = mouse_gene_sets, cells = cells)
GSD.res <- as.data.frame(GSD)
write.table(GSDensity.res,'xxx/PaaSc_Data/TMS/TMS_facs_selected_alltraits_GSDensity_score.tsv',sep='\t',col.names=T,row.names=T,quote=F)


###CelliD
set.seed(123)
HGT <-RunCellHGT(X = obj,
            reduction = "mca",
            pathways = mouse_gene_sets,
            log.trans = T, minSize = 5,
            dims=1:20
        )
CelliD.res = HGT%>%t()%>%as.data.frame()
write.table(CelliD.res,'xxx/PaaSc_Data/TMS/TMS_facs_selected_alltraits_CelliD_score.tsv',sep='\t',col.names=T,row.names=T,quote=F)



