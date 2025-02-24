library(Seurat)
library(tidyverse)
library(PaaSc)
library(CelliD)
library(GSDensity)
library(reshape2)

#-----load data-----#
mouse_gene_sets = readRDS('xxx/PaaSc_Data/Traits/mouse_traits_geneset.RDS')
keggs = readRDS('xxx/PaaSc_Data/Traits/mouse_kegg_geneset.RDS')
obj = readRDS('xxx/PaaSc_Data/TMS/TMS_facs_selected.RDS')
gs_anno = read.table('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Traits/gs_anno.txt',sep = '\t', header=T)
ct_anno = read.table('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Traits/ct_anno.txt',sep = '\t', header=T)

#-----Compute Score-----#
# PaaSc
gene_rate = getGeneRate(background.geneset = keggs, pathway.geneset = mouse_gene_sets,mode ='single')
regression_data = doRegression(obj, gene.rate = gene_rate)
PaaSc.res = computeScore(obj, regression.data = regression_data, pvalue = 0.05, weight = FALSE, normalize = "z-score")
write.table(PaaSc.res,'xxx/PaaSc_Data/TMS/TMS_facs_selected_alltraits_PaaSc_score.tsv',sep='\t',col.names=T,row.names=T,quote=F)


# GSDensity
ce <- compute.mca(object = obj)
cells <- colnames(obj)
el <- compute.nn.edges(coembed = ce, nn.use = 300)
GSD <- run.rwr(el = el, gene_set = mouse_gene_sets, cells = cells)
GSD.res <- as.data.frame(GSD)
write.table(GSDensity.res,'xxx/PaaSc_Data/TMS/TMS_facs_selected_alltraits_GSDensity_score.tsv',sep='\t',col.names=T,row.names=T,quote=F)


# CelliD
set.seed(123)
HGT <-RunCellHGT(X = obj,
            reduction = "mca",
            pathways = mouse_gene_sets,
            log.trans = T, minSize = 5,
            dims=1:20
        )
CelliD.res = HGT%>%t()%>%as.data.frame()
write.table(CelliD.res,'xxx/PaaSc_Data/TMS/TMS_facs_selected_alltraits_CelliD_score.tsv',sep='\t',col.names=T,row.names=T,quote=F)


#-----Fisher's Test-----#

# Fisher test after binarization
set.seed(1)
t='UKB_460K.blood_LYMPHOCYTE_COUNT'
label_data = doBinarization(PaaSc.res, n.cluster = 3)
ft = lapply(ct_anno$CellType,function(i){   
    cell = obj@meta.data%>%filter(cell_ontology_class == i) %>%rownames()
    label = factor(label_data[[paste0(t, ".label")]],levels = c('positive','negative'))
    group = factor(ifelse(rownames(label_data) %in% cell,i,'others'),levels = c(i,'others'))
    table_data <- table(label,group)
    test = fisher.test(table_data,alternative = 'greater')
    p_value = test$p.value
    if (is.infinite(test$estimate)) {
        odds_ratio = (table_data[1,1] / table_data[1,2]) / (1 / table_data[2,2])
    } else {
        odds_ratio = test$estimate
    }
    
    res = data.frame(CellType=ct_anno$abbreviation[match(i, ct_anno$CellType)],
                     Pvalue=p_value,
                     Odds_ratio= odds_ratio)
    return(res)
})
df=bind_rows(ft)

write.table(label_data,'xxx/PaaSc_Data/TMS/LYMPHOCYTE_COUNT_binarization.tsv',sep='\t',col.names=T,row.names=T,quote=F)
write.table(df,'xxx/PaaSc_Data/TMS/LYMPHOCYTE_COUNT_fishertest.tsv',sep='\t',col.names=T,row.names=F,quote=F)


# For all traits
fisher_test <- function(method.res){
    Score_all=melt(method.res)
    threshold = quantile(Score_all$value, 0.95)
    #mean score
    df_summary <- method.res[,gs_anno$Traits] %>% mutate(CellType = obj@meta.data$cell_ontology_class)%>%
                  group_by(CellType) %>%
                  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))%>%
                  pivot_longer(col=-1,values_to = 'Mean',names_to = 'Traits')
    #fisher.test
    fisher.res = data.frame()
    for (t in gs_anno$Traits){
        for (i in ct_anno$CellType){
            cell = obj@meta.data%>%filter(cell_ontology_class == i) %>%rownames()
            label = factor(ifelse(PaaSc.res[[t]] > threshold ,'positive','negative') ,levels = c('positive','negative'))
            group = factor(ifelse(rownames(PaaSc.res) %in% cell,i,'others'),levels = c(i,'others'))
            table_data <- table(label,group)
            test = fisher.test(table_data,alternative = 'greater')
            if (is.infinite(test$estimate)) {
               odds_ratio = (table_data[1,1] / table_data[1,2]) / (1 / table_data[2,2])
            } else {
                
              odds_ratio = test$estimate
                
            }
            res = data.frame(CellType=i,Traits=t,Pvalue=test$p.value,Odds_ratio= odds_ratio)
            fisher.res = rbind(fisher.res,res) 
        }
    }
    df <- merge(df_summary, fisher.res[,1:4], by = c("CellType", "Traits"))%>%
            mutate(CellType=factor(ct_anno$abbreviation[match(CellType, ct_anno$CellType)],levels=ct_anno$abbreviation))%>%
            mutate(Traits=factor(gs_anno$abbreviation[match(Traits, gs_anno$Traits)],levels=rev(gs_anno$abbreviation)))%>%
            mutate(FDR=p.adjust(Pvalue,method = 'fdr'))
    return(df)
}

Fig5D_PaaSc_df = fisher_test(PaaSc.res)
Fig5D_GSDensity_df = fisher_test(GSDensity.res)
Fig5D_CelliD_df = fisher_test(CelliD.res)

write.table(Fig5D_PaaSc_df,'xxx/PaaSc_Data/TMS/TMS_facs_selected_alltraits_PaaSc_fishertest.tsv',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Fig5D_GSDensity_df,'xxx/PaaSc_Data/TMS/TMS_facs_selected_alltraits_GSDensity_fishertest.tsv',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Fig5D_CelliD_df,'xxx/PaaSc_Data/TMS/TMS_facs_selected_alltraits_CelliD_fishertest.tsv',sep='\t',col.names=T,row.names=F,quote=F)
