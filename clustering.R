

rm(list = ls()) # clean environment

library(Seurat)
library(data.table)
library(scales)
library(ggbeeswarm)
library(tidyverse)
library(gplots)
library(dplyr)
library(ggpubr)
library(ggsci)
library(reshape2)
library(ggplot2)
library(heatmap3)

###### WZ 04/16/2020 #######
####get variable feature####
library(cowplot)
library(matrixStats)
library(ggpointdensity)


hypervar <- function(data, span = 0.5, showplot = TRUE, font_size = 14){
  
  gene_mean_all <- rowMeans(data)
  gene_var_all <- rowVars(data)
  
  data_filter <- data[gene_mean_all > 0 & gene_var_all > 0,]
  
  gene_mean <- rowMeans(data_filter)
  gene_var <- rowVars(data_filter)
  
  data_fit <- data.frame(X=gene_mean,Y=gene_var)
  fit_model <- loess(formula = log2(x=Y) ~ log2(x=X),
                     data = data_fit,
                     span = span)
  
  gene_var_expect <- fit_model$fitted
  gene_hyper_var <- log2(gene_var) - gene_var_expect
  
  result <- data.frame(feature=row.names(data_filter), mean=gene_mean, var=gene_var,
                       var_expect_log2=gene_var_expect,hypervar_log2=gene_hyper_var)
  
  p1 <- ggplot(result, aes(log2(mean), log2(var))) +
    geom_pointdensity() +
    scale_color_viridis_c() +
    geom_point(data=result,aes(log2(mean),var_expect_log2),color="red") +
    theme_bw() +
    theme(axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))
  
  p2 <- ggplot(result, aes(log2(mean), hypervar_log2)) +
    geom_pointdensity() +
    scale_color_viridis_c() +
    theme_bw() +
    theme(axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))
  
  combined_plot <- plot_grid(p1, p2, labels = c('A', 'B'))
  
  if(showplot){
    print(combined_plot)
  }
  
  return(list(data=result,plot=combined_plot))
}


theme_set(theme_classic())

dir.create("path")
path<-"path"
ser.count<-readRDS("D:\\TCR\\102018-Neountx\\Data\\Processed\\Scseq\\ser.count.rds") # all cells passed qc
# genes to be excluded from clustering
ifn<-read.csv("D:\\TCR\\102018-Neountx\\Data\\Raw\\Published gene list\\ifn.csv")
ig<-read.csv("D:\\TCR\\102018-Neountx\\Data\\Raw\\Published gene list\\ig.csv")
mito<-read.csv("D:\\TCR\\102018-Neountx\\Data\\Raw\\Published gene list\\mito.csv")
cluster.exclude<-c(ifn$IFN_module,ig$Ig_module,mito$mito)

ser<-CreateSeuratObject(ser.count,min.cells =5,min.features = 250,meta=meta)



ser.list<-SplitObject(ser,split.by = 'sample')

for (i in 1:length(ser.list)) {
  ser.list[[i]] <- NormalizeData(ser.list[[i]] , verbose = FALSE)
  
  ###### WZ 04/16/2020 #######
  #ser.list[[i]]  <- FindVariableFeatures(ser.list[[i]] , selection.method="vst",nfeatures=3000,verbose = FALSE)
  #########get variable feature##########
  gene_count <- as.matrix(ser.list[[i]]@assays$RNA@counts)
  gene_count_norm <- sweep(gene_count,2,colSums(gene_count),FUN="/")*1e4
  gene_hypervar <- hypervar(gene_count_norm,showplot=FALSE)
  gene_hypervar_sort <- gene_hypervar$data %>% arrange(.,desc(hypervar_log2))
  VariableFeatures(ser.list[[i]]) <- gene_hypervar_sort$feature[1:3000]
  VariableFeatures(ser.list[[i]])<-setdiff(VariableFeatures(ser.list[[i]]),cluster.exclude)
  #######################################
}



###### integrate data #######
# use sample with top N as reference group
temp<-ser[[]]
top.sample<-temp %>% group_by(sample.n,tissue,response) %>% summarise(n=n()) %>%
  group_by(tissue,response) %>% top_n(1,wt=n)
ref<-which(names(ser.list) %in% top.sample$sample.n)
ref<-as.numeric(ref)

ser.list[ref]

anchors <- FindIntegrationAnchors(object.list = ser.list,reference=ref,
                                  dims = 1:30)  # specify anchoring dataset here

ser.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)


ser.integrated <- ScaleData(ser.integrated, verbose = FALSE) #vars.to.regress = c("cc_score1","S.Score","G2M.Score")
DefaultAssay(ser.integrated)<-"integrated"

ser.integrated <- RunPCA(object = ser.integrated,npcs =30,verbose = FALSE,
                         features =setdiff(VariableFeatures(object = ser.integrated),cluster.exclude))

ElbowPlot(ser.integrated, ndims = 100)
print(ser.integrated[["pca"]], dims = 1:30, nfeatures = 5)
ser.integrated <- FindNeighbors(ser.integrated, reduction = "pca", dims = 1:30, nn.eps = 0.5)
ser.integrated <- FindClusters(ser.integrated, resolution = 0.7, n.start = 10)
ser.integrated <- FindClusters(ser.integrated, resolution = 0.6, n.start = 10)
ser.integrated <- FindClusters(ser.integrated, resolution = 0.5, n.start = 10)
ser.integrated <- FindClusters(ser.integrated, resolution = 0.3, n.start = 10)
ser.integrated <- RunUMAP(object = ser.integrated,reduction='pca',dims = 1:30)

top3000 <- head(VariableFeatures(ser.integrated), 3000) 

# find diffrential markers
ser.integrated.small<-subset(ser.integrated,downsample=5000)
markers <- FindAllMarkers(ser.integrated.small, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
top10<-markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top20<-markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
tail20<-markers %>% group_by(cluster) %>% top_n(n = 20, wt = -avg_logFC)
write.csv(markers,paste0(path,"/markers.rawcount.csv"),row.names = F)




DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.3' ,label = F)
ggsave(paste0(path, "/cluster.0.3.tiff"),height = 6,width =7,dpi = 300)

DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.5' ,label = F)
ggsave(paste0(path, "/cluster.0.5.tiff"),height = 6,width =7,dpi = 300)

DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.7' ,label = F)
ggsave(paste0(path, "/cluster.0.7.tiff"),height = 6,width =7,dpi = 300)


DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.3' ,label = T)
ggsave(paste0(path, "/cluster.0.3.label.tiff"),height = 6,width =7,dpi = 300)

DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.5' ,label = T)
ggsave(paste0(path, "/cluster.0.5.label.tiff"),height = 6,width =7,dpi = 300)

DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.6' ,label = T)
ggsave(paste0(path, "/cluster.0.6.label.tiff"),height = 6,width =7,dpi = 300)


DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.7' ,label = T)
ggsave(paste0(path, "/cluster.0.7.label.tiff"),height = 6,width =7,dpi = 300)


cluster<-Idents(ser.integrated)
umap<-ser.integrated@reductions$umap@cell.embeddings

saveRDS(top10, paste0(path,"/top10.rds"))
saveRDS(ser.integrated,paste0(path,"/ser.integrate.rds"))
saveRDS(cluster,paste0(path,"/cluster.rds"))
saveRDS(umap, paste0(path,"/umap.rds"))
