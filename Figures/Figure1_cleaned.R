
# 11192020
# JZ

library(ggpubr)
library(ggsci)
library(reshape2)
library(ggplot2)
library(heatmap3)
library(Seurat)
library(tidyverse)

# clean environment
rm(list = ls())  

##################
### load data ####
##################
path<-"D:\\TCR\\102018-Neountx\\Result\\neoonly.multilane.new\\"
ser.integrated<-readRDS("D:\\TCR\\102018-Neountx\\Result\\neoonly.multilane.new\\ser.integrate.rds") 
umap<-ser.integrated@reductions$umap@cell.embeddings
trab.wide<-readRDS("D:\\TCR\\102018-Neountx\\Data\\Processed\\Scseq\\trab.wide.public.rds")


# remove cluster with 2 cells 
test<-ser.integrated[[]]
test %>% group_by(seurat_clusters) %>% summarise(n=n()) %>% arrange(n)
test<-test %>% filter(seurat_clusters!=16&seurat_clusters!=15)
ser.integrated<-subset(ser.integrated,subset=barcode %in% test$barcode)


### annotate clusters ###
cell.type<-read.csv("D:\\TCR\\102018-Neountx\\Result\\neoonly.multilane.new\\celltype.csv")
new.cluster.ids<-as.character(cell.type$CellType)
names(new.cluster.ids) <- levels(ser.integrated)
ser.integrated <- RenameIdents(ser.integrated, new.cluster.ids)
ser.integrated$CellType <- Idents(ser.integrated)
DefaultAssay(ser.integrated)<-'integrated'


umap<-as.data.frame(umap)
umap$barcode<-rownames(umap)
newmeta<-ser.integrated[[]]
unique(newmeta$CellType)
newmeta$CellType<-factor(newmeta$CellType,levels = c("CD8-Effector(1)",
                                                     "CD8-Effector(2)",
                                                     "CD8-Effector(3)",
                                                     "CD8-TRM(1)",
                                                     'CD8-TRM(2)',
                                                     "CD8-Proliferating",
                                                     "CD8-MHCII",
                                                     "MAIT",
                                                     'Stem-like memory',
                                                     "CD4-Th(1)", 
                                                     "CD4-Th(2)",
                                                     "CD4-Th(3)",
                                                     "CD4-Tfh(1)",
                                                     "CD4-Tfh(2)",
                                                     "CD4-Treg"))

newmeta<-newmeta %>% left_join(umap)

trab.wide<-readRDS("D:\\TCR\\102018-Neountx\\Data\\Processed\\Scseq\\trab.wide.public.rds")
newmeta<-newmeta %>% left_join(trab.wide)
rownames(newmeta)<-newmeta$barcode


ser.integrated <- AddMetaData(
  object = ser.integrated,
  metadata = newmeta
)



### Figure 1 ####
# Fig1b
DimPlot(ser.integrated,reduction = 'umap',group.by = 'CellType',label = F,pt.size =0.01)+NoLegend()
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig1\\umap.nolabel.pdf",height =10,width =11)

DimPlot(ser.integrated,reduction = 'umap',group.by = 'CellType',repel = T,label = T)
ggsave(paste0(path, "umap.dimplot.lab.pdf"),height = 5,width =8,dpi = 300)
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig1\\umap.pdf",height =5,width =8)


# Fig1c
ser.integrated.small<-subset(ser.integrated,downsample=5000)
markers<-read.csv(paste0(path,"/markers.rawcount.csv"))
top5<-markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_logFC)  
write.csv(top5,"D:\\Dropbox\\Scdata\\result\\top5.csv",row.names = F)
top5<-read.csv("D:\\Dropbox\\Scdata\\result\\fig1.gene.show.csv")

DefaultAssay(ser.integrated.small)<-'integrated'
DoHeatmap(ser.integrated.small,group.by="CellType",raster = F ,features = top5$gene,angle=90)
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig1\\genes.rawcount.heatmap.top5.with.legend.pdf", width =4, height =6)

DoHeatmap(ser.integrated.small,group.by="CellType",raster = F ,features = top5$gene,angle=90) +NoLegend()
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig1\\genes.rawcount.heatmap.top5.without.legend.pdf", width =4, height =6)




# Figure 1D expression profile by genes 
DefaultAssay(ser.integrated.small)<-'RNA'
genes<-c("CD8A","CD4","FOXP3", 'GZMK','ITGAE',
         'ZNF683', 'TCF7','CXCL13', 'SLC4A10',"PDCD1", 'CTLA4',"HAVCR2","TIGIT",
         'ENTPD1','LAG3' )

p3 <- FeaturePlot(ser.integrated.small, features = genes,order=T,
                  reduction = "umap",pt.size = 1,  cols = c("grey", "red"),combine = FALSE) 
p3 <- lapply(X = p3, FUN = function(x) AugmentPlot(x + NoLegend()))
CombinePlots(plots = p3,ncol=3)
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig1\\genes.pdf",height =15,width =9)


# PCA
dt = readRDS(paste0("D:\\Dropbox\\Scdata\\result\\pca\\pca_patienttissue\\all\\dat_hvg_m2.rds")) # read in file
dt$response[dt$resi_tumor>0.1]<-'NR'
dt$response[is.na(dt$response)==T]<-'R'
min = min(dt$PC1,dt$PC2)
max = max(dt$PC1,dt$PC2)


ggplot(dt,aes(x=PC1,y=PC2,color =response)) +
  geom_point(size=3) +
  geom_point(shape = 1,size = 3,colour = "black")+
  theme_classic() +
  scale_color_npg()+
  xlim(min,max) + ylim(min,max) +
  theme(axis.text.x = element_text(angle =45, hjust = 1,size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        strip.text.x = element_text(size = 14,face = "bold"),
        axis.title =element_blank(),legend.position = 'none')
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig1\\pca.response.pdf",height =3,width =3.5)


ggplot(dt, aes(x=PC1,y=PC2,color =tissue)) +
  geom_point(size=3) +
  geom_point(shape = 1,size = 3,colour = "black")+
  theme_classic() +
  scale_color_jco()+
  xlim(min,max) + ylim(min,max)+
  theme(axis.text.x = element_text(angle =45, hjust = 1,size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        strip.text.x = element_text(size = 14,face = "bold"),
        axis.title =element_blank(),legend.position = 'none')
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig1\\pca.tissue.pdf",height =3,width =3.5)




