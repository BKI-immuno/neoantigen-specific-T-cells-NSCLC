# Figure 4
# 05292021




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



dir.create("D:\\TCR\\102018-Neountx\\Result\\blood\\")
path<-"D:\\TCR\\102018-Neountx\\Result\\blood\\"
annotate<-read.csv("D:\\TCR\\102018-Neountx\\Data\\Processed\\Scseq\\annotate.new.dge.csv")
annotate$sample.n<-paste0(annotate$patient_id,':',annotate$tissue,'-',annotate$batch)
ser.count<-readRDS("D:\\TCR\\102018-Neountx\\Data\\Processed\\Scseq\\ser.count.rds")
cells<-readRDS("D:\\TCR\\102018-Neountx\\Data\\Processed\\Scseq\\cells.pass.qc.rds")
qc.meta<-readRDS("D:\\TCR\\102018-Neountx\\Data\\Processed\\Scseq\\qc.rds")
qc.meta<-qc.meta %>% left_join(unique(annotate[,c('sample','sample.n')]))
rownames(qc.meta)<-qc.meta$barcode
file<-annotate %>% filter(blood.mana!=0)

pbmc<-ser.integrated@assays$RNA@counts

hobit.cluster<-subset(ser.integrated,subset=integrated_snn_res.0.3==5)

hobit<-hobit.cluster[[]] %>% group_by(TRB_aa_1,antigen) %>% summarise(n=n())


# genes to be excluded from clustering
ifn<-read.csv("D:\\TCR\\102018-Neountx\\Data\\Raw\\Published gene list\\ifn.csv")
ig<-read.csv("D:\\TCR\\102018-Neountx\\Data\\Raw\\Published gene list\\ig.csv")
mito<-read.csv("D:\\TCR\\102018-Neountx\\Data\\Raw\\Published gene list\\mito.csv")
cluster.exclude<-c(ifn$IFN_module,ig$Ig_module,mito$mito)

# include all cells passed qc 

pbmc.cells<-readRDS("D:\\TCR\\102018-Neountx\\Data\\Processed\\Scseq\\pbmc.include.rds")
cd8<-readRDS("D:\\TCR\\102018-Neountx\\Data\\Processed\\Scseq\\DGE\\saver\\cutoff\\CD8A.cutoff.rds")
cells<-intersect(pbmc.cells$barcode,cd8)


ser.count<-ser.count[,sub('_.*','',colnames(ser.count)) %in% file$sample&colnames(ser.count) %in% cells]

ser.count$orig.ident<-sub('_.*','',colnames(ser.count))

ser.count <- AddMetaData(
  object = ser.count,
  metadata = qc.meta
)

ser.list<-SplitObject(ser.count,split.by = 'sample.n')
rm(ser.count)
#rm(ser)

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

anchors <- FindIntegrationAnchors(object.list = ser.list,
                                  dims = 1:30)  # specify anchoring dataset here

ser.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)


ser.integrated <- ScaleData(ser.integrated, verbose = FALSE) #vars.to.regress = c("cc_score1","S.Score","G2M.Score")
DefaultAssay(ser.integrated)<-"integrated"

ser.integrated <- RunPCA(object = ser.integrated,npcs =30,verbose = FALSE,
                         features =setdiff(VariableFeatures(object = ser.integrated),cluster.exclude))


ser.integrated<-RunTSNE(object = ser.integrated,npcs =30,verbose = FALSE,
                        features =setdiff(VariableFeatures(object = ser.integrated),cluster.exclude))

ElbowPlot(ser.integrated, ndims = 100)
print(ser.integrated[["pca"]], dims = 1:30, nfeatures = 5)
ser.integrated <- FindNeighbors(ser.integrated, reduction = "pca", dims = 1:30, nn.eps = 0.5)

ser.integrated <- FindClusters(ser.integrated, resolution = 0.6, n.start = 10)
ser.integrated <- FindClusters(ser.integrated, resolution = 0.5, n.start = 10)
ser.integrated <- FindClusters(ser.integrated, resolution = 0.3, n.start = 10)
ser.integrated <- FindClusters(ser.integrated, resolution = 0.7, n.start = 10)
ser.integrated <- RunUMAP(object = ser.integrated,reduction='pca',dims = 1:30)


top3000 <- head(VariableFeatures(ser.integrated), 3000) 
dir.create(paste0(path,'/'))
saveRDS(top3000,paste0(path,"/top3000.hvg.rds"))


FeatureScatter(ser.integrated, feature1 = "PC_1", feature2 = "PC_2",
               group.by = 'integrated_snn_res.0.3')


# find diffrential markers
ser.integrated.small<-subset(ser.integrated,downsample=10000)
markers <- FindAllMarkers(ser.integrated.small, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
top10<-markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top20<-markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
tail20<-markers %>% group_by(cluster) %>% top_n(n = 20, wt = -avg_logFC)
write.csv(markers,paste0(path,"/markers.rawcount.csv"),row.names = F)
write.csv(top20,paste0(path,"/top20.rawcount.csv"),row.names = F)
write.csv(tail20,paste0(path,"/tail20.rawcount.csv"),row.names = F)
write.csv(top10,paste0(path,"/top10.rawcount.csv"),row.names = F)


rownames(umap)<-umap$barcode
ser.integrated <- AddMetaData(
  object = ser.integrated,
  metadata = umap
)

ser.integrated$IL7R<-ser.integrated@assays$RNA@data['IL7R',]
trab.wide<-readRDS("D:\\TCR\\102018-Neountx\\Data\\Processed\\Scseq\\trab.wide.public.rds")
ser.integrated<-readRDS(paste0(path,"/ser.integrate.rds"))
umap<-readRDS(paste0(path,"/umap.rds"))
umap<-as.data.frame(umap)
umap$barcode<-rownames(umap)


meta <-ser.integrated[[]]
umap<-umap %>% left_join(meta) %>% left_join(trab.wide) 

umap$tissue<-factor(umap$tissue,levels=c('W2','W4','M3'))

ggplot(umap,aes(x=UMAP_1,y=UMAP_2)) +
  geom_point(aes(color=integrated_snn_res.0.3),alpha=0.5,size=1)+
   theme_classic() + 
  scale_color_npg()+
  theme( legend.position = "none",axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank()) 

ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig4\\cluster.pdf",width = 5,height = 3.5)


ser.integrated.small<-subset(ser.integrated,downsample=100)
DefaultAssay(ser.integrated)<-'RNA'
markers<-read.csv(paste0(path,"/markers.rawcount.csv"))
top5<-markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_logFC)

hist(top5$avg_logFC,breaks=100)
DefaultAssay(ser.integrated.small)<-'integrated'
pdf("D:\\Dropbox\\Scdata\\result\\ms\\Fig4\\top4.gene.",width = 5,height = 7)
DoHeatmap(ser.integrated.small, raster = F,features = top5$gene,angle=90,group.by = 'integrated_snn_res.0.3') 
dev.off()
dggsave(paste0(path, "/top4.gene.0.3.tiff"),height = 7,width =5,dpi = 300)
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig4\\top4.gene.pdf",width = 5,height = 7)


write.csv(top5,paste0(path,"/top5.csv"),row.names = F)
gene<-read.csv(paste0(path,"/gene.show.csv"))
DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.3' ,label = F)
ggsave(paste0(path, "/cluster.0.3.reduced.tiff"),height = 6,width =7,dpi = 300)




ggplot(umap,aes(x=UMAP_1,y=UMAP_2)) +
  geom_point(aes(color=integrated_snn_res.0.3),alpha=0.3,size=1)+
  geom_point(data=umap %>% filter(antigen %in% c('MANA','mana_score')),aes(col=integrated_snn_res.0.3),
             size=1.5,shape=17)+
  geom_point(data=umap %>% filter(antigen %in% c('MANA','mana_score')),
             shape=2,size=1.5,col='black')+
  theme_classic() + 
  scale_color_npg()+
  facet_wrap(~tissue)+
  theme(plot.title = element_text(size = 30, face = "bold"),
        legend.position = "none",axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank()) 
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig4\\all.mana.pdf",width = 8,height =2.5 )




umap %>% filter(antigen %in% c('MANA','mana_score')) %>% group_by(tissue,integrated_snn_res.0.3) %>%
  summarise(n=n()) %>% group_by(tissue) %>% mutate(total=sum(n),p=n/total) %>%
  #  filter(integrated_snn_res.0.6==6) %>%
  ggplot(aes(x=tissue,y=p,fill=integrated_snn_res.0.3))+
  #geom_bar(stat = 'identity',width=1,position = position_dodge(preserve = "single")) +
  geom_bar(stat = 'identity', color = "black") +
  # geom_jitter(position=position_dodge(0.8))+
  theme_classic()+
  scale_fill_manual(values = c("#DC0000B2","#00A087B2", "#4DBBD5B2","#8491B4B2"),
                    breaks = c("0", "2", "1","5"),name='')+
  scale_y_continuous(labels = scales::percent)+
  labs(y='',x='')+
  #  stat_compare_means(label ='p.signif')+
  theme(axis.text.x = element_text(angle =40, hjust = 1,size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        strip.text.x = element_text(size = 14,face = "bold"),
        axis.title.x =element_text(size = 14,face = "bold"),
        axis.title.y=element_text(size = 12,face = "bold"),
        legend.position = 'none'
        
  )
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig4\\prop.pdf",width = 3,height =2 )

