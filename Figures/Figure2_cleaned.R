# Figure 2
# JZ
# 05292021

library(ggpubr)
library(ggsci)
library(reshape2)
library(ggplot2)
library(heatmap3)
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library('RColorBrewer')
library(pheatmap)
library(ggrepel)

# clean environment
rm(list = ls())  

path<-"D:\\TCR\\102018-Neountx\\Result\\neoonly.multilane.cd8.new\\"
ser.integrated<-readRDS("D:\\TCR\\102018-Neountx\\Result\\neoonly.multilane.cd8.new\\ser.integrate.rds")
umap<-readRDS("D:\\TCR\\102018-Neountx\\Result\\neoonly.multilane.cd8.new\\umap.rds")
trab.wide<-readRDS("D:\\TCR\\102018-Neountx\\Data\\Processed\\Scseq\\trab.wide.public.rds")



cell.type<-read.csv("D:\\TCR\\102018-Neountx\\Result\\neoonly.multilane.cd8.new\\celltype.csv")
Idents(ser.integrated)<-'integrated_snn_res.0.5'
new.cluster.ids<-as.character(cell.type$CellType)
names(new.cluster.ids) <-0:13
ser.integrated <- RenameIdents(ser.integrated, new.cluster.ids)
ser.integrated$CellType <- Idents(ser.integrated)


newmeta<-ser.integrated[[]]

newmeta$CellType<-factor(newmeta$CellType,levels=c(
  "Stem-like memory",
  "Effector(I)",
  "Effector(II)",
  "Effector(III)",
  "TRM(I)",
  "TRM(II)",
  "TRM(III)",
  "TRM(IV)",
  "TRM(V)",
  "TRM(VI)",
  "Proliferating",
  "MAIT",
  "CD4CD8(I)",
  "CD4CD8(II)"
))
umap<-as.data.frame(umap)
umap$barcode<-rownames(umap)
newmeta<-newmeta %>% left_join(trab.wide,by=c('barcode','patient_id')) %>% left_join(umap) 
newmeta$antigen[newmeta$antigen=='Viral (CMV, EBV, Influenza A)'&newmeta$Type=='EBV']<-'EBV'
newmeta$antigen[newmeta$antigen=='Viral (Influenza)']<-'InfluenzaA'

rownames(newmeta)<-newmeta$barcode
ser.integrated <- AddMetaData(
  object = ser.integrated,
  metadata = newmeta
)

saveRDS(newmeta,"D:\\Dropbox\\Scdata\\result\\ms\\Fig2\\meta.rds")


#  umap for all T cells from included pts
DimPlot(ser.integrated,reduction = 'umap',label = F)+NoLegend()
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig2\\umap.nolabel.pdf",height =10,width =11)



#################################
#### antigen specific cells #####
#################################
newmeta %>% group_by(antigen,Type) %>% summarise(n=n())
ggplot(newmeta,aes(x=UMAP_1,y=UMAP_2)) +
  geom_point(aes(color=seurat_clusters),alpha=0.2,size=0.2)+
  geom_point(data=newmeta %>% filter(antigen=='InfluenzaA'),col='blue',size=2,shape=16)+
  geom_point(data=newmeta %>% filter(antigen=='EBV'),col='purple',size=2,shape=16)+
  geom_point(data=newmeta %>% filter(antigen=='MANA'),col='red',size=2,shape=17)+
  theme_classic() + 
  theme(plot.title = element_text(size = 30, face = "bold"),
        legend.position = "none",axis.title = element_blank(),
        axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank()) 
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig2\\mana.cef.cluster.colored.pdf",width =7.5, height =7)



#######################################
#### Differential genes by antigen ####
#######################################
ser<-readRDS(paste0(path,"ser.antigen.tumor.without.mana.score.rds"))


##### heatmap for MANA vs EBV vs Flu ####
group.marker <- FindAllMarkers(ser, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
group.marker.sig<-group.marker %>% filter(p_val_adj<0.05) 
write.csv(group.marker.sig,paste0(path,'mana.ebv.flu.tumor.csv'),row.names = F)

ser[[]] %>% group_by(antigen) %>% summarise(n=n())
###  volcano plot ###
plot<-group.marker %>% filter(p_val_adj<0.05&abs(avg_logFC)>0.8) 

### heatmap of differential gene ###
expr<-ser@assays$RNA@data
colnames(expr)<-paste0(sub('-.*','',Idents(ser)),'_',
                       ser$TRM_gene_set_score1,'_',colnames(expr))

# add selective genes
marker.genes<-c('EOMES','TBX21','TOX','TOX2','PRDM1','HLA-DRA','GZMB','GZMH','GZMK','NKG7','TNFRSF9','IFNG','ZNF683','ITGAE','TCF7','IL7R','PDCD1','CTLA4','ENTPD1','LAG3','TIGIT','HAVCR2')
marker.genes.sig<-intersect(marker.genes,group.marker.sig$gene)
genes<-unique(c(marker.genes.sig,plot$gene))
mat<-expr[rownames(expr) %in% genes,]
#mat<-expr[rownames(expr) %in% plot$gene,]

mat<-as.matrix(mat)
mat<-mat[, order(colnames(mat))]


side = do.call(rbind,strsplit(colnames(mat),"_")) %>% data.frame()
rownames(side)<-colnames(mat)
side<-side[,1:2]
colnames(side)<-c('antigen','TRM.signature')
c1 <- c('purple','blue','red')
names(c1) <- c(unique(side$antigen))
side$TRM.signature<-as.numeric(side$TRM.signature)

library(pheatmap)
pdf("D:\\Dropbox\\Scdata\\result\\ms\\Fig2\\genes_antigen_add_marker.pdf", width =6, height = 6)
pheatmap(mat,cluster_rows = T,
         color=colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols = T,
         show_colnames =F,
         annotation_col = side,
         cutree_rows = 4,
         annotation_colors = list(antigen=c1[1:3]),
         scale = 'row')

dev.off()


### flowplot for selective gene panel ###
RidgePlot(ser,features = c('ZNF683','ITGAE','TCF7','IL7R','EOMES','TBX21','TOX','TOX2'
                           ,'PRDM1','PDCD1','CTLA4','ENTPD1','LAG3','TIGIT','HAVCR2',
                           'HLA-DRA','GZMB','GZMK','NKG7','TNFRSF9','IFNG'),ncol=6,
          cols=c('red','blue','purple'),slot = 'counts')
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig2\\flowplot.genes.pdf",width=22,height=12)





#### flu vs MANA specific T cells ####
group.marker <- FindMarkers(ser, only.pos = F,ident.1 = 'MANA',ident.2 = 'InfluenzaA', min.pct = 0.25, logfc.threshold = 0.25)
group.marker$gene<-rownames(group.marker)
group.marker.sig<-group.marker %>% filter(p_val_adj<0.05) 

mana<-group.marker.sig %>% top_n(15,wt = avg_logFC) 
flu<-group.marker.sig %>% top_n(15,wt = -avg_logFC) 
mana$antigen<-"Enriched in MANA"
flu$antigen<-"Enriched in flu"

plotdata<-rbind(mana,flu)


# plot difference  
ggbarplot(plotdata, x = "gene", y = "avg_logFC",
          fill = "antigen",          
          color = "antigen",            
          palette = "jco",            
          sort.val = "asc",         
          sort.by.groups = T,    
          x.text.angle = 45,         
          ylab = "",
          xlab = "",
          legend.title = "",
          rotate = F
          
)
ggsave('D:\\Dropbox\\Scdata\\result\\ms\\Fig2\\flu.mana.pdf',width = 5,height = 3.5)







### IL7 functional experiment  ###

ser.saver.antigen<-readRDS('D:\\TCR\\102018-Neountx\\Result\\IL7R.multilane\\ser.saver.antigen.rds')

pub1<-c(
  'IKZF4',
  'SYNE3',
  'CCR5',
  'ETS1',
  'BHLHE40',
  'UGCG',
  'FAM101B',
  'ST3GAL5',
  'SLC37A3',
  'IRF4',
  'AFAP1',
  'SLC4A10',
  'IGFBP3',
  'FKBP5',
  'PITRM1',
  'CYLD',
  'KLF7',
  'ADAM19',
  'TXK',
  'MB21D2',
  'AHR',
  'TAF4B',
  'FAM13A',
  'PTGER2',
  'BCL2',
  'CISH',
  'AP3M2',
  'PDE4B',
  'DPP4',
  'CD8B',
  'CMTM6',
  'GNPDA1',
  'CMAHP',
  'ITGA4',
  'SOCS2',
  'TLR1',
  'CCR2',
  'LTA',
  'RGS1',
  'SOS1',
  'CDK6',
  'IL2RA',
  'FRMD4B',
  'MEOX1',
  'HSPA1L',
  'IFNG')


genelist<-list(pub1)
ser.saver.antigen <- AddModuleScore(
  object = ser.saver.antigen,
  features = genelist,
  ctrl = 5,
  name = 'IL7.gene.list'
)


meta<-ser.saver.antigen[[]]
flu<-meta$barcode[grepl('Flu',meta$sample.n)]
ma<-meta$barcode[grepl('MANA',meta$sample.n)]


stat<-ser.saver.antigen[[]]
se <- function(x, ...) sqrt(var(x, ...)/length(x))
stat.plot<-stat  %>% group_by(patient_id,antigen,pep,dose) %>% 
  summarise(mean=mean(IL7.gene.list1), median=median(IL7.gene.list1),
            se=se(IL7.gene.list1),n=n())

pd <- position_dodge(0.1)
stat.plot$antigen.condition[(stat.plot$antigen=="MANA"&stat.plot$pep=='MANA')]<-'MANA specific T cell+MANA peptide'
stat.plot$antigen.condition[(stat.plot$antigen=="MANA"&stat.plot$pep=='Flu')]<-'MANA specific T cell+Flu peptide'
stat.plot$antigen.condition[(stat.plot$antigen=="InfluenzaA"&stat.plot$pep=='MANA')]<-'Flu specific T cell+MANA peptide'
stat.plot$antigen.condition[(stat.plot$antigen=="InfluenzaA"&stat.plot$pep=='Flu')]<-'Flu specific T cell+Flu peptide'
stat.plot$match[stat.plot$antigen.condition %in% c('MANA specific T cell+MANA peptide','Flu specific T cell+Flu peptide')]<-'Antigen matched'
stat.plot$match[stat.plot$antigen.condition %in% c('MANA specific T cell+Flu peptide','Flu specific T cell+MANA peptide')]<-'Antigen mismatched'

ggplot(stat.plot, aes(x=factor(dose), y=mean, colour=antigen, group=antigen.condition)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd,size=2,aes(linetype=match)) +
  geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
  xlab("") +
  ylab("") +
  theme_classic() +
  scale_color_manual(name="Antigen specific T cells",    # Legend label, use darker colors
                     breaks=c("MANA", "InfluenzaA"),
                     labels=c("MANA", "Influenza A"),values = c('red','blue'))+
  theme(axis.text.x = element_text(angle =45, hjust = 1,size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        strip.text.x = element_text(size = 10,face = "bold"),
        axis.title.x =element_text(size = 14,face = "bold"),
        axis.title.y=element_text(size = 14,face = "bold"),
        legend.position = 'none')
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig2\\IL7.gene.score.pdf",width=5,height = 3)



