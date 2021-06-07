# Figure 3
# R vs NR


library(ggpubr)
library(ggsci)
library(reshape2)
library(ggplot2)
library(heatmap3)
library(Seurat)
library(tidyverse)
library(pheatmap)

dat<-readRDS("D:\\TCR\\102018-Neountx\\Data\\Processed\\Scseq\\DGE\\saver\\saver.rds")
newmeta<-readRDS("D:\\TCR\\102018-Neountx\\Result\\neoonly.multilane.cd8.new\\meta.rds")
validated<-read.csv("D:\\TCR\\102018-Neountx\\Data\\Raw\\Antigen\\mana_11.13.2020.csv")



an<-newmeta %>% filter(antigen %in% c('mana_score','MANA')&tissue %in% c('tumor'))
an %>% group_by(antigen) %>% summarise(n=n())
an %>% group_by(antigen,TRB_aa_1,patient_id) %>% summarise(n=n()) %>% group_by(antigen) %>% summarise(n.clones=n())
dat1<-dat[,colnames(dat) %in% an$barcode]
#dat2<-dat[,colnames(dat) %in% test.include$barcode]
rownames(newmeta)<-newmeta$barcode
ser<-CreateSeuratObject(dat1,meta.data = newmeta)
ser<-NormalizeData(ser, verbose = FALSE)
ser <- ScaleData(ser, verbose = FALSE) #vars.to.regress = c("cc_score1","S.Score","G2M.Score")


# antigen type
Idents(ser) <- "response"
group.marker <- FindMarkers(ser, only.pos = F,ident.1 = 'R',ident.2 = 'NR',logfc.threshold = 0.25)
group.marker$gene<-rownames(group.marker)
group.marker.sig<-group.marker %>% filter(p_val_adj<0.05) 
write.csv(group.marker.sig,paste0(path,'mana.r.nr.tumor.add.mana.score.csv'),row.names = F)
write.csv(group.marker,paste0(path,'mana.r.nr.tumor.all.csv'),row.names = F)



plot<-group.marker %>% filter(p_val_adj<0.05&abs(avg_logFC)>0.8) 
expr<-ser@assays$RNA@data
# add panel genes 
marker.genes<-c('EOMES','TBX21','TOX','TOX2','PRDM1','HLA-DRA','GZMB','GZMH','GZMK','NKG7','TNFRSF9','IFNG','ZNF683','ITGAE','TCF7','IL7R','PDCD1','CTLA4','ENTPD1','LAG3','TIGIT','HAVCR2')
marker.genes.sig<-intersect(marker.genes,group.marker.sig$gene)
genes<-unique(c(marker.genes.sig,plot$gene))
mat<-expr[rownames(expr) %in% genes,]
mat<-as.matrix(mat)
mat<-mat[, order(colnames(mat))]


### heatmap of differential gene ###
colnames(mat)<-paste0(ser$response,'_',ser$TRB_aa_1,'_',
                      ser$imid.x,'_',colnames(mat))


side = do.call(rbind,strsplit(colnames(mat),"_")) %>% data.frame()
rownames(side)<-colnames(mat)
side<-side[,1:3]
colnames(side)<-c('response','Clonotype','patientID')
c1 <- c("#E64B35B2","#4DBBD5B2","#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2",
         "#91D1C2B2", "#DC0000B2" )
names(c1) <- c(unique(side$response),unique(side$patientID))



my_gene_col<-group.marker %>% filter(gene %in% genes) %>% select(p_val_adj)
my_gene_col$p_val_adj<--log10(my_gene_col$p_val_adj)
names(my_gene_col)<-'-log10(FDR)'


pdf("D:\\Dropbox\\Scdata\\result\\ms\\Fig3\\mprvsnonmpr.pdf", width =7, height = 7)
pheatmap(mat,cluster_rows = T,
         color=colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols = T,
         show_colnames =F,
         annotation_col = side,
         annotation_row = my_gene_col,
         cutree_rows = 2,
         annotation_colors = list(response=c1[1:2],patientID=c1[3:8]),
         scale = 'row'
)
dev.off()

tiff("D:\\Dropbox\\Scdata\\result\\ms\\Fig3\\RvsNR.gene.add.mamascore.tiff", width =8, height = 10, units = 'in', res = 300)
pheatmap(mat,cluster_rows = T,
         color=colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols = T,
         show_colnames =F,
         annotation_col = side[c(1)],
         annotation_row = my_gene_col,
         cutree_rows = 2,
         annotation_colors = list(response=c1[1:2],Clone.type=c1[3:5]),
         scale = 'row'
)
dev.off()


ser<-subset(ser,subset=tissue %in% c('tumor')&antigen %in% c('MANA','mana_score')&patient_id!='JS10')
Idents(ser)<-'response'


RidgePlot(ser,features = c( 'ZNF683','ITGAE','TCF7','IL7R','HLA-DRA','GZMB','PDCD1','CTLA4','GZMK','NKG7','ENTPD1','LAG3','TNFRSF9','IFNG',
                               'TIGIT','HAVCR2','EOMES','TBX21','TOX','TOX2','PRDM1'),
          ncol=4,cols=c("#E64B35B2","#4DBBD5B2"),slot = 'counts')
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig3\\flowplot.genes.pdf",width=10,height=15,dpi=300)




# antigen specific T cells: MANA vs EBV vs Flu
newmeta$response[is.na(newmeta$response)==T]<-'NR'
an<-newmeta %>% filter(antigen %in% c('MANA','mana_score')&tissue %in% c('tumor'))
an %>% group_by(antigen,tissue,response) %>% summarise(n=n())
an %>% group_by(antigen,TRB_aa_1,validated) %>% summarise(n=n()) %>% group_by(antigen,validated) %>% summarise(n.clones=n())

genes<-c('IL7R','LAYN')
expr<-t(as.matrix(dat[genes,colnames(dat) %in% an$barcode]))
expr<-as.data.frame(expr)
expr$barcode<-rownames(expr)
dat.long<-gather(expr, genes, expression, IL7R:LAYN, factor_key=TRUE)

library(ggpubr)
dat.long.meta<-dat.long %>% left_join(newmeta)
ggplot(dat.long.meta %>% filter(genes=='IL7R') ,
       aes(x=response, y= expression,fill=response)) +
  geom_violin(position=position_dodge(0.8)) +
  labs(y="",x="")+
  theme_classic()+
  scale_fill_npg()+
  stat_compare_means(label.y = 4.5)+
  theme(axis.text.x =  element_blank(),
        legend.position = 'none',
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x =element_text(size = 14,face = "bold"),
        axis.title.y=element_text(size = 14,face = "bold"))
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig3\\IL7R.cell.response.mana.cells.pdf",height = 4,width =4)


#### T cell immunce checkpoint score #####

an<-newmeta %>% filter(antigen %in% c('mana_score','MANA')&
                         tissue %in% c('tumor'))
an %>% group_by(antigen) %>% summarise(n=n())
an %>% group_by(antigen,TRB_aa_1,patient_id) %>% summarise(n=n()) %>% group_by(antigen) %>% summarise(n.clones=n())
dat1<-dat[,colnames(dat) %in% an$barcode]
rownames(newmeta)<-newmeta$barcode
ser<-CreateSeuratObject(dat1,meta.data = newmeta)
ser<-NormalizeData(ser, verbose = FALSE)
ser <- ScaleData(ser, verbose = FALSE) #vars.to.regress = c("cc_score1","S.Score","G2M.Score")

ex<-c('PDCD1','CTLA4','LAG3','HAVCR2','TIGIT','ENTPD1')

genelist<-list(as.character(ex))
ser <- AddModuleScore(
  object = ser,
  features = genelist,
  ctrl = 5,
  name = 'immune.checkpoint.score'
)

plotdata<-ser[[]]
ggplot(plotdata ,
         aes(x=response, y= immune.checkpoint.score1,fill=response)) +
  geom_violin(position=position_dodge(0.8)) +
  labs(y="",x="")+
  theme_classic()+
  scale_fill_npg()+
  stat_compare_means(label.y = 1)+
  theme(axis.text.x =  element_blank(),
        legend.position = 'none',
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x =element_text(size = 14,face = "bold"),
        axis.title.y=element_text(size = 14,face = "bold"))
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig3\\immune.checkpoint.score.response.mana.cells.pdf",height = 4,width =4)



# correlation analysis with immune checkpoint score
ser<-FindVariableFeatures(ser, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ser), 3000)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ser)



ser.r<-subset(ser,subset=response=='R')
mprexpr.tumor <- ser.r@assays$RNA@scale.data
ser.nr<-subset(ser,subset=response=='NR'&antigen %in% c('MANA','mana_score'))
nonmprexpr.tumor <- ser.nr@assays$RNA@scale.data


corvec <- apply(mprexpr.tumor,1,cor,ser.r$immune.checkpoint.score1)
mprcorvec <- sort(corvec,decreasing = T)
mprcorvec<-as.data.frame(mprcorvec)
library(data.table)
setDT(mprcorvec, keep.rownames = TRUE)[]
mpr30<-mprcorvec %>% top_n(30) %>% filter(rn %in% VariableFeatures(ser))

corvec <- apply(nonmprexpr.tumor,1,cor,ser.nr$immune.checkpoint.score1)
nonmprcorvec <- sort(corvec,decreasing = T)
nonmprcorvec<-as.data.frame(nonmprcorvec)
setDT(nonmprcorvec, keep.rownames = TRUE)[]
nonmpr30<-nonmprcorvec %>% top_n(30) %>% filter(rn %in% VariableFeatures(ser))
plotgene<-c(mpr30$rn,nonmpr30$rn)
nonmpr30<-nonmprcorvec %>% filter(rn %in% plotgene)
mpr30<-mprcorvec %>% filter(rn %in% plotgene)

nonmpr30$mpr<-"Non-MPRs"
colnames(nonmpr30)[2]<-"Cof"
mpr30$mpr<-"MPRs"
colnames(mpr30)[2]<-"Cof"


# correlation of top 30 genes by responders vs non-responders
cof<-rbind(nonmpr30,mpr30)

cof.wide<-dcast(cof, rn~mpr,value.var="Cof" ) 
cof.wide$dif<-cof.wide$`Non-MPRs`-cof.wide$MPRs


cof.wide$dif.label[cof.wide$dif>0]<-'Higher correlation in non-MPR'
cof.wide$dif.label[cof.wide$dif<=0]<-'Higher correlation in MPR' 
cof.wide<-cof.wide %>% arrange(dif)

# plot difference  
ggbarplot(cof.wide %>% filter(!rn %in% ex), x = "rn", y = "dif",
          fill = "dif.label",           # change fill color by mpg_level
          color = "dif.label",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = T,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "",
          xlab = "",
          legend.title = "",
          #facet.by ="dif.label",
          rotate = F
          
)
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig3\\cor.dif.immune.checkpoints.pdf",width = 10,height = 5)

