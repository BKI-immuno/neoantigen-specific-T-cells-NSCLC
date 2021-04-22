library(tidyverse)
library(sva)
library(ggplot2)

setwd('~/scratch/TCR/')

##############################
## Gene expression PCA (concatenate features): global UMAP or global CD8 UMAP
##############################

## meta: sample meta
## s: gene expression list (each element is a pseudobulk at sample-level for a specific celltype)
## output.dir: output directory

run_pca_concat_features <- function(meta,s,output.dir){
    method = "m2"
    if(!dir.exists(output.dir)) dir.create(output.dir)
        
d <- sapply(1:length(s),function(i) {
    
    d <- s[[i]][!grepl('MT-',rownames(s[[i]])),]
    d <- d[rowMeans(d > 0.01) > 0.1,]
    #print(identical(meta$orig.ident,colnames(d)))
    new = inner_join(data.frame(patient.tissue = colnames(d)),meta)
    batch = new$center
    residual.tumor = new$resi_tumor
    tissue = new$tissue
    
    # further exclude genes have 0 var within each batch
    genes.sd = lapply(unique(batch),function(b){
        
        d.sub = d[,which(batch==b),drop=F]
        if(ncol(d.sub)>1){
            sd = rowSds(d.sub)
            g = rownames(d.sub)[sd > 0]
        }else{
            g = rownames(d.sub)#[d.sub>0]
        }
        g
    })
    genes.tbl = table(genes.sd %>% unlist)
    genes.include = names(genes.tbl)[genes.tbl==length(unique(batch))] # not too many
    
    ## subset d
    d = d[genes.include,]
    
    # combat
    d <- ComBat(d,batch,mod = model.matrix(~residual.tumor))
    
    # select highly variable genes
    cm <- rowMeans(d)
    csd <- sqrt((rowMeans(d*d) - cm^2) / (ncol(d) - 1) * ncol(d))
    mod <- loess(csd~cm)
    rid <- which(resid(mod) > 0)
    d <- d[rid,]
    rownames(d) <- paste0(rownames(d),':',names(s)[i])
    
    # standardize
    (d-rowMeans(d))/apply(d,1,sd)
    
})

## concatenate features
allsamp = which(sapply(d,ncol) %>% unlist == max(sapply(d,ncol) %>% unlist))
d <- do.call(rbind,d[allsamp])
pseudo_bulk = d
# saveRDS(d, paste0(output.dir,"combat_pb_hvg_",method,".rds"))

## PCA
pr = prcomp(t(pseudo_bulk),scale. = T)$x
pseudo_bulk_pca <- (pr[,1:min(20,ncol(d))])
# saveRDS(pseudo_bulk_pca, paste0(output.dir,"pca_hvg_",method,".rds"))

pca.hvg = data.frame(pseudo_bulk_pca[,1:2],patient.tissue = rownames(pseudo_bulk_pca))
dat.hvg = inner_join(pca.hvg,meta) # %>% mutate(log10.num.cells = log10(num.cells))
saveRDS(dat.hvg,paste0(output.dir,'dat_hvg_',method,'.rds'))

return(dat.hvg)
}

plot_pca <- function(dt,column){
    min = min(dt$PC1,dt$PC2)
    max = max(dt$PC1,dt$PC2)
    if(column == 'tissue'){
        
        (ggplot(dt,aes_string(x='PC1',y='PC2',color = column,label = 'patient.tissue')) +
             geom_point() +
             theme_classic() +
             scale_color_manual(breaks=c('normal','tumor'),
                                values=c('#F8766D','#00BFC4')) +
             xlim(min,max) + ylim(min,max)) %>% print
        
    }else{
        (ggplot(dt,aes_string(x='PC1',y='PC2',color = column)) +
             geom_point() +
             theme_classic() +
             xlim(min,max) + ylim(min,max)) %>% print
    }
}

##############
## Example
##############
## (1) Global UMAP: limit to normal tissue
output.dir = './result/1global_treated/pca/pb_patient/resubmission/normal/'
pseudobulk.dat = 'pb_norm_patienttissue.rds'
m = readRDS('./data/1global_treated/meta.rds') %>% select(barcode,orig.ident,batch,patient_id,center,cohort,tissue,patient.tissue,resi_tumor,sample.n,group,response) %>% filter(tissue %in% c('normal'))
clu = readRDS("./data/1global_treated/cluster.rds")
clu = clu[names(clu) %in% m$barcode]
meta = m %>% select(patient.tissue,center,tissue,resi_tumor,patient_id,response) %>% unique

s = readRDS(paste0("./data/1global_treated/resubmission/",pseudobulk.dat))
s = lapply(s,function(dt) dt[,which(colnames(dt) %in% meta$patient.tissue)])
sampsize = sapply(s,ncol)
print(sampsize)  # 12 normal

dat_pca = run_pca_concat_features(meta,s,output.dir)
# plot
pdf(paste0(output.dir,'pca_plot.pdf'))
plot_pca(dat_pca,'response')
dev.off()


## (2) Global UMAP: limit to tumor tissue
output.dir = './result/1global_treated/pca/pb_patient/resubmission/tumor/'
pseudobulk.dat = 'pb_norm_patienttissue.rds'
m = readRDS('./data/1global_treated/meta.rds') %>% select(barcode,orig.ident,batch,patient_id,center,cohort,tissue,patient.tissue,resi_tumor,sample.n,group,response) %>% filter(tissue %in% c('tumor'))
clu = readRDS("./data/1global_treated/cluster.rds")
clu = clu[names(clu) %in% m$barcode]
meta = m %>% select(patient.tissue,center,tissue,resi_tumor,patient_id,response) %>% unique

s = readRDS(paste0("./data/1global_treated/resubmission/",pseudobulk.dat))
s = lapply(s,function(dt) dt[,which(colnames(dt) %in% meta$patient.tissue)])
sampsize = sapply(s,ncol)
print(sampsize)  # 15 tumor

dat_pca = run_pca_concat_features(meta,s,output.dir)
# plot
pdf(paste0(output.dir,'pca_plot.pdf'))
plot_pca(dat_pca,'response')
dev.off()

## (3) Global UMAP: tumor and normal
output.dir = './result/1global_treated/pca/pb_patient/resubmission/overall/'
pseudobulk.dat = 'pb_norm_patienttissue.rds'
m = readRDS('./data/1global_treated/meta.rds') %>% select(barcode,orig.ident,batch,patient_id,center,cohort,tissue,patient.tissue,resi_tumor,sample.n,group,response) %>% filter(tissue %in% c('tumor','normal'))
clu = readRDS("./data/1global_treated/cluster.rds")
clu = clu[names(clu) %in% m$barcode]
meta = m %>% select(patient.tissue,center,tissue,resi_tumor,patient_id,response) %>% unique

s = readRDS(paste0("./data/1global_treated/resubmission/",pseudobulk.dat))
s = lapply(s,function(dt) dt[,which(colnames(dt) %in% meta$patient.tissue)])
sampsize = sapply(s,ncol)
print(sampsize)  # 27

dat_pca = run_pca_concat_features(meta,s,output.dir)
# plot
pdf(paste0(output.dir,'pca_plot.pdf'))
plot_pca(dat_pca,'tissue')
dev.off()

## (4) CD8 UMAP: limit to tumor tissue
output.dir = './result/1global_treated/pca/pb_patient/resubmission/CD8/tumor/'
pseudobulk.dat = 'finalbyclu.tumor.pb_norm_patienttissue.rds'
m = readRDS('./data/1global_treated/meta.cd8.rds') %>% select(barcode,orig.ident,batch,patient_id,center,cohort,tissue,patient.tissue,resi_tumor,sample.n,response,CellType) %>% filter(tissue %in% c('tumor'))
meta = m %>% select(patient.tissue,center,tissue,resi_tumor,patient_id,response) %>% unique

s = readRDS(paste0("./data/1global_treated/resubmission/",pseudobulk.dat))
s = lapply(s,function(dt) dt[,which(colnames(dt) %in% meta$patient.tissue)])
sampsize = sapply(s,ncol)
print(sampsize)  # 15 tumor

dat_pca = run_pca_concat_features(meta,s,output.dir)
# plot
pdf(paste0(output.dir,'pca_plot.pdf'))
plot_pca(dat_pca,'response')
dev.off()
