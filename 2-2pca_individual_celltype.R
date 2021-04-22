library(tidyverse)
library(ggplot2)
library(sva)

setwd('~/scratch/TCR/')

##############################
## Gene expression PCA
##############################
## meta: sample meta
## s: gene expression matrix of 1 cell type
## output.dir: output directory
run_pca <- function(meta,s,output.dir){
    
    method = "m2"
    if(!dir.exists(output.dir)) dir.create(output.dir)

d <- s[!grepl('MT-',rownames(s)),]
d <- d[rowMeans(d > 0.01) > 0.1,]
#print(identical(meta$orig.ident,colnames(d)))
new = inner_join(data.frame(patient.tissue = colnames(d)),meta)
batch = new$center
residual.tumor = new$resi_tumor
tissue = new$tissue

## further exclude genes have 0 var within each batch
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
genes.include = names(genes.tbl)[genes.tbl==length(unique(batch))]
## subset d
d = d[genes.include,]
## combat
d <- ComBat(d,batch,mod = model.matrix(~residual.tumor))
## select highly variable genes
cm <- rowMeans(d)
csd <- sqrt((rowMeans(d*d) - cm^2) / (ncol(d) - 1) * ncol(d))
mod <- loess(csd~cm)
rid <- which(resid(mod) > 0)
d <- d[rid,]
## standardize
d = (d-rowMeans(d))/apply(d,1,sd)
pseudo_bulk = d
## PCA
pr = prcomp(t(pseudo_bulk),scale. = T)$x
pseudo_bulk_pca <- (pr[,1:min(20,ncol(d))])
## summarize
pca.hvg = data.frame(pseudo_bulk_pca[,1:2],patient.tissue = rownames(pseudo_bulk_pca))
print(dim(pca.hvg))
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

###############
## Example
###############
## (1) PCA for combined MANA enriched cluster
meta.dat = './data/1global_treated/meta.cd8.rds'
pseudobulk.dat = 'finalcomb4clu.tumor.bulk.pb_norm_patienttissue.rds'
output.dir = './result/1global_treated/pca/pb_patient/resubmission/CD8/final_comb4clu_tumor/'

m = readRDS(meta.dat) %>%
    select(barcode,orig.ident,batch,patient_id,center,cohort,tissue,response,resi_tumor,sample.n) %>%
    mutate(patient.tissue = paste0(patient_id,':',tissue)) %>% filter(tissue=='tumor')
aa = m %>% group_by(patient.tissue) %>% summarise(num.cells = n_distinct(barcode))
m = left_join(m,aa)
meta = m %>% select(patient.tissue,center,tissue,resi_tumor,patient_id,response,num.cells) %>% unique

s = readRDS(paste0("./data/1global_treated/",pseudobulk.dat))
s = s[,which(colnames(s) %in% meta$patient.tissue)]
print(ncol(s))  # 15 CD8+ tumor samples

dat_pca = run_pca(meta,s,output.dir)
# plot
pdf(paste0(output.dir,'pca_plot.pdf'))
plot_pca(dat_pca,'response')
dev.off()

## (2) PCA for other non MANA enriched clusters (evaluate separately)
unq.ct = c("Effector(I)", "TRM(I)", "Effector(II)", "MAIT", "TRM(VI)", "CD4CD8(I)", "Effector(III)", "TRM(III)","Stem-like memory", "CD4CD8(II)")

cancor.ls = c()
for(ctname in unq.ct){
    print(ctname)
    
    output.dir = paste0('./result/1global_treated/pca/pb_patient/resubmission/CD8/',gsub('\\/','or',ctname),'/')
    pseudobulk.dat = 'finalbyclu.tumor.pb_norm_patienttissue.rds'

    ## Limit to tumor samples and a specific celltype
    m = readRDS('./data/1global_treated/meta.cd8.rds') %>% select(barcode,orig.ident,batch,patient_id,center,cohort,tissue,patient.tissue,resi_tumor,sample.n,response,CellType) %>%
    mutate(CellType = as.character(CellType)) %>%
    filter(tissue %in% c('tumor')) %>% filter(CellType==ctname)
    clu = m$CellType
    names(clu) = m$barcode
    meta = m %>% select(patient.tissue,center,tissue,resi_tumor,patient_id,response) %>% unique

    s = readRDS(paste0("./data/1global_treated/",pseudobulk.dat))
    s = s[[which(names(s) %in% unique(clu))]] # extract 1 cluster
    print(ncol(s))  # 15 tumor

    dat_pca = run_pca(meta,s,output.dir) %>% mutate(response = ifelse(response=='NR',0,1))

    cancor.ls = c(cancor.ls,cancor(dat_pca[,1:2],dat_pca$response)$cor)
}

cancor.dat = data.frame(CellType = unq.ct,cancor = cancor.ls)
