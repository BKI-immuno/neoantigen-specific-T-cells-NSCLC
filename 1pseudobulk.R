library(Matrix)
library(parallel)
library(dplyr)
setwd('~/scratch/TCR/data/1global_treated/') # set working directory

## ====== Extract info ===== ##
ser = readRDS('/scratch/users/zji4@jhu.edu/tmp/ser.integrate.rds')
count = ser@assays$RNA@counts
saveRDS(count,'count.rds')
meta = readRDS('/scratch/users/zji4@jhu.edu/tmp/meta.rds')
meta$group <- paste0(meta$batch,meta$center)
meta$comb.lane <- paste0(sub('-[0-9]*_.*','',rownames(meta)),'_',meta$group)
aa <- meta %>% group_by(sample.n) %>% tally
meta.new = left_join(meta,aa)
saveRDS(meta.new,'meta.rds')
# 
cluster = meta$CellType
names(cluster) = meta$barcode
saveRDS(cluster,'cluster.rds')


## ====== GLOBAL PART ===== ##
output.dir = './resubmission/'
##############
## 15 refined clusters
##############
s <- readRDS('count.rds')
clu <- readRDS('cluster.rds')
uc <- unique(clu) %>% sort
meta <- readRDS('meta.rds') %>% select(barcode,orig.ident,patient_id,tissue,comb.lane,sample.n) %>% mutate(patient.tissue = paste0(patient_id,':',tissue))

## pb - count
m3 <- mclapply(uc,function(sc) {
    print(sc)
    tmp <- s[,names(clu)[which(clu==sc)]]
    new <- inner_join(data.frame(barcode = colnames(tmp)),meta)
    p <- new$patient.tissue
    sapply(unique(p),function(sp) rowSums(tmp[,p==sp,drop=F]))
    
},mc.cores = min(detectCores(),length(uc)))
names(m3) <- uc
# saveRDS(m3, 'pb_count_patienttissue.rds')

## pb - norm
d3 <- lapply(m3,function(i) {
    rc <- colSums(i)/1e6
    log2(t(t(i)/rc + 1))
})
saveRDS(d3, paste0(output.dir,'pb_norm_patienttissue.rds'))

## ====== CD8+ UMAP ===== ##
##############
## MANA combined - combine 4 MANA enriched cell subsets (tumor)
###############
## Identify MANA enriched clusters
tt = mana.pos %>% filter(tissue=='tumor') %>% group_by(CellType,tissue) %>% summarize(n.cluster = n_distinct(barcode))
tt$ctprop = tt$n.cluster/sum(tt$n.cluster)
tt.mana = mana.pos %>% filter(tissue=='tumor' & antigen=='MANA') %>% group_by(CellType,tissue,antigen)%>% summarise(obs = n_distinct(barcode))
tt.comb = left_join(tt,tt.mana) %>% mutate(num.mana = sum(tt.mana$obs,na.rm=T),expect = ctprop*num.mana,foldchange = obs/expect) %>% arrange(-foldchange)
saveRDS(tt.comb,paste0(output.dir,'manaprop.cutoff.rds'))

##
cutoff = readRDS(paste0(output.dir,'manaprop.cutoff.rds')) %>% filter(foldchange>2)
s <- readRDS('count.rds')
mana.pos = readRDS('meta.cd8.rds') %>% filter(CellType %in% cutoff$CellType) %>% filter(tissue=='tumor')
mana.pos$patient.tissue = paste0(mana.pos$patient_id,':',mana.pos$tissue)

tmp <- s[,colnames(s) %in% mana.pos$barcode]
new <- inner_join(data.frame(barcode = colnames(tmp)),mana.pos)
p <- new$patient.tissue
m6 <- sapply(unique(p),function(sp) rowSums(tmp[,p==sp,drop=F]))
# saveRDS(m6, 'finalcomb4clu.tumor.bulk.pb_count_patienttissue.rds')

## normalize
rc2 <- colSums(m6)/1e6
d6 <- log2(t(t(m6)/rc2 + 1))
saveRDS(d6, paste0(output.dir,'finalcomb4clu.tumor.bulk.pb_norm_patienttissue.rds'))

##############
## CD8+ individual cell subsets (tumor)
##############
## pb - count
s <- readRDS('count.rds')
mana.pos =  readRDS('meta.cd8.rds') %>% filter(tissue=='tumor')
mana.pos$patient.tissue = paste0(mana.pos$patient_id,':',mana.pos$tissue)

s <- s[,colnames(s) %in% mana.pos$barcode]
clu <- mana.pos$CellType
names(clu) <- mana.pos$barcode
uc <- unique(mana.pos$CellType) %>% sort
m5 <- mclapply(uc,function(sc) {
    print(sc)
    tmp <- s[,names(clu)[which(clu==sc)]]
    new <- inner_join(data.frame(barcode = colnames(tmp)),mana.pos)
    p <- new$patient.tissue
    sapply(unique(p),function(sp) rowSums(tmp[,p==sp,drop=F]))
    
},mc.cores = min(detectCores(),length(uc)))
names(m5) <- uc
# saveRDS(m5, 'finalbyclu.tumor.pb_count_patienttissue.rds')

## pb - norm
d5 <- lapply(m5,function(i) {
    rc <- colSums(i)/1e6
    log2(t(t(i)/rc + 1))
})
saveRDS(d5, paste0(output.dir,'finalbyclu.tumor.pb_norm_patienttissue.rds'))

