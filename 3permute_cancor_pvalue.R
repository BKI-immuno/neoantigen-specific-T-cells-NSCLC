library(tidyverse)

## dat: input data (after running 2-xpca_xx.R)
## var: response variable (cancor(gene expr, var))
## B: permutation times
## permute: indicator whether conduct permutation test or not

permute_axis <- function(dat,
                         var = 'resi_tumor',
                         B=1000,
                         permute=T){
    
    ## canonical correlation
    cancor = cancor(dat %>% select(PC1,PC2),dat[,var])
    obs.cor = cancor$cor
    
    if(permute){
        
        new.idx = t(sapply(1:B,function(b){
            set.seed(b)
            vec = sample(1:nrow(dat),replace = F)
            if(identical(vec,1:nrow(dat))) vec = NULL
            
            return(vec)
            
        })) %>% unique
        print(paste0('#unique permuted rows: ',nrow(new.idx)))
        
        new.cor = sapply(1:nrow(new.idx),function(i){
                              
                              newdat = dat
                              newdat[,var] = dat[new.idx[i,],var]
                              cancor(newdat %>% select(PC1,PC2),newdat[,var])$cor
                              
                          })
        
        p.value = sum(new.cor>=obs.cor)/nrow(new.idx)
        print(paste0('=== Observed cancor: ',obs.cor))
        print(paste0('=== P-value: ',p.value))
        
        list(cancor = cancor, obs.cor = obs.cor,permute.cor = new.cor,p.value = p.value)
  
    }
}

##############
## Example
##############
setwd('~/scratch/TCR/')

## (1) Global UMAP: limit to tumor tissue
output.dir = './result/1global_treated/pca/pb_patient/resubmission/tumor/'
file = 'dat_hvg_m2.rds'
dat = readRDS(paste0(output.dir,file)) # result file by running 2-xpca_xx.R (dat_hvg_m2.rds)
dat$response = ifelse(dat$response=='NR',0,1)
dat$tissue = ifelse(dat$tissue=='normal',0,1)
res.tumor = permute_axis(dat, var = 'response', B = 10000, permute = T)

## (2) Global UMAP: limit to normal tissue
output.dir = './result/1global_treated/pca/pb_patient/resubmission/normal/'
file = 'dat_hvg_m2.rds'
dat = readRDS(paste0(output.dir,file)) # result file by running 2-xpca_xx.R (dat_hvg_m2.rds)
dat$response = ifelse(dat$response=='NR',0,1)
dat$tissue = ifelse(dat$tissue=='normal',0,1)
res.norm = permute_axis(dat, var = 'response', B = 10000, permute = T)

# (3) Global UMAP: limit to tumor and normal
output.dir = './result/1global_treated/pca/pb_patient/resubmission/overall/'
file = 'dat_hvg_m2.rds'
dat = readRDS(paste0(output.dir,file)) # result file by running 2-xpca_xx.R (dat_hvg_m2.rds)
dat$response = ifelse(dat$response=='NR',0,1)
dat$tissue = ifelse(dat$tissue=='normal',0,1)
res.tis = permute_axis(dat, var = 'tissue', B = 10000, permute = T)

# (4) CD8 UMAP: limit to tumor tissue
output.dir = './result/1global_treated/pca/pb_patient/resubmission/CD8/tumor/'
file = 'dat_hvg_m2.rds'
dat = readRDS(paste0(output.dir,file)) # result file by running 2-xpca_xx.R (dat_hvg_m2.rds)
dat$response = ifelse(dat$response=='NR',0,1)
dat$tissue = ifelse(dat$tissue=='normal',0,1)
res.cd8 = permute_axis(dat, var = 'response', B = 10000, permute = T)

# (5) Combined MANA cluster
output.dir = './result/1global_treated/pca/pb_patient/resubmission/CD8/final_comb4clu_tumor/'
file = 'dat_hvg_m2.rds'
dat = readRDS(paste0(output.dir,file)) # result file by running 2-xpca_xx.R (dat_hvg_m2.rds)
dat$response = ifelse(dat$response=='NR',0,1)
dat$tissue = ifelse(dat$tissue=='normal',0,1)
res.mana = permute_axis(dat, var = 'response', B = 10000, permute = T)
