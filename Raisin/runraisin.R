library(parallel)
library(Matrix)
source('/home-4/zji4@jhu.edu/scratch/raisin/software/raisin/raisin.R')
meta <- readRDS('/home-4/zji4@jhu.edu/scratch/TCR/integrate/meta.rds')
meta <- meta[sub('-.*','',sub('.*:','',meta$sample.n)) %in% c('tumor','normal'),]
meta$sample.n <- paste0(meta$sample.n,meta$center)
meta$CellType <- sub('\\)','',sub('\\(','-',gsub(' ','',sub('/','',meta$CellType))))
restumor <- unique(meta[,c('patient_id','response')])
clu <- as.character(meta$CellType)
names(clu) <- meta$barcode
combclu <- clu
combclu[grep('CD4',combclu)] <- 'CD4'
combclu[grep('CD8',combclu)] <- 'CD8'
print(length(clu))

s <- readRDS('/home-4/zji4@jhu.edu/scratch/TCR/saver/full/combine/tumornormalsaver.rds')
s <- s[,names(clu)]
s <- s[rowSums(s) > 0,]

patient <- sub(':.*','',meta$barcode)
patlist <- unique(data.frame(meta$response,patient))
tumornormal <- meta$tissue

for (sc in c(unique(clu),'CD4','CD8')) {
  print(sc)
  for (i in unique(tumornormal)) {
    if (sc %in% c('CD4','CD8')) {
      cellid <- names(which(combclu == sc & tumornormal==i))  
    } else {
      cellid <- names(which(clu == sc & tumornormal==i))  
    }
    sample <- paste0(sub('-[0-9]*_.*','',cellid))
    samplename <- unique(sample)
    
    pb <- rowsum(t(s[,cellid]),sample)
    pb <- pb/as.vector(table(sample)[rownames(pb)])
    selgn <- names(which(colSums(pb > 0.1) > 0))
    
    tg <- restumor[match(sub(':.*','',samplename),restumor[,1]),2]
    design <- data.frame(sample=samplename,feature=tg)
    print(design)
    
    if (length(unique(design[,'feature'])) > 1) {
      res <- RAISINtest(RAISINfit(s[selgn,cellid],sample,testtype='unpaired',design=design,filtergene=F,ncores=15))
      
      write.csv(res,file=gsub(' ','',paste0('/home-4/zji4@jhu.edu/scratch/TCR/diff/res/full_RNR/',sc,':',i,'.csv')),quote=F)
    }
  }  
}


