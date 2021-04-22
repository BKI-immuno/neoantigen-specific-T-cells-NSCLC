## Fig4e (MANA prop differs w2,w4)
## Null hypothesis: ratio = MANA_prop/ctprop is the same across time w2 and w4
library(dplyr)
## count_data: data frame of number of cells listed in 4 columns: each row represents a cell type; column names: w2_mana, w2_all, w4_mana, w4_all
construct_stat <- function(count_data){
    
    prop_data = apply(count_data,2,function(x) x/sum(x)) %>% data.frame
    # calculate ratio
    ratio_data = prop_data %>% mutate(w2_ratio = w2_mana/w2_cell,
                                      w4_ratio = w4_mana/w4_cell) %>% select(w2_ratio,w4_ratio)
    # stat
    stat = sum((ratio_data$w2_ratio-ratio_data$w4_ratio)^2)
    return(stat)
}

## process data: pbmc.stat.rds
raw = read.csv('~/Downloads/Fig4e.csv')

## observed
obs_count = raw %>% select(-celltype)
obs_prop = apply(obs_count,2,function(x) x/sum(x)) %>% data.frame
obs_stat = construct_stat(obs_count)

## null distribution
count_all = obs_count %>% mutate(mana = w2_mana + w4_mana,
                                  cell = w2_cell + w4_cell) %>% select(mana,cell)
prop_all = apply(count_all,2,function(x) x/sum(x)) %>% data.frame
# calculate ratio
common_ratio = prop_all$mana/prop_all$cell
# calculate under null hypothesis: mana prop = ratio * ctprop
num_mana_w2 = sum(obs_count$w2_mana)
num_mana_w4 = sum(obs_count$w4_mana)

null_w2_mana = common_ratio * obs_prop$w2_cell; null_w2_mana = null_w2_mana/sum(null_w2_mana)
null_w4_mana = common_ratio * obs_prop$w4_cell; null_w4_mana = null_w4_mana/sum(null_w4_mana)


get_null <- function(B = 1000){
    
    null_stat = c()
    for(b in 1:B){
        
        set.seed(b)
        tmp_w2 = rmultinom(1,num_mana_w2,null_w2_mana)
        set.seed(b)
        tmp_w4 = rmultinom(1,num_mana_w4,null_w4_mana)
        
        tmp_stat = construct_stat(data.frame(w2_mana = tmp_w2, w2_cell = obs_count$w2_cell,
                                             w4_mana = tmp_w4, w4_cell = obs_count$w4_cell))
        null_stat = c(null_stat,tmp_stat)
        
    }
    null_stat
}

## p-value
B=10^4
null_vector = get_null(B)
# two-sided
two.sided = sum(abs(null_vector) >= abs(obs_stat))/B # p-value < 1e-4
two.sided
