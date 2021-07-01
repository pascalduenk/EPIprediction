# combine results into tables and figures
setwd("/lustre/nobackup/WUR/ABGC/duenk002/EPI_prediction")

library(dplyr)
library(data.table)
library(tidyr)

# scenarios
scens <- c("chr1_log", "rand_log")

# read model performance
mp <- data.frame(metric = NA, model = NA, value = NA, scen = NA)

for(s in scens){
    res <- read.csv(paste0("results_alldata/", s, "/model_performance.txt"),                            
                    col.names = c('metric', 'model', 'value'), sep = ",")
                 
    res$scen <- s
    mp <- rbind(mp, res)
}

mp <- mp %>% filter(metric %in% c("Specificity", "Precision", "Recall", "F1", "Matthews Correlation"))

mp_table <- mp %>% 
    pivot_wider(names_from = model, values_from = value) %>% 
    mutate_if(is.numeric, ~ round(., 3)) %>%
    as.data.frame()

write.table(mp_table, "table_model_performance.csv", col.names=T, row.names=F, quote=F, sep=",")