options(stringsAsFactors = F)

library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)

theme_set(theme_classic(base_size = 10) + 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                plot.title = element_text(hjust = 0.5),
                plot.margin = unit(c(1, 1, 1, 1), "lines")))
              
# Read data --------------------------------------------
mods <- c("rf","gbm","svm")

for(m in mods){
    vi <- read.csv(paste0("var_importance_",m,".txt"), header=F, skip=1, sep=" ")
    colnames(vi) <- c("feature", "importance")
    vi$imp_round <- round(vi$importance, 0)

    ggplot(vi, aes(x=reorder(feature, importance), y=importance)) +
        geom_col(fill = "lightgrey")+
        coord_flip()+
        labs(x="Feature",y="Importance")+
        geom_text(aes(label = imp_round), hjust = -0.5)+
        ylim(0, 115)

        ggsave(paste0("vi_",m,".pdf"), width = 85, height = 75, scale = 1.5, units = "mm")
}

theme_set(theme_classic(base_size = 10) + theme(plot.margin = unit(c(1,1,1,1), "lines"))