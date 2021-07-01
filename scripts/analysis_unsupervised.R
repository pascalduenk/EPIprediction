setwd("/lustre/nobackup/WUR/ABGC/duenk002/EPI_prediction")

library(dplyr)
library(data.table)
library(tidyr)
library(caret)
library(ggplot2)
library(PRROC)

set.seed(110789)

theme_set(theme_classic(base_size = 10) + 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                plot.title = element_text(hjust = 0.5),
                plot.margin = unit(c(1, 1, 1, 1), "lines")))

Matt_Coef <- function (conf_matrix)
{
  TP <- conf_matrix$table[1,1]
  TN <- conf_matrix$table[2,2]
  FP <- conf_matrix$table[1,2]
  FN <- conf_matrix$table[2,1]
  
  mcc_num <- (TP*TN - FP*FN)
  mcc_den <- as.double((TP+FP))*as.double((TP+FN))*as.double((TN+FP))*as.double((TN+FN))
  
  mcc_final <- mcc_num/sqrt(mcc_den)
  return(mcc_final)
}

# Options -----------------------------------------------------------------
options(stringsAsFactors = F)

# Read data -----------------------------------------------------------------
dat <- fread("dataTrainTest.csv", header=T, sep="\t")

# select predictors 
vec_predictors <- paste0(c("H3K27ac.enhancer", "H3K27ac.promoter", "H3K27ac.window",
                    "H3K27me3.enhancer", "H3K27me3.promoter", "H3K27me3.window",
                    "H3K4me1.enhancer", "H3K4me1.promoter", "H3K4me1.window",
                    "Methylation.enhancer", "Methylation.promoter", "Methylation.window",
                    "atac.enhancer","atac.promoter","atac.window","gene.expression"), 
                    ".log")
# Preprocess -----------------------------------------------------------------
dat_scaled <- dat %>%
    select(all_of(vec_predictors)) %>%
    mutate_all(scale)

# optimize k 
ks <- seq(1,200,10)   # look at what happens with 1:10
all_ttws <- vector("numeric", length=length(ks))
for(i in 1:length(ks)){
    k_opt <- kmeans(dat_scaled, centers = ks[i], nstart=25, iter.max = 500)
    all_ttws[i] <- k_opt$tot.withinss
}

ttws <- data.frame(ks, all_ttws)

ggplot(ttws, aes(x=ks, y=all_ttws))+
    geom_line()+
    geom_point()
ggsave("optimal_k.pdf", width = 75, height = 75, scale = 1.5, units = "mm")

# k-means
chosen_k <- 30
km <- kmeans(dat_scaled, centers = chosen_k, nstart=100, iter.max = 500)

pred <- km$cluster

mc_percluster <- vector("numeric", chosen_k)
F1_percluster <- vector("numeric", chosen_k)

for(i in 1:chosen_k){
    cl_pred <- km$cluster == i
    cl_cf <- confusionMatrix(factor(cl_pred), factor(dat$label))
    mc_percluster[i] <- Matt_Coef(cl_cf)
    F1_percluster[i] <- cl_cf$byClass[7]
}

res_kmeans <- data.frame(mc_percluster, F1_percluster)
write.table(res_kmeans, "k-means_result.txt", row.names=F, quote=F)

# hierarchical clustering 
dm <- dist(dat_scaled)
hc <- hclust(dm)

allfreqs <- matrix(0, 20, 20)

for(i in 1:20){
    ct <- cutree(hc, i)

    freqtable <- vector("numeric", i)

    for(j in 1:i){
        clustj <- dat %>% slice(which(ct == j))
        ntrue <- nrow(clustj[clustj$label == 1,])
        nfalse <- nrow(clustj[clustj$label == 0,])
        freqtable[j] <- ntrue / (ntrue + nfalse)
    }
    allfreqs[i,1:i] <- freqtable
}

write.table(allfreqs, "cluster_frequencies.txt", quote=F, col.names=F, row.names=F)

pdf("dendogram.pdf", width = 21, height = 14)
plot(hc, labels=FALSE)
dev.off()

table(is.na(dat_scaled))