setwd("/lustre/nobackup/WUR/ABGC/duenk002/EPI_prediction")

library(dplyr)
library(data.table)
library(tidyr)
library(caret)
library(randomForest)
library(kernlab)
library(ggplot2)
library(GGally)
library(gridExtra)
library(PRROC)
library(gbm)
library(parallel)
library(doParallel)

set.seed(110789)

# Options -----------------------------------------------------------------
options(stringsAsFactors = F)

script_args = commandArgs(trailingOnly=TRUE)

if(script_args[1] == "all"){
  selchrom <- paste0("chr", seq(1:22))
} else {
  selchrom <- paste0("chr", seq(1:6))
}

fold <- ifelse(script_args[1] == "all", "results_alldata", "results")
transf <- script_args[2]
cvmethod <- script_args[3]
target_var <- script_args[4]

# Functions ---------------------------------------------------------------

# Function to compute AU-PRC for model training
auprcSummary <- function(data, lev = NULL, model = NULL){
  index_class2 <- data$obs == "true"
  index_class1 <- data$obs == "false"
  
  the_curve <- pr.curve(data$true[index_class2], data$false[index_class1], curve = FALSE)
  out <- the_curve$auc.integral
  names(out) <- "AUPRC"
  
  out
}

# Matthews correlation
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

# density plot of each variable 
histplot <- function(df, vbl){
  ggp <- ggplot(df, aes_string(x = vbl))+
    geom_histogram(color="darkblue", fill="lightblue")+
    ggtitle(vbl)+
    labs(x="")
  return(ggp)
}

theme_set(theme_classic(base_size = 10) + 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                plot.title = element_text(hjust = 0.5)))
              

# choose number of cores to use for parallel computing
cl <- makeCluster(6)
registerDoParallel(cl)
clusterEvalQ(cl, library(PRROC))   # load package needed in each core

# Read data ---------------------------------------------------------------
dat <- fread("data_complete.csv", header=T) %>%
  filter(promoter_chrom %in% selchrom)

dat <- dat %>%
  mutate(!!target_var := recode_factor(factor(dat %>% pull(target_var)),
         `1` = "true", `0` = "false"))

# change column names
cnames <- colnames(dat)
cnames2 <- sapply(cnames, function(x){gsub("\\(","",x)})
cnames2 <- sapply(cnames2, function(x){gsub("\\)","",x)})
cnames2 <- sapply(cnames2, function(x){gsub("_",".",x)})
cnames2 <- sapply(cnames2, function(x){gsub(" ",".",x)})

colnames(dat) <- cnames2

# select predictors 
vec_predictors <- c("H3K27ac.enhancer", "H3K27ac.promoter", "H3K27ac.window",
                    "H3K27me3.enhancer", "H3K27me3.promoter", "H3K27me3.window",
                    "H3K4me1.enhancer", "H3K4me1.promoter", "H3K4me1.window",
                    "Methylation.enhancer", "Methylation.promoter", "Methylation.window",
                    "atac.enhancer","atac.promoter","atac.window","gene.expression")

# Remove gene expression values above a certain level 
dat <- dat %>%
  filter(gene.expression < 1000) 

# Summary statistics and plots --------------------------------------------

sum_predictors <- dat %>%
  select(all_of(vec_predictors)) %>%
  summarise_all(list(mean = function(x){mean(x, na.rm=T)}, 
                     sd = function(x){sd(x, na.rm=T)},
                     min = function(x){min(x, na.rm=T)}, 
                     max = function(x){max(x, na.rm=T)}, 
                     n = function(x){sum(!is.na(x))},
                     n.unique = function(x){length(unique(x))},
                     nNA = function(x){sum(is.na(x))})) %>%
  gather(stat,val) %>%
  separate(stat, into = c("var", "stat"), sep = "_") %>%
  spread(stat, val) 

write.csv(sum_predictors, "summary_predictors.csv", row.names = F, quote = F)
write.table(summary(dat %>% select(all_of(target_var))), "summary_target.csv", row.names = T, quote = F)

# numbers 
sum_info <- dat %>%
  select(all_of(c('promoter.name','enhancer.name'))) %>%
  summarise_all(list(n.unique = function(x){length(unique(x))}))

n_row_per_promoter <- dat %>%
  group_by(promoter.name, .drop = FALSE) %>%
  tally()

n_row_per_promoter %>%
  histplot(., "n") +
  ggtitle("")+
  labs(x="N enhancers per promoter")+
  ggsave("hist_promoters.pdf", width = 75, height = 75, scale = 1.5, units = "mm")

n_row_per_promoter %>%
  with(table(n))

# number of true and false EPIs per promoter
n_per_promoter <- dat %>%
  group_by(promoter.name, label, .drop = FALSE) %>%
  count(.drop = FALSE) %>% 
  group_by(label,n, .drop = FALSE) %>%
  count(name = "number", .drop=FALSE) 

n_per_promoter %>%
  write.table(., "n_pairs_per_promoter.txt", col.names=T, row.names=F, quote=F)

n_per_promoter %>%
  ggplot(aes(x=n, y=number))+
    geom_bar(stat = 'identity', color="darkblue", fill="lightblue")+
    facet_grid(.~label)+
    ggsave("hist_truefalse.pdf", width = 150, height = 75, scale = 1.5, units = "mm")


# # correlation plot
# corplot <- dat %>%
#   select(all_of(vec_predictors)) %>%
#   ggpairs(.) + 
#   theme_bw()
# 
# corplot
# ggsave("corplot.png", width=75,  height=75, units = "mm", scale=5)

plots <- list()
for(i in 1:length(vec_predictors)){
  plots[[i]] <- histplot(dat, vec_predictors[i])
}

ggsave("histograms.pdf", do.call(arrangeGrob, c(plots, list(nrow=4, ncol=4))),
       width=150, height=150, scale=1.5, units="mm")


# Pre processing ----------------------------------------------------------

# add a small number to values of variables that have a minimum of 0, so that we can use log transformation
min_vals <- dat %>% 
  select(all_of(vec_predictors)) %>% 
  sapply(., function(x){min(x[x>0])})

for(i in 1:length(vec_predictors)){
  v <- vec_predictors[i]
  nv <- as.numeric(dat %>% pull(v)) + 0.5*(min_vals[i])
  dat <- dat %>% mutate("{v}" := nv)
}

# log-transform values that are skewed
dat <- dat %>%
  mutate(across(all_of(vec_predictors), .fns = list(log = ~ log(.)), .names = "{.col}.log"))

# plot log transformed variables
vec_predictors_log <- paste0(vec_predictors, ".log")

plots <- list()
for(i in 1:length(vec_predictors_log)){
  plots[[i]] <- histplot(dat, vec_predictors_log[i])
}

ggsave("histograms_log.pdf", do.call(arrangeGrob, c(plots, list(nrow=4, ncol=4))),
       width=150, height=150, scale=1.5, units="mm")

# choose transformed or untransformed features
if(transf == "log"){
  vec_predictors <- vec_predictors_log
} else if(transf == "all"){
  vec_predictors <- c(vec_predictors, vec_predictors_log)
} else {
  vec_predictors <- c(vec_predictors)
}

# Model arguments and parameter tuning 
models <- c("rf", "gbm", "gbm2", "svm")

methods <- list(rf = "rf", gbm = "gbm", gbm2 = "gbm", svm = "svmLinear")

grids <- list(rf = list(expand.grid(mtry = c(5,10,20,50,100,500,1000))),
              gbm = list(expand.grid(interaction.depth = c(1, 2, 4, 6),
                    n.trees = c(500,1000,5000), 
                    shrinkage = c(0.01,0.05,0.1),
                    n.minobsinnode = c(5, 10, 20))),
              gbm2 = NULL,
              svm = list(expand.grid(C = c(0.75, 0.9, 1, 1.1, 1.25, 1.5)))
)

model.args <- list(rf = list(n.tree = 500),
                   gbm = list(NULL),
                   gbm2 = list(NULL),
                   SVM = list(NULL)
)

# write final data file 
dat %>%
  fwrite(., "dataTrainTest.csv", quote=F, row.names=F, col.names=T, sep="\t")

# Model function ----------------------------------------------------------

run_models <- function(trainIx, rf){
  
  # create train and test sets
  explanatoryTrain <- dat %>%
    slice(trainIx) %>%
    select(setdiff(vec_predictors, target_var)) %>%
    mutate_all(scale)
    
  explanatoryTest <- dat %>%
    slice(-trainIx) %>%
    select(setdiff(vec_predictors, target_var)) %>%
    mutate_all(scale) 

  targetTrain <- dat %>%
    slice(trainIx) %>%
    pull(target_var)

  targetTest <- dat %>%
    slice(-trainIx) %>%
    pull(target_var)

  # complete train and test sets
  dataTrain <- explanatoryTrain %>%
    mutate(!!target_var := targetTrain)
  
  dataTest <- explanatoryTest %>%
    mutate(!!target_var := targetTest)
    
  # plot histograms of scaled features
  plots <- list()
  for(i in 1:length(vec_predictors)){
    plots[[i]] <- histplot(explanatoryTrain, vec_predictors[i])
  }

  ggsave(paste0(rf,"/histograms_scaled2.pdf"), 
         do.call(arrangeGrob, 
                 c(plots, list(nrow=round(sqrt(length(vec_predictors))), 
                              ncol=round(sqrt(length(vec_predictors)))))),
         width=150, height=150, scale=1.5, units="mm")

  # create index for subsampling (per chromosome)
  chrom <- dat %>% slice(trainIx) %>% pull(promoter.chrom)
  u.chrom <- unique(chrom)
  subsamp.list <- vector(mode="list", length = length(u.chrom) )
  for(c in 1:length(u.chrom)){
    subsamp.list[[c]] <- which(chrom != u.chrom[c])
  }
  
  # Models ------------------------------------------------------------------
  modelFormula <- as.formula(paste(target_var, '~ .'))
  
  ctrl <- trainControl(index = subsamp.list,  
                       classProbs = TRUE,
                       summaryFunction = auprcSummary,
                       allowParallel = TRUE)
  for(m in models){
    # get arguments to pass to train function
    all.args <- c(list(modelFormula, data = dataTrain, 
                       trControl = ctrl, metric = "AUPRC", method = methods[[m]]), 
                       tuneGrid = grids[[m]], model.args[[m]]) 
    
    # train the model
    model.train <- do.call(train, all.args)

    # predictions and evalutaions
    predTest <- model.train %>% predict(explanatoryTest)
    cm <- confusionMatrix(predTest, factor(targetTest))
    mc <- Matt_Coef(cm)

    predTrain <- model.train %>% predict(explanatoryTrain)
    cm.T <- confusionMatrix(predTrain, as.factor(targetTrain))
    mc.T <- Matt_Coef(cm.T)
  
    vi.rf <- varImp(model.train)

    # write results
    sink(paste0(rf,"/summary_",m,".txt"))
    print(model.train)
    sink()

    write.table(cm$table,paste0(rf,"/cm_",m,".txt"), quote=F, row.names = F, col.names = F)
    write.table(cm.T$table,paste0(rf,"/cm_",m,"_Train.txt"), quote=F, row.names = F, col.names = F)

    perfmet <- as.data.frame(cbind(m,c(cm$byClass, mc)))
    rownames(perfmet)[nrow(perfmet)] <- "Matthews Correlation"
    write.table(perfmet, paste0(rf,"/model_performance.txt"), quote=F, row.names = T, col.names = F, append=T, sep=",")

    write.table(vi.rf['importance'], paste0(rf, "/var_importance_",m,".txt"), quote=F, row.names=T, col.names=T)

    dat %>%
      slice(-trainIx) %>%
      cbind(., predTest) %>%
      write.csv(., paste0(rf,"/predictions_", m, ".csv"), quote=F, row.names=F)
  }

}

# Create data partition and run models --------------------------------------------------
if(cvmethod == "chromosome"){

  foldvec <- as.numeric(as.factor(dat$promoter.chrom))
  chrs <- unique(foldvec)
  chrs <- 1:4
  for(i in seq(1,length(chrs),2)){ # run for every chromosome 
    print(i)
    valchrom <- chrs[c(i, i+1)]
    trainIndex <- as.matrix(which(!foldvec %in% valchrom))
    resfold <- paste0(fold,"/chr",i,"_", transf)
    if(file.exists(paste0(resfold,"/model_performance.txt"))){file.remove(paste0(resfold,"/model_performance.txt"))}
    dir.create(resfold, recursive=T, showWarnings = F)
    run_models(trainIndex, resfold)
  }

} else if(cvmethod == "promoter"){
  
  trainFraction <- 0.8
  proms <- unique(dat$promoter.name)
  nproms <- length(proms)
  trainproms <- proms[sample(1:nproms, nproms * trainFraction, replace=FALSE)]
  trainIndex <- as.matrix(which(dat$promoter.name %in% trainproms))

  resfold <- paste0(fold,"/prom_", transf)
  if(file.exists(paste0(resfold,"/model_performance.txt"))){file.remove(paste0(resfold,"/model_performance.txt"))}
  run_models(trainIndex, resfold)

} else {

  trainFraction <- 0.8
  trainIndex <- createDataPartition(dat[,label], p = trainFraction, list = FALSE)
  resfold <- paste0(fold,"/rand_", transf)
  dir.create(resfold, recursive=T, showWarnings = F)
  if(file.exists(paste0(resfold,"/model_performance.txt"))){file.remove(paste0(resfold,"/model_performance.txt"))}
  run_models(trainIndex, resfold)

}

stopCluster(cl)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6941312/

