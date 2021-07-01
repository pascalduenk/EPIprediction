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
dat <- fread("../../dataTrainTest.csv", header=T, sep="\t")
  
pred <- fread("predictions_rf.csv", header=F)
pred <- pred %>%
    mutate(V1 := paste0(V1,", ", V2)) %>%
    select(-V2)
colnames(pred) <- c(colnames(dat),"predicted")

# get training data 
train <- anti_join(dat, pred, by=c("promoter.name","enhancer.name"))

# correct or wrong prediction 
pred$prom.intrain <- pred$promoter.name %in% train$promoter.name

pred %>%
    group_by(predicted, label, prom.intrain) %>%
    count() %>%
    write.table(., "correct_wrong_train.txt", row.names=F, quote=F)

# check correct true predictions. Promoter in train all true?
pos.prom <- pred %>%
    filter(prom.intrain == TRUE & (predicted == label) & predicted == 1) %>%
    pull(promoter.name)

train %>%
    filter(promoter.name %in% pos.prom) %>%
    group_by(promoter.name, label) %>%
    count()

false.prom <- pred %>%
    filter(prom.intrain == TRUE & (predicted != label) & predicted == 1) %>%
    pull(promoter.name)

train %>%
    filter(promoter.name %in% false.prom) %>%
    select(promoter.name, label) %>%
    with(.,table(promoter.name, label))

train %>%
    filter(promoter.name %in% pos.prom) %>%
    select(promoter.name, label) %>%
    with(.,table(promoter.name, label))