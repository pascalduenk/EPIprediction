setwd("C:/Users/duenk002/OneDrive - WageningenUR/Projects/11 EPI prediction")
library(dplyr)
library(rtracklayer)
library(data.table)
library(tidyr)
library(caret)

options(stringsAsFactors = F)
selchrom <- "6"

# ATAC data 
cnames_atac <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")
atac <- fread("ATACseq/IMR90/ENCFF830JMV_37.bed", header=F, col.names = cnames_atac) %>%
  #filter(chrom == selchrom) %>%
  setorder(chrom, chromStart) %>%
  mutate(chrom := paste0("chr", chrom))

# Training data of TargetFinder
tf_train <- fread("TargetFinder/training.csv", header=T)
tf_train <- tf_train %>%
  select(bin, enhancer_chrom, enhancer_distance_to_promoter, enhancer_end, enhancer_name, enhancer_start,
         label, promoter_chrom, promoter_end, promoter_name, promoter_start, window_end, window_start, window_chrom, window_name,
         interactions_in_window, active_promoters_in_window, `H3K27ac (enhancer)`, `H3K27ac (promoter)`, `H3K27ac (window)`,
         `H3K27me3 (enhancer)`, `H3K27me3 (promoter)`,`H3K27me3 (window)`, `H3K4me1 (enhancer)`,`H3K4me1 (promoter)`,`H3K4me1 (window)`,
         `Methylation (enhancer)`, `Methylation (promoter)`, `Methylation (window)`) %>%
  #filter(promoter_chrom == paste0('chr',selchrom)) %>%
  setorder(promoter_chrom, promoter_start)

# write files for Ole #
# tf_train %>%
#   select(promoter_chrom, promoter_start, promoter_end, promoter_name) %>%
#   distinct(promoter_chrom, promoter_name, promoter_start, promoter_end) %>%
#   write.table("promoters.bed", row.names = F, col.names = F, quote = F, sep = "\t")
# 
# tf_train %>%
#   select(enhancer_chrom, enhancer_start, enhancer_end, enhancer_name) %>%
#   distinct(enhancer_chrom, enhancer_name, enhancer_start, enhancer_end) %>%
#   write.table("enhancers.bed", row.names = F, col.names = F, quote = F, sep = "\t")

# Gene expression (Total RNA)
rna_seq <- fread("RNAseq/IMR90/ENCFF353SBP.tsv")

# reference genome
refg <- rtracklayer::import("Reference/Homo_sapiens.GRCh37.87.gtf")
ref_df <- as.data.table(refg) %>% filter(type == "gene")
rm(refg)
gc(reset=T)

# combine rna expression with reference genome
rna_ref <- ref_df %>%
  select(seqnames, start, end, width, strand, gene_id, gene_version, gene_name) %>%
  mutate(gene_id2 = paste(gene_id, gene_version, sep=".")) %>%
  left_join(rna_seq, by=c("gene_id2" = "gene_id"))

# count rows 
rna_ref %>%
  filter(!is.na(TPM)) %>%
  nrow()

# write file for Ole 
rna_ref %>%
  write.table("rna_seq.bed", row.names = F, col.names = F, quote = F, sep = "\t")

######
# combine atac seq with training data 
atac <- atac %>% 
  setorder(chrom, chromStart)

# write file for Ole
atac %>%
  write.table("atac.bed", row.names = F, col.names = F, quote = F, sep = "\t")

# make the regions of promoters and enhancers a bit larger 
tf_train <- tf_train %>%
  mutate(promstart := promoter_start - 100, promend := promoter_end + 100,
         enhstart := enhancer_start - 500, enhend := enhancer_end + 500,
         windstart := ifelse(promstart > enhstart, enhend + 1,
                             promend + 1),
         windend := ifelse(promstart > enhstart, promstart -1,
                             enhstart - 1))

# make tables with unique enhancers, promoters and windows #
unique_enh <- tf_train %>%
  distinct(enhancer_chrom, enhancer_name, enhstart, enhend) %>%
  setkey(enhancer_chrom, enhstart, enhend)

unique_prom <- tf_train %>%
  distinct(promoter_chrom, promoter_name, promstart, promend) %>%
  setkey(promoter_chrom, promstart, promend)

unique_wind <- tf_train %>%
  distinct(window_chrom, window_name, windstart, windend) %>%
  setkey(window_chrom, windstart, windend)

# combine with ATAC data #
atac_enh <- atac %>% 
  select(chrom, chromStart, chromEnd, name, signalValue) %>% 
  setkey(chrom, chromStart, chromEnd) %>%
  foverlaps(unique_enh, ., type="any", nomatch = NA) %>%
  group_by(enhancer_name) %>%
  summarize(atac_enhancer = sum(signalValue, na.rm=T))

atac_prom <- atac %>%
  select(chrom, chromStart, chromEnd, name, signalValue) %>% 
  setkey(chrom, chromStart, chromEnd) %>%
  foverlaps(unique_prom, ., type="any", nomatch = NA) %>%
  group_by(promoter_name) %>%
  summarize(atac_promoter = sum(signalValue, na.rm=T))

atac_wind <- atac %>%
  select(chrom, chromStart, chromEnd, name, signalValue) %>% 
  setkey(chrom, chromStart, chromEnd) %>%
  foverlaps(unique_wind, ., type="any", nomatch = NA) %>%
  group_by(window_name) %>%
  summarize(atac_window = sum(signalValue, na.rm=T))

# combine with the training data 
tf_train <- tf_train %>%
  left_join(atac_enh, by="enhancer_name") %>%
  left_join(atac_prom, by="promoter_name") %>%
  left_join(atac_wind, by="window_name")

# adjust atac values according to size #
tf_train <- tf_train %>%
  mutate(windsize := windend - windstart,
         promsize := promend - promstart,
         enhsize := enhend - enhstart) %>%
  mutate(atac_window := (atac_window/windsize),
         atac_promoter := (atac_promoter/promsize),
         atac_enhancer := (atac_enhancer/enhsize))

# # check how many promoters have ATAC values
# tf_train %>%
#   filter(is.na(atac_promoter)) %>%
#   with(unique(promoter_name)) %>%
#   length()
# 
# # check how many enhancers have ATAC values
# tf_train %>%
#   filter(is.na(atac_enhancer)) %>%
#   with(unique(enhancer_name)) %>%
#   length()
# 
# # check how many windows have ATAC values
# tf_train %>%
#   filter(is.na(atac_window)) %>%
#   with(unique(window_name)) %>%
#   length()

# add RNA seq to the training data
rna_prom <- tf_train %>%
  distinct(promoter_chrom, promoter_name, promoter_start, promoter_end)

rna_prom$gene_id <- character()
rna_prom$gene_start <- integer()
rna_prom$gene_end <- integer()
rna_prom$gene_strand <- character()
rna_prom$gene_expression <- numeric()

for(i in 1:nrow(rna_prom)){
  # print(i)
  # find closest gene on plus strand
  ps <- rna_prom[i, promoter_start]
  gn_plus <- rna_ref %>%
    filter(start >= ps & strand == "+") %>%
    slice_min(start) %>%
    filter(row_number()==1) %>%
    select(gene_id, start, end, strand, TPM)
  
  # find closest gene on min strand
  pe <- rna_prom[i, promoter_end]
  gn_min <- rna_ref %>%
    filter(end <= pe & strand == "-") %>%
    slice_max(end) %>%
    filter(row_number()==1) %>%
    select(gene_id, start, end, strand, TPM)
  
  d_plus <- gn_plus$start - pe
  d_min <- ps - gn_min$end
  
  if(d_plus < d_min){
    rna_prom[i,c('gene_id','gene_start','gene_end','gene_strand','gene_expression')] <- gn_plus
  } else {
    rna_prom[i,c('gene_id','gene_start','gene_end','gene_strand','gene_expression')] <- gn_min
  }
}

# compute distance between start site of promoter and start site of gene #
rna_prom <- rna_prom %>%
  mutate(d_prom_gene := case_when(gene_strand == "+" ~ gene_start - promoter_start,
                                  gene_strand == "-" ~ promoter_end - gene_end)) %>%
  mutate(gene_expression := case_when(is.na(gene_expression) ~ 0.00,
                                      !is.na(gene_expression) ~ gene_expression)) %>%
  select(promoter_name, gene_expression, d_prom_gene) 
  
plot(density(rna_prom$d_prom_gene, na.rm=T))

# combine with the training data 
tf_train <- tf_train %>%
  left_join(rna_prom, by="promoter_name")

# write results 
tf_train %>%
  write.table("data_complete.csv", row.names = F, col.names = T, quote = F, sep = "\t")


tf_train %>%
  group_by(label) %>%
  summarise(int = mean(interactions_in_window),
            int_sd = sd(interactions_in_window))
