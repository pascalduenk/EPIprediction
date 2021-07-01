# original CHiP-seq data from ENCODE - cell line K562 #
cnames_chipseq <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")

K4me1 <- fread("CHiPseq/ENCFF590NGQ.bed", header=F, col.names = cnames_chipseq) %>%
  setorder(chrom, chromStart)
K4me3 <- fread("CHiPseq/ENCFF616DLO.bed", header=F, col.names = cnames_chipseq)  %>%
  setorder(chrom, chromStart)
K27me3 <- fread("CHiPseq/ENCFF031FSF.bed", header=F, col.names = cnames_chipseq) %>%
  setorder(chrom, chromStart)
K27ac <- fread("CHiPseq/ENCFF038DDS.bed", header=F, col.names = cnames_chipseq) %>%
  setorder(chrom, chromStart)

# combine active enhancers with active promoters
setkey(K4me1, chrom, chromStart, chromEnd)
setkey(K4me3, chrom, chromStart, chromEnd)
setkey(K27me3, chrom, chromStart, chromEnd)
setkey(K27ac, chrom, chromStart, chromEnd)

act_enh <- foverlaps(K27ac, K4me1, type="any", nomatch = NULL)
act_prom <- foverlaps(K27ac, K4me3, type="any", nomatch = NULL)

# both <- foverlaps(K4me1, K4me3, type="any", nomatch=NULL)

# atac seq #
cnames_atac <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")
atac <- fread("ATACseq/K562/ENCFF678GBM.bed", header=F, col.names = cnames_atac) %>%
  setorder(chrom, chromStart)

setkey(atac, chrom, chromStart, chromEnd)
setkey(act_prom, chrom, chromStart, chromEnd)
setkey(act_enh, chrom, chromStart, chromEnd)

prom_atac <- foverlaps(act_prom[,1:10], atac, type="any")
enh_atac <- foverlaps(act_enh[,1:10], atac, type="any")

prom_atac %>%
  filter(chrom == "chr6") %>%
  with(unique(promoter))
