#### GARBAGE BIN ######
## function to add atac seq data to the training data #
add_atac <- function(i){
  print(i)
  # enhancer
  ind <- which(atac$chrom == tf_train[i, enhancer_chrom] & 
                 ( (atac$chromStart >= tf_train[i, enhstart] & atac$chromStart <= tf_train[i, enhend]) | 
                     (atac$chromEnd >= tf_train[i, enhstart] & atac$chromEnd <= tf_train[i, enhend]) |
                     (atac$chromStart <= tf_train[i, enhstart] & atac$chromEnd >= tf_train[i, enhend])
                 )
  )
  
  tf_train[i, atac_enhancer := mean(atac[ind, signalValue], na.rm=T) ]
  
  # promoter
  ind2 <- which(atac$chrom == tf_train[i, promoter_chrom] & 
                  ( (atac$chromStart >= tf_train[i, promstart] & atac$chromStart <= tf_train[i, promend]) | 
                      (atac$chromEnd >= tf_train[i, promstart] & atac$chromEnd <= tf_train[i, promend]) |
                      (atac$chromStart <= tf_train[i, promstart] & atac$chromEnd >= tf_train[i, promend])
                  )
  )
  
  tf_train[i, atac_promoter := mean(atac[ind2, signalValue], na.rm=T) ]
  
  # window 
  ind3 <- which(atac$chrom == tf_train[i, window_chrom] & 
                  ( (atac$chromStart >= tf_train[i, windstart] & atac$chromStart <= tf_train[i, windend]) | 
                      (atac$chromEnd >= tf_train[i, windstart] & atac$chromEnd <= tf_train[i, windend]) |
                      (atac$chromStart <= tf_train[i, windstart] & atac$chromEnd >= tf_train[i, windend])
                  )
  )
  
  tf_train[i, atac_window := mean(atac[ind3, signalValue], na.rm=T) ]
}

sapply(1:nrow(tf_train), add_atac)