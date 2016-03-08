library(XML)
library(plyr)

chond_rediff_SAS_agg <- data.frame()
chond_rediff_SAS_agg <- rbind(chond_rediff_SAS_agg, c(0, 1, 5.4453, 0.03338, 325, 163.15, '<.0001', 231.66, 7.7319))
chond_rediff_SAS_agg[,c(1:9)] <- sapply(chond_rediff_SAS_agg[,c(1:9)], as.character)
chond_rediff_SAS_agg <- rbind(chond_rediff_SAS_agg, c(5, 1, 5.6317, 0.03723, 353, 151.28, '<.0001', 279.13, 10.3911))
chond_rediff_SAS_agg <- rbind(chond_rediff_SAS_agg, c(0, 14, 5.4175, 0.03565, 325, 151.97, '<.0001', 225.32, 8.0327))
chond_rediff_SAS_agg <- rbind(chond_rediff_SAS_agg, c(5, 14, 5.0266, 0.05285, 353, 95.10, '<.0001', 152.41, 8.0552))
names(chond_rediff_SAS_agg) <- c('Passage', 'Days_CM', 'Estimate', 'Standard Error', 'DF', 't Value', 'Pr > |t|', 'Mean', 'StandardErrorMean')
chond_rediff_SAS_agg <- cbind(chond_rediff_SAS_agg, 'experiment' = 'Chondrocyte gel - glmm with random residuals -- aggrecan')






files <- list.files(path = '~/Dropbox/cow fishing/cowFISH_figures/Statistics',
                  pattern="^(?i).{1,}(?:random residuals).{1,}.html$",recursive=TRUE,full.names=TRUE)
filecodes <- str_split( files, pattern = '/|\\.')

allSAS <- data.frame()
for(i in seq(along = files) ) {
  
  x <- readHTMLTable(files[i])
  x <- x[[12]]
  experiment <- filecodes[[i]][length(filecodes[[i]])-1]
  x <- cbind( x, experiment )
  rm(experiment)
#   allSAS[[i]] <- x 
  allSAS <- rbind.fill(allSAS, x)
}

allSAS <- rbind.fill(allSAS, chond_rediff_SAS_agg)

allSAS$Mean <- as.numeric(as.character(allSAS$Mean))
allSAS$StandardErrorMean <- as.numeric(as.character(allSAS$StandardErrorMean))
allSAS$Days_CM <- as.numeric(as.character(allSAS$Days_CM))
allSAS$Passage <- as.numeric(as.character(allSAS$Passage))

msc_correlation_SAS <- readHTMLTable('~/Dropbox/cow fishing/cowFISH_figures/Statistics/MSC gel correlations - linear mixed model with random intercept.html')
msc_correlation_SAS <- msc_correlation_SAS[[11]]
msc_correlation_SAS$Days_CM <- as.numeric(as.character(msc_correlation_SAS$Days_CM))
msc_correlation_SAS$Estimate <- as.numeric(as.character(msc_correlation_SAS$Estimate))
msc_correlation_SAS <- rename(msc_correlation_SAS, c('Standard Error' = 'Standard_Error'))
msc_correlation_SAS$Standard_Error <- as.numeric(as.character(msc_correlation_SAS$Standard_Error))

chond_correlation_SAS <- readHTMLTable('~/Dropbox/cow fishing/cowFISH_figures/Statistics/Chondrocyte gel correlations - linear mixed model with random intercept.html')
chond_correlation_SAS <- chond_correlation_SAS[[13]]
chond_correlation_SAS$Days_CM <- as.numeric(as.character(chond_correlation_SAS$Days_CM))
chond_correlation_SAS$Passage <- as.numeric(as.character(chond_correlation_SAS$Passage))
chond_correlation_SAS$Estimate <- as.numeric(as.character(chond_correlation_SAS$Estimate))
chond_correlation_SAS <- rename(chond_correlation_SAS, c('Standard Error' = 'Standard_Error'))
chond_correlation_SAS$Standard_Error <- as.numeric(as.character(chond_correlation_SAS$Standard_Error))

rm(x, i, files, filecodes)

#### ----- single import
# SAS_MSC_gel_glmm_aggrecan <- readHTMLTable('~/Dropbox/Cow FISHing/cowFISH_figures/Statistics/MSC gel - glmm with random residuals -- Aggrecan.html')
# SAS_MSC_gel_glmm_aggrecan <- SAS_MSC_gel_glmm_aggrecan[[12]]
# SAS_MSC_gel_glmm_aggrecan$Days_CM <- as.numeric(as.character(SAS_MSC_gel_glmm_aggrecan$Days_CM))
# SAS_MSC_gel_glmm_aggrecan$Mean <- as.numeric(as.character(SAS_MSC_gel_glmm_aggrecan$Mean))
# SAS_MSC_gel_glmm_aggrecan$StandardErrorMean <- as.numeric(as.character(SAS_MSC_gel_glmm_aggrecan$StandardErrorMean))

