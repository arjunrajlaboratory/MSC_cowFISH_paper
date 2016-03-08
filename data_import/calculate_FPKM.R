##calculating FPKM for cow genome RNA-seq results
library(plyr)
library(dplyr)
library(data.table)

countTable <- read.table('./all_csvs/CowFishClone_HTSeqCounts_2015-10-08_20-02.tsv', header = TRUE)
gene_info <- read.table('./all_csvs/2015-10-07_complete_feature_file.symbolsonly.gtf')
gene_info_mod <- select(gene_info, V1, V3, V4, V5, V7, V8, V10, V13,  V19, V31)
names(gene_info_mod) <- c('Chromosome', 'Feature', 'StartPosition', 'EndPosition', 'Strand', 'Codon', 'transcript_number',
                          'gene_id', 'transcript_id', 'gene_name')

# gene_info_mod <- subset(gene_info_mod, select = c(Chromosome, Source, Feature, StartPosition, EndPosition, Strand, Codon, gene_id, 
#                                                   transcript_id, exon_number, exon_id, gene_name)) #remove empty columns
gene_info_mod <- subset(gene_info_mod, Feature == 'exon') #select only exons
gene_info_mod$feature.length <- c(gene_info_mod$EndPosition - gene_info_mod$StartPosition)


## ------ calculating gene length

calcTotalExonLength <- function(start,end) {
  #subDat <- data[,c('start','end')]
  subDat <- cbind(start,end)
  subDat <- data.frame(subDat)
  colnames(subDat) <- c('start','end')
  
  orderDat <- subDat[order(subDat$start),]
  orderDat <- unique(orderDat)
  
  geneLength <- 0
  for (i in 1:nrow(orderDat)) {
    if (i == 1) {
      geneLength <- orderDat[1,2] - orderDat[1,1]
    } else if (orderDat[i,1] >= orderDat[i-1,2]) { # this exon does not overlap previous at all
      geneLength <- geneLength + (orderDat[i,2] - orderDat[i,1])
    } else if (orderDat[i,2] > orderDat[i-1, 2]) { # this exon starts in the middle of previous, but extends beyond it
      geneLength <- geneLength + (orderDat[i,2] - orderDat[i-1,2])
    }
  }
  
  return(geneLength)
}

gene_info_mod <- data.table(gene_info_mod)
gene_length <- gene_info_mod[ , list(exontotal = calcTotalExonLength(StartPosition,EndPosition)), by = 'gene_name']
gene_length <- as.data.frame(gene_length)
gene_length$exontotal <- gene_length$exontotal/1000

###---------

mapped_size <- ddply(countTable, .(sampleID), summarize, million.reads = sum(counts/10^6))

FPKM_table <- merge(countTable, mapped_size, by = 'sampleID')
FPKM_table <- merge(x = FPKM_table, y = gene_length, by.x = 'gene_id',by.y = 'gene_name', all.x = TRUE )
FPKM_table$FPKM <- c(FPKM_table$counts/FPKM_table$exontotal/FPKM_table$million.reads)

rm(calcTotalExonLength, gene_info, gene_info_mod, gene_length, mapped_size)
