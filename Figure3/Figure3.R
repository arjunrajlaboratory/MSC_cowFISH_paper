library(stringr)
library(plyr)
library(data.table)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(RColorBrewer)
library(DESeq)
library(gplots)
library(ggfortify)


if(!exists('FPKM_table')) source ('./data_import/calculate_FPKM.R')
if(!exists('gag_per_pellet')) source ('./data_import/loading_functional_data.R')
theme_set(theme_bw(base_size = 8))

gag_per_dna <- filter(gag_per_dna, Media == 'CM+')
gag_per_dna <- arrange(gag_per_dna, Sample)
gag_per_dna <- arrange(gag_per_dna, Mean)
FPKM_table <- filter(FPKM_table, sampleID %in% gag_per_dna$Sample)

FPKM_forclustering <- dcast(FPKM_table, 'gene_id ~ sampleID', value.var = 'FPKM')
FPKM_forclustering <- filter(FPKM_forclustering, B != 'NA')
row.names(FPKM_forclustering) <- FPKM_forclustering$gene_id
FPKM_forclustering$gene_id <- NULL
FPKM_forclustering <- t(FPKM_forclustering)
FPKM_forclustering <- FPKM_forclustering[,apply(FPKM_forclustering, 2, function(x) any(x>=1))]
FPKM_forclustering <- t(FPKM_forclustering)

FPKM_selectgenes <- mutate(FPKM_table, log2FPKM = log2(FPKM + 0.25))
FPKM_selectgenes <- dcast(FPKM_selectgenes, 'sampleID ~ gene_id', value.var = 'log2FPKM')


### heatmap --------------------
pearsonDist <- function(matdata){
  return(as.dist(1-cor(t(matdata))))
}
rowDistance = as.dist(1-cor((FPKM_forclustering)))
colDistance = as.dist(1-cor(t(FPKM_forclustering)))
rowCluster = hclust(rowDistance)
colCluster = hclust(colDistance)

pearsonCompleteClustering <- function(x){
  hclust(pearsonDist(x))
}
geneClusterer <- pearsonCompleteClustering
geneClustering <- geneClusterer(FPKM_forclustering)
color4heatmap <- brewer.pal(11,'RdYlBu')

seq_heatmap_nochondrocytes <- 
heatmap.2(as.matrix(FPKM_forclustering), 
          scale = 'row', 
          dendrogram = c('column'),
          Colv = as.dendrogram(rowCluster),
#           Colv =  c('J', 'Q', 'C', 'Het', 'P', 'D', 'L', 
#                     'G', 'E', 'M', 'O', 'K', 
#                     'I', 'B', 'H', 'F'),
          Rowv = as.dendrogram(colCluster),
          col = color4heatmap,
          # labCol=colnames(FPKM_forclustering),
          labRow=FALSE,
          density.info='none', 
          symkey=FALSE, trace='none',
          margins = c(2,2), cexRow=1,
          srtCol = 0, adjCol = c(.5, 0))

### PCA --------------------------------
FPKM_PCA <- prcomp(t(FPKM_forclustering), center = TRUE, scale. = TRUE)
# summary(FPKM_PCA)

PCAplot <- 
  autoplot(FPKM_PCA, shape = FALSE, label.size = 8, data = gag_per_dna, colour = 'Mean') + 
  scale_colour_gradient2(mid = 'sienna1',  high = 'cornflowerblue', name = 'Mean\nGAG/DNA') +
  theme(legend.position = c(1,1), legend.justification = c(1,1))

### functional data ----------------------------------------
gag_per_dna_graph <- 
  ggplot(gag_per_dna, aes(x = as.factor(Sample), y = Mean)) + 
  geom_bar(stat = 'identity', color = '#000000', show_guide = FALSE, fill = 'mediumaquamarine') + 
  geom_errorbar(aes(x = as.factor(Sample), ymax = Mean + StDev, ymin = Mean - 0.0001), width = 0.5) +
  xlab('Sample') + ylab('GAG/DNA') + scale_x_discrete(limits = c('M', 'Het', 'F', 'B', 'E', 'D', 'K', 
                                                                    'P', 'Q', 'L', 'G', 'C', 'O', 'J', 'H', 'I'))

plot(gag_per_dna_graph)

### selected genes vs. functional data --------------------
select_genes <- c('ACAN', 'COMP', 'SOX9', 'SPP1',
                  'COL2A1', 'COL1A2', 'RUNX2', 'ALPL', 'MMP13', 'TGFBR1', 
                  'TEK', 'ITGAV', 'THY1', 'CD63', 'CD200', 'PTPRC', 
                  'FGF2', 'SOX2', 'ENG', 
                  # 'POU5F1', 
                  'CDK2',
                  'ACTB', 'LDHA', 'GAPDH',
                  'CCNA2')
# missing from annotation: 
# alt names: CD105 = ENG, CD45 = PTPRC, TIE2 = TEK

seq_genes_of_interest <- select(FPKM_selectgenes, one_of(c(select_genes, 'sampleID')))
seq_genes_of_interest <- merge(x = seq_genes_of_interest, y = gag_per_dna, by.x = 'sampleID', by.y = 'Sample')
seq_genes_of_interest <- melt(seq_genes_of_interest, id.vars = c('sampleID', 'Mean', 'StDev', 'Media'))


seq_genes_of_interest_cor_test <- 
  seq_genes_of_interest %>% group_by(variable) %>% summarize(testvalue = cor.test(Mean, value)$p.value)

lm_eqn <- function(df){
  m <- lm(value ~ Mean, df);
  eq <- substitute(~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

genes_with_functional_data <- 
ggplot(seq_genes_of_interest, aes(x = Mean, y = value)) + geom_point(color = 'mediumaquamarine') + 
  facet_wrap(~variable, scale = 'free_y', ncol = 4) + xlab('Mean GAG/DNA') + ylab('log2(FPKM)') +
  # geom_smooth(method = 'lm', color = 'black', se = FALSE) +
  theme(strip.background = element_rect(colour="white", fill="white")) 
# +
#   geom_text(data=ddply(seq_genes_of_interest, .(variable), lm_eqn),
#             aes(x = 200, y = 0, label=V1), parse = TRUE, inherit.aes=FALSE)
plot(genes_with_functional_data)









### final plotting --------------
plots_square <- c('gag_per_dna_graph', 'PCAplot')
plots_rectangle_portrait <- c('seq_heatmap_nochondrocytes')
plots_rectangle_landscape <- c('genes_with_functional_data')

directory <- paste('./Figure3/graphs/', Sys.Date(),'/', sep = '')
if(directory %in% dir() == FALSE) dir.create(directory) 
for (i in seq(along = plots_square)){
  savefile <- paste(Sys.Date(),'_',as.name(plots_square[i]),'.pdf',sep = '')
  ggsave(filename = savefile, plot = get(plots_square[i]), device = pdf,
         path = directory, scale = 1, 
         width = 10, height = 10, units = "cm", dpi = 300)
}

# for (i in seq(along = plots_rectangle_portrait)){
   savefile <- paste(Sys.Date(),'_','heatmap','.pdf',sep = '')
   pdf(paste(directory,savefile, sep = ''))
    # get(plots_rectangle_portrait[i])  
   heatmap.2(as.matrix(FPKM_forclustering), 
             scale = 'row', 
             dendrogram = c('column'),
             Colv = as.dendrogram(rowCluster),
             Rowv = as.dendrogram(colCluster),
             col = color4heatmap,
             labCol=colnames(FPKM_forclustering),labRow=FALSE,
             density.info='none', 
             symkey=FALSE, trace='none',
             margins = c(2,2), cexRow=1,
             srtCol = 0, adjCol = c(.5, 0))
   
   dev.off()
   rm(savefile)
 # }

for (i in seq(along = plots_rectangle_landscape)){
  savefile <- paste(Sys.Date(),'_',as.name(plots_rectangle_landscape[i]),'.pdf',sep = '')
  ggsave(filename = savefile, plot = get(plots_rectangle_landscape[i]), device = pdf,
         path = directory, scale = 1, 
         width = 20, height = 30, units = "cm", dpi = 300)
}


# transposed heatmap:
#   heatmap.2(t(as.matrix(FPKM_forgraphs)), 
#             scale = 'column', 
#             Colv = as.dendrogram(colCluster),
#             Rowv = as.dendrogram(rowCluster),
#             col = color4heatmap,
#             labCol=FALSE,labRow=colnames(FPKM_forgraphs),
#             density.info='none', 
#             symkey=FALSE, trace='none',
#             margins = c(2,10), cexCol=1, 
#             main='pearson row scaled')
