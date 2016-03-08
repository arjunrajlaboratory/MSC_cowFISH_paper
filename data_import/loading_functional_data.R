# loading in functional data and separating into groups
# library(stringr)
# library(ggplot2)
# library(reshape2)
# library(RCurl)
# library(plyr)


gag_per_pellet <- read.csv('./all_csvs/gag_per_pellet.csv', header = TRUE)
# gag_per_pellet <- transform(gag_per_pellet, Sample = reorder(Sample, Media, Mean))
# gag_per_pellet <- gag_per_pellet[with(gag_per_pellet, order(Media, Mean)), ]

gag_per_dna <- read.csv('./all_csvs/gag_per_DNA.csv', header = TRUE)
# gag_per_dna <- transform(gag_per_dna, Sample = reorder(Sample, Mean))
# gag_per_dna_melt <- melt(gag_per_dna, id.vars = c('Sample', 'Media'))
# gag_per_dna_cast <- dcast(gag_per_dna_melt, 'Sample ~ Media + variable', value.var = 'value')
# gag_per_dna_cast <- rename(gag_per_dna_cast, c('Sample' = 'Sample','CM-_Mean' = 'CMminus_Mean', 'CM-_StDev' ='CM_minus_StDev',
#                                                'CM+_Mean'  = 'CMplus_Mean', 'CM+_StDev' = 'CMplus_StDev'))
# gag_per_dna_cast <- transform(gag_per_dna_cast, Sample = reorder(Sample, CMplus_Mean - CMminus_Mean))
# 
# 
# 
# PCA_clusters <- data.frame(Sample = c('A', 'B', 'C', 'chondrocyte', 'D', 'E', 'F', 'G', 'H', 'Het', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q'), 
#                            cluster = c('cluster5',
#                    'cluster5', 
#                    'cluster3', 
#                    'chondrocyte', 
#                    'cluster5', 
#                    'cluster5', 
#                    'cluster5', 
#                    'cluster3', 
#                    'cluster5',
#                    'Het', 
#                    'cluster5', 
#                    'cluster5', 
#                    'cluster5', 
#                    'cluster3', 
#                    'cluster5', 
#                    'cluster2', 
#                    'cluster3', 
#                    'cluster2', 
#                    'cluster2'))
# 
# gag_per_dna_cluster <- merge(gag_per_dna, PCA_clusters)
# gag_per_dna_cluster <- gag_per_dna_cluster[order(gag_per_dna_cluster$cluster),]
# 
# gag_per_pellet_cluster <- merge(gag_per_pellet, PCA_clusters)
# gag_per_pellet_cluster <- gag_per_pellet_cluster[order(gag_per_pellet_cluster$cluster),]
# 
# 
# 
# theme_set(theme_bw(base_size = 18))
# pdf('~/Dropbox/Personal/Ally Lab/Ally Presentation/2014-10-24_CBP_RIP/rnaseq_function.pdf')
# ggplot(gag_per_dna, aes(x = as.factor(Sample), y = Mean, fill = Media)) + 
#   geom_bar(stat = 'identity') + 
#   geom_bar(stat = 'identity', color = '#000000', show_guide = FALSE) + 
#   geom_errorbar(aes(x = as.factor(Sample), ymax = Mean + StDev, ymin = Mean - 0.0001), width = 0.5) + 
#   facet_grid(Media ~ .) + xlab('Sample') + ylab('GAG/DNA') + 
#   theme(legend.position = c(1,1), legend.justification = c(1,1)) + 
#   scale_fill_manual(values = c('mediumaquamarine', 'turquoise4'), labels = c(expression(paste('-TGF', beta)), 
#                                                                              expression(paste('+TGF', beta))))
# 
#   ggplot(gag_per_pellet, aes(x = as.factor(Sample), y = Mean, fill = Media)) + 
#   geom_bar(stat = 'identity') + 
#   geom_bar(stat = 'identity', color = '#000000', show_guide = FALSE) + 
#   geom_errorbar(aes(x = as.factor(Sample), ymax = Mean + StDev, ymin = Mean - 0.0001), width = 0.5) + 
#   theme(legend.position = c(1,1), legend.justification = c(1,1)) + 
#   scale_fill_manual(values = c('mediumaquamarine', 'turquoise4'), labels = c(expression(paste('-TGF', beta)), 
#                                                                              expression(paste('+TGF', beta)))) +
#   facet_grid(Media ~ .) + xlab('Sample') + ylab('GAG/pellet')
# dev.off()
# 
# 
# ggplot(gag_per_dna_cast, aes(x = as.factor(Sample), y = CMplus_Mean - CMminus_Mean)) + 
#   geom_bar(stat = 'identity')
# 
# ggplot(gag_per_dna_cast, aes(x = as.factor(Sample), y = CMplus_Mean/CMminus_Mean)) + 
#   geom_bar(stat = 'identity')
# 
# 
# ggplot(subset(gag_per_pellet, Media == 'CM+'), aes(x = as.factor(Sample), y = Mean)) + 
#   geom_bar(stat = 'identity') +
#   geom_errorbar(aes(x = as.factor(Sample), ymax = Mean + StDev, ymin = Mean - (StDev)))
# 
# ggplot(gag_per_pellet, aes(x = as.factor(Sample), y = Mean)) + 
#   geom_bar(stat = 'identity') + facet_grid(Media ~ .) +
#   geom_errorbar(aes(x = as.factor(Sample), ymax = Mean + StDev, ymin = Mean - (StDev)))
# 
# 
# ## plots with clustering information-------------------------
# theme_set(theme_bw(base_size = 18))
# clustercolors <- c('palegreen', 'mediumaquamarine', 'seagreen4', 'darkslategray')
# 
# pdf('~/Dropbox/Personal/Ally Lab/Ally Presentation/2014-11-07_labmeeting/rna_seq_function_cluster.pdf')
# ggplot(gag_per_dna_cluster, aes(x = Sample, y = Mean, fill = cluster)) + 
#   geom_bar(stat = 'identity') + 
#   geom_bar(stat = 'identity', color = '#000000', show_guide = FALSE) + 
#   geom_errorbar(aes(x = Sample, ymax = Mean + StDev, ymin = Mean - 0.0001), width = 0.5) + 
#   facet_grid(Media ~ .) + xlab('Sample') + ylab('GAG/DNA') + 
#   theme(legend.position = c(1,1), legend.justification = c(1,1)) + scale_x_discrete(limits = gag_per_dna_cluster$Sample) +
#   scale_fill_manual(name = 'PCA Cluster', labels = c('Group 1', 'Group 2', 'Group 3', 'Het'), values = clustercolors)
# 
# ggplot(gag_per_pellet_cluster, aes(x = as.factor(Sample), y = Mean, fill = cluster)) + 
#   geom_bar(stat = 'identity') + 
#   geom_bar(stat = 'identity', color = '#000000', show_guide = FALSE) + 
#   geom_errorbar(aes(x = as.factor(Sample), ymax = Mean + StDev, ymin = Mean - 0.0001), width = 0.5) + 
#   facet_grid(Media ~ .) + xlab('Sample') + ylab('GAG/pellet') + 
#   theme(legend.position = c(1,1), legend.justification = c(1,1)) + scale_x_discrete(limits = gag_per_pellet_cluster$Sample) +
#   scale_fill_manual(name = 'PCA Cluster', labels = c('Group 1', 'Group 2', 'Group 3', 'Het'), values = clustercolors)
# dev.off()

# if(!exists('microdata3')) source('~/code/cowfish/Old_Data_Import_Scripts/micromechanical_data_compile.R')
# microdata_clustering <- merge(x = microdata3, y = PCA_clusters, by.x = 'clone', by.y = 'Sample')
# microdata_clustering <- microdata_clustering[order(microdata_clustering$cluster),]
# ggplot(subset(microdata_clustering,day==7 & `30per`>=1 & `0per`>=1),aes(x=as.factor(clone),y=(`30per`-`0per`), fill = cluster))+geom_boxplot(width = 1)+
#   xlab('Clone')+ylab('AR at 30% Strain - AR at 0% Strain') + scale_x_discrete(limits = microdata_clustering$clone)
# 

