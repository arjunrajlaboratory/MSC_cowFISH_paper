#Figures for paper
if(!exists('forscattersGAPDH10')) source('./data_import/calculate_all_stats_forpaper.R')
theme_set(theme_bw(base_size = 18))
systemType <- 'Mac'

# MSC GAPDH Expression
msc_mean_gapdh <- 
  ggplot(subset(allstats, Purpose == 'MSCs het differentiation in gels core' & Days.CM %in% c(1, 4, 7, 14, 21)), 
         aes(x = as.factor(Days.CM), y = GAPDH.RNA.mean, fill = Donor)) + 
  geom_errorbar(aes(ymin = GAPDH.RNA.mean - GAPDH.RNA.mean/20, ymax = GAPDH.RNA.mean + GAPDH.RNA.se, y = NULL), 
                width = 0.3, position=position_dodge(.5)) + 
  geom_bar(width = 0.5, stat = 'identity', position = position_dodge()) +
  scale_fill_manual(values = donor_colors_gray, labels = c('A', 'B', 'C')) +
  facet_grid(Media ~ .) + 
  theme(legend.position = c(1,1), legend.justification = c(1, 1), strip.background = element_blank(), strip.text = element_blank()) +
  geom_bar(data = mean_mscdiff, aes(x = as.factor(Days.CM), y = GAPDH.RNA.mean.mean), 
           fill = 'mediumseagreen', alpha = 0.6, stat = 'identity', color = '#000000', width = 0.7) +
  ylab('Mean GAPDH RNA') + xlab('Days in Agarose')

# MSC Live Dead Analysis and Graph
MSC_LD$Replicate.Number <- as.factor(MSC_LD$Replicate.Number)
MSC_LD$Day <- as.factor(MSC_LD$Day)
MSC_LD$liveFraction <- as.numeric(MSC_LD$X.Live...Filter)

viability_means<-ddply(MSC_LD, .(Replicate.Number, Media, Day), summarise, 
                       live=mean(liveFraction))
viability_means.means<-ddply(viability_means, .(Media, Day), summarise, 
                             liveM=mean(live),
                             sem=sd(live)/sqrt(length(live)))

msc_viability <- 
  ggplot()+
  geom_bar(data=viability_means,  
         aes(x = as.factor(Day), y = live, fill = as.factor(Replicate.Number)),
         stat='identity', position = position_dodge()) + 
  scale_fill_manual(values = donor_colors_gray, labels = c('A', 'B', 'C')) +
  geom_bar(data=viability_means.means,  
           aes(x = as.factor(Day), y = liveM),
           stat='identity', alpha=0.4, fill="mediumseagreen", color="black" )+
  geom_errorbar(data=viability_means.means,
                aes(x=as.factor(Day), ymin = liveM - sem, ymax = liveM + sem, y = NULL), 
                width = 0.3, position=position_dodge(.5)) + 
  facet_grid(Media ~ .) + 
  theme(legend.justification = c(1, 1), strip.background = element_blank()) +
  ylab('Fraction of Cells Alive') + xlab('Days in Agarose')


## fixable dead gapdh > 10 graph

isDead_plot_log <-  
  ggplot(subset(forscatters, experimentid %in% c('A68','A69')), aes(x = as.factor(isDead), y = GAPDH.RNA)) + geom_jitter(color = '#FF0000') + 
  facet_grid(Media~.) + scale_y_log10() +
  xlab('IsDead?') + ylab('GAPDH RNA')
isDead_plot <-  
  ggplot(subset(forscatters, experimentid %in% c('A68','A69')), aes(x = as.factor(isDead), y = GAPDH.RNA)) + geom_jitter(color = '#FF0000') + 
  facet_grid(Media~.) +
  xlab('IsDead?') + ylab('GAPDH RNA')


##printing plots
plots_square <- c('isDead_plot_log', 'isDead_plot')
plots_rectangle <- c('msc_viability', 'msc_mean_gapdh')


for (i in seq(along = plots_square)){
  savefile <- paste(Sys.Date(),'_',as.name(plots_square[i]),'.pdf',sep = '')
  ggsave(filename = savefile, plot = get(plots_square[i]), device = pdf,
         path = directory, scale = 1, 
         width = 10, height = 10, units = "cm", dpi = 300)
}

  
  if(systemType=="Mac"){
    
    directory <- paste('./FigureS1/graphs/', Sys.Date(),'/', sep = '')
    
    if(directory %in% dir() == FALSE) dir.create(directory) 
    for (i in seq(along = plots_square)){
      savefile <- paste(Sys.Date(),'_',as.name(plots_square[i]),'.pdf',sep = '')
      ggsave(filename = savefile, plot = get(plots_square[i]), device = pdf,
             path = directory, scale = 1, 
             width = 10, height = 10, units = "cm", dpi = 300)
    }
    
    for (i in seq(along = plots_rectangle)){
      savefile <- paste(Sys.Date(),'_',as.name(plots_rectangle[i]),'.pdf',sep = '')
      ggsave(filename = savefile, plot = get(plots_rectangle[i]), device = pdf,
             path = directory, scale = 1, 
             width = 15, height = 10, units = "cm", dpi = 300)
    }
    
  }


if(systemType=="Windows"){
  setwd("FigureS1\\")
  directory <- paste('graphs\\', Sys.Date(),'\\', sep = "")
  if(directory %in% dir() == FALSE) dir.create(directory)
  for (i in seq(along = plots_square)){
    savefile <- paste(Sys.Date(),'_',as.name(plots_square[i]),'.pdf',sep = '')
    ggsave(filename = savefile, plot = get(plots_square[i]), device = pdf,
           path = directory, scale = 1, 
           width = 10, height = 10, units = "cm", dpi = 300)
  }
  
  for (i in seq(along = plots_rectangle)){
    savefile <- paste(Sys.Date(),'_',as.name(plots_rectangle[i]),'.pdf',sep = '')
    ggsave(filename = savefile, plot = get(plots_rectangle[i]), device = pdf,
           path = directory, scale = 1, 
           width = 15, height = 10, units = "cm", dpi = 300)
  }
  
  setwd("..\\")
}
