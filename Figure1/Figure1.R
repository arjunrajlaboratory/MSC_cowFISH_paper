#Figures for paper
if(!exists('forscattersGAPDH10')) source('./data_import/calculate_all_stats_forpaper.R')
systemType <- 'Mac'
theme_set(theme_bw(base_size = 18))

donor_colors_gray <- c('grey77', 'grey51', 'grey36')
##fig1
fig1_histograms_MSC_chondrocytes <- 
  ggplot(subset(forscattersGAPDH10, experimentid %in% c('G66', 'G64')),
         aes(x = cy.RNA)) + geom_histogram(color = '#000000', fill = 'lightseagreen') + 
  facet_grid(Cell.Type ~ ., scales = 'free_y') + 
#   theme(legend.position = 'none') + 
#   scale_fill_manual(breaks = c('MSC', 'chondrocyte'),values = c('lightseagreen', 'seagreen3')) + 
#   scale_fill_hue(l = 45) +
  xlab('Aggrecan RNA') + ylab('# of Cells') 


fig1_histograms_MSC_time <- 
  ggplot(subset(forscattersGAPDH10, Purpose == 'MSCs het differentiation in gels core' & Days.CM %in% c(1,21) & Media == 'CM+' & Donor == 3),
       aes(x = cy.RNA, fill = as.factor(Days.CM))) + 
  geom_histogram() + 
  geom_histogram(color = '#000000', alpha = 0, show_guide = FALSE) + 
  facet_grid(Days.CM ~ ., scales = 'free_y') + 
  scale_fill_manual(values = c('mediumaquamarine', 'turquoise4'), name = 'Days in \nDifferentiation\nMedia') + 
#   scale_fill_hue(l = 45) +
  theme(strip.text = element_blank(), strip.background = element_blank(),
        legend.title.align = 0.5, legend.text.align = 0.5, legend.position = c(1,1), legend.justification = c(1,1)) + 
  xlab('Aggrecan RNA') + ylab('# of Cells')
  


osteopontin_correlation <- 
  ggplot(subset(forscattersGAPDH10, Purpose == 'MSCs het differentiation in gels core' & Days.CM %in% c(1,21) & Media == 'CM+' & Donor == 3),
       aes(x = cy.RNA, y = tmr.RNA, color = as.factor(Days.CM))) + 
  geom_point(size = 4) + 
#   geom_smooth(method = lm, se = FALSE, fullrange = TRUE, size = 2) +
  scale_color_manual(values = c('mediumaquamarine', 'turquoise4'), name = 'Days in \nDifferentiation\nMedia') +   
  xlab('Aggrecan RNA') + ylab('Osteopontin RNA') + theme(legend.position = 'none')

lpl_correlation <- 
  ggplot(subset(forscattersGAPDH10, Purpose == 'MSCs het differentiation in gels core' & Days.CM %in% c(1,21) & Media == 'CM+' & Donor == 3),
       aes(x = cy.RNA, y = alexa.RNA, color = as.factor(Days.CM))) + 
  geom_point(size = 4) + 
#   geom_smooth(method = lm, se = FALSE, fullrange = TRUE, size = 2, show_guide = FALSE) +
  scale_color_manual(values = c('mediumaquamarine', 'turquoise4'), name = 'Days in \nDifferentiation\nMedia') +   
  xlab('Aggrecan RNA') + ylab('LPL RNA') + theme(legend.position = c(1, 1), legend.justification = c(1,1))

GAPDH_correlation <- 
  ggplot(subset(forscattersGAPDH10, Purpose == 'MSCs het differentiation in gels core' & Days.CM %in% c(1,21) & Media == 'CM+' & Donor == 3),
         aes(x = cy.RNA, y = GAPDH.RNA, color = as.factor(Days.CM))) + 
  geom_point(size = 4) + 
#   geom_smooth(method = lm, se = FALSE, fullrange = TRUE, size = 2, show_guide = FALSE) +
  scale_color_manual(values = c('mediumaquamarine', 'turquoise4'), name = 'Days in \nDifferentiation\nMedia') +   
  xlab('Aggrecan RNA') + ylab('GAPDH RNA') + theme(legend.position = c(1, 1), legend.justification = c(1,1))

GAPDH_OPN_correlation <- 
  ggplot(subset(forscattersGAPDH10, Purpose == 'MSCs het differentiation in gels core' & Days.CM %in% c(1,21) & Media == 'CM+' & Donor == 3),
         aes(x = GAPDH.RNA, y = tmr.RNA, color = as.factor(Days.CM))) + 
  geom_point(size = 4) + 
  #   geom_smooth(method = lm, se = FALSE, fullrange = TRUE, size = 2, show_guide = FALSE) +
  scale_color_manual(values = c('mediumaquamarine', 'turquoise4'), name = 'Days in \nDifferentiation\nMedia') +   
  xlab('GAPDH RNA') + ylab('Osteopontin RNA') + theme(legend.position = c(1, 1), legend.justification = c(1,1))



donor_variability_aggrecan_jitter <-
  ggplot(subset(forscattersGAPDH10, Days.CM == 7 & Purpose == 'MSCs het differentiation in gels core' & Media == 'CM+'), 
         aes(x = as.factor(Donor), y = cy.RNA)) + 
  geom_jitter(color = 'mediumaquamarine', position = position_jitter(width = 0.1), size = 3) + 
  geom_boxplot(alpha = 0.8, fill = 'white', width = 0.4, color = 'mediumaquamarine') +
  scale_x_discrete(labels = c('A', 'B', 'C')) + xlab('Donor') + ylab('Aggrecan RNA\n(in MSCs after 7 Days in Agarose)')

donor_variability_aggrecan_dotplot <-
  ggplot(subset(forscattersGAPDH10, Days.CM == 7 & Purpose == 'MSCs het differentiation in gels core' & Media == 'CM+'), 
         aes(x = as.factor(Donor), y = cy.RNA)) + 
  geom_dotplot(binaxis = 'y',binwidth = 1, stackdir = 'center', 
               dotsize = 25, color = '#666666', fill = 'mediumaquamarine') + 
  geom_boxplot(alpha = 0.8, fill = 'white', width = 0.4, color = 'mediumaquamarine') +
  scale_x_discrete(labels = c('A', 'B', 'C')) + xlab('Donor') + ylab('Aggrecan RNA (in MSCs after 7 Days in Agarose)')

msc_mean_aggrecan <- 
  ggplot(subset(allstats, Purpose == 'MSCs het differentiation in gels core' & Days.CM %in% c(1, 4, 7, 14, 21) & Donor %in% c(3, 5, 7)), 
         aes(x = as.factor(Days.CM), y = Aggrecan.RNA.mean, fill = Donor)) + 
  geom_errorbar(aes(ymin = Aggrecan.RNA.mean - Aggrecan.RNA.mean/20, ymax = Aggrecan.RNA.mean + Aggrecan.RNA.se, y = NULL), 
                width = 0.3, position=position_dodge(.5)) + 
  geom_bar(width = 0.5, stat = 'identity', position = position_dodge()) +
  scale_fill_manual(values = donor_colors_gray, labels = c('A', 'B', 'C')) +
  facet_grid(Media ~ .) + 
  theme(legend.position = c(1,1), legend.justification = c(1, 1), strip.background = element_blank(), strip.text = element_blank()) +
  geom_bar(data = mean_mscdiff, aes(x = as.factor(Days.CM), y = Aggrecan.RNA.mean.mean), 
           fill = 'mediumseagreen', alpha = 0.6, stat = 'identity', color = '#000000', width = 0.7) +
  scale_x_discrete(limits = c('1', '4', '7', '14', '21')) +
  ylab('Mean Aggrecan RNA') + xlab('Days in Agarose')




#calculate pearson correlations for graphs
corD1<-subset(forscattersGAPDH10, Purpose == 'MSCs het differentiation in gels core' & Days.CM %in% c(1) & Media == 'CM+' & Donor == 3)
corD21<-subset(forscattersGAPDH10, Purpose == 'MSCs het differentiation in gels core' & Days.CM %in% c(21) & Media == 'CM+' & Donor == 3)

trackCorrelations<-c()
#agg-opn
trackCorrelations["Day 1 AGG-OPN"]<-cor(corD1$cy.RNA, corD1$tmr.RNA)
trackCorrelations["Day 21 AGG-OPN"]<-cor(corD21$cy.RNA, corD21$tmr.RNA)

cor.test(corD1$cy.RNA, corD1$tmr.RNA)
cor.test(corD21$cy.RNA, corD21$tmr.RNA)

#agg-lpl
trackCorrelations["Day 1 AGG-LPL"]<-cor(corD1$cy.RNA, corD1$alexa.RNA)
trackCorrelations["Day 21 AGG-LPL"]<-cor(corD21$cy.RNA, corD21$alexa.RNA)

cor.test(corD1$cy.RNA, corD1$alexa.RNA)
cor.test(corD21$cy.RNA, corD21$alexa.RNA)

#lpl-opn
trackCorrelations["Day 1 LPL-OPN"]<-cor(corD1$tmr.RNA, corD1$alexa.RNA)
trackCorrelations["Day 21 LPL-OPN"]<-cor(corD21$tmr.RNA, corD21$alexa.RNA)

cor.test(corD1$tmr.RNA, corD1$alexa.RNA)
cor.test(corD21$tmr.RNA, corD21$alexa.RNA)

R2<-trackCorrelations^2

trackCorrelations
R2
#lineage marker correlation statistics for experimentid A12 (Day 1): OPN-AGG correlation R^2=0.2409, LPL-AGG R^2=7.5802e-05, LPL-OPN R^2= 2.00479e-05; n=105
#lineage marker correlation statistics for experimentid A20 (Day 21): OPN-AGG correlation R^2=0.118, LPL-AGG R^2=.0022, LPL-OPN R^2= .045; n=79

##printing plots
plots_square <- c('lpl_correlation', 'osteopontin_correlation', 'GAPDH_correlation',
           'donor_variability_aggrecan_dotplot', 'donor_variability_aggrecan_jitter', 'GAPDH_OPN_correlation')
plots_rectangle <- c('fig1_histograms_MSC_chondrocytes', 'fig1_histograms_MSC_time', 'msc_mean_aggrecan')

directory <- paste('./Figure1/graphs/', Sys.Date(),'/', sep = '')
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
