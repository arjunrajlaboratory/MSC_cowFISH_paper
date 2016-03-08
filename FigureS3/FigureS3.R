#Figures for paper
if(!exists('forscattersGAPDH10')) source('./data_import/calculate_all_stats_forpaper.R')
#library(XLConnect)
theme_set(theme_bw(base_size = 14))
systemType <- 'Mac'

## graphing -------
sistercell_fluctuations_AGG_GAPDH <-
  ggplot(all_div, aes(x = Daughter1_cy.RNA-Daughter2_cy.RNA, y = Daughter1_GAPDH.RNA-Daughter2_GAPDH.RNA, color = TimeSinceDivision)) + geom_point(size = 4) + 
  scale_color_continuous(low = 'sienna1',  high = 'cornflowerblue', name = 'Time Since Divisions\n(hours)') +
  geom_segment(aes(x = -600, y = -600, xend = 600, yend = 600), color = 'grey', linetype = 'dashed', alpha = 0.5) +
  theme(legend.position = c(1,0), legend.justification = c(1,0)) +
  xlim(-600, 600) + ylim(-600, 600)+
  xlab('Aggrecan Difference (Sister #1-Sister #2)') + ylab('GAPDH Difference (Sister #1-Sister #2)')

cor.test((all_div$Daughter1_cy.RNA-all_div$Daughter2_cy.RNA), (all_div$Daughter1_GAPDH.RNA-all_div$Daughter2_GAPDH.RNA))

sistercell_fluctuations_AGG_OPN <-
  ggplot(all_div, aes(x = Daughter1_cy.RNA-Daughter2_cy.RNA, y = Daughter1_tmr.RNA-Daughter2_tmr.RNA, color = TimeSinceDivision)) + geom_point(size = 4) + 
  scale_color_continuous(low = 'sienna1',  high = 'cornflowerblue', name = 'Time Since Divisions\n(hours)') +
  geom_segment(aes(x = -750, y = -750, xend = 750, yend = 750), color = 'grey', linetype = 'dashed', alpha = 0.5) +
  theme(legend.position = c(1,0), legend.justification = c(1,0)) +
  #xlim(-750, 750) + ylim(-750, 750)+
  xlab('Aggrecan Difference (Sister #1-Sister #2)') + ylab('Osteopontin Difference (Sister #1-Sister #2)')  

cor.test((all_div$Daughter1_cy.RNA-all_div$Daughter2_cy.RNA), (all_div$Daughter1_tmr.RNA-all_div$Daughter2_tmr.RNA))

sistercell_GAPDH <-
ggplot(all_div, aes(x = minDaughter_GAPDH.RNA, y = maxDaughter_GAPDH.RNA, color = TimeSinceDivision)) + geom_point(size = 4) + 
  scale_color_continuous(low = 'sienna1',  high = 'cornflowerblue', name = 'Time Since Divisions\n(hours)') + xlab('Sister Cell # 1 GAPDH RNA') +
  ylab('Sister Cell # 2 GAPDH RNA') + 
  geom_segment(aes(x = 0, y = 0, xend = 900, yend = 900), color = 'grey', linetype = 'dashed', alpha = 0.5) +
  theme(legend.position = c(1,0), legend.justification = c(1,0)) + xlim(0, 925) + ylim(0, 925)

sistercell_aggperGAPDH <-
  ggplot(all_div, aes(x = minDaughter_aggpergapdh, y = maxDaughter_aggpergapdh, color = TimeSinceDivision)) + geom_point(size = 4) + 
  scale_color_continuous(low = 'sienna1',  high = 'cornflowerblue', name = 'Time Since Divisions\n(hours)') + xlab('Sister Cell # 1 AGG/GAPDH RNA') +
  ylab('Sister Cell # 2 AGG/GAPDH RNA') + 
  geom_segment(aes(x = 0, y = 0, xend = 3, yend = 3), color = 'grey', linetype = 'dashed', alpha = 0.5) +
  theme(legend.position = c(1,0), legend.justification = c(1,0)) + xlim(0, 3) + ylim(0, 3)

gapdh_difference_with_category <- 
  ggplot(all_div, aes(x = as.factor(TimeCategory), y = (maxDaughter_GAPDH.RNA-minDaughter_GAPDH.RNA))) + 
  geom_jitter(stackdir = 'center', dotsize = 4, color = 'cornflowerblue', position = position_jitter(width = .1)) + 
  ylab('Difference in \nGAPDH RNA Number \nBetween Sister Cells') + xlab('Time Since Division (in hours)') + 
  theme(legend.position = "none") + theme(text = element_text(size = 20) ) + 
  geom_boxplot(alpha = 0.8, width = 0.5, fill = 'white',color = 'royalblue4')

aggpergapdh_difference_with_category <- 
  ggplot(all_div, aes(x = as.factor(TimeCategory), y = aggpergapdh_difference)) + 
  geom_jitter(stackdir = 'center', dotsize = 4, color = 'cornflowerblue', position = position_jitter(width = .1)) + 
  ylab('Difference in \nAggrecan/GAPDH \nBetween Sister Cells') + xlab('Time Since Division (in hours)') + 
  theme(legend.position = "none") + theme(text = element_text(size = 20) ) + 
  geom_boxplot(alpha = 0.8, width = 0.5, fill = 'white',color = 'royalblue4')


##printing plots
plots_square <- c(
   'sistercell_fluctuations_AGG_GAPDH',
   'sistercell_fluctuations_AGG_OPN',
   'sistercell_GAPDH',
   'sistercell_aggperGAPDH'   
  )

plots_rectangle<-c(
  'gapdh_difference_with_category',
  'aggpergapdh_difference_with_category'
  )

  
  directory <- paste('./FigureS3/graphs/', Sys.Date(),'/', sep = '')
  
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
  

