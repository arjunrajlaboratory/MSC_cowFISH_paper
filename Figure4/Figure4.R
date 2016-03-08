#Figures for paper
if(!exists('forscattersGAPDH10')) source('./data_import/calculate_all_stats_forpaper.R')
library(XLConnect)
theme_set(theme_bw(base_size = 14))
systemType <- 'Mac'

#small colony stuff
## specific data import -------
# colonypositions<-getURL("https://docs.google.com/spreadsheets/d/1BmlsjFPeXhd7bJfhXWdPwDiHt_k7BABHQDBcUcmk0Kw/export?format=csv")
# colonypositions<-read.csv(textConnection(colonypositions),header=TRUE,sep=",")
if (systemType=="Mac"){
colonypositions <- read.csv('~/Dropbox/Cow FISHing/all_csvs/colony_positions.csv', header = TRUE)
colonyspotcounts<-read.table("~/Dropbox/Cow FISHing/all_csvs/G3.csv",header=TRUE,sep=",")
}


colonypositions$Colony.Number <- as.factor(colonypositions$Colony.Number)
levels(colonypositions$Colony.Number) <- c('A','B','C','D','E')


colonyspotcounts<-subset(colonyspotcounts, isGood==1)
colonyspotcounts<-merge(colonyspotcounts,colonypositions,'dataFile')
colonyspotcounts<-subset(colonyspotcounts, Colony.Number != 'E')
colonyspotcounts.summary <- ddply(colonyspotcounts, .(Colony.Number), summarize, 
                                  Aggrecan.RNA.mean = mean(cy.RNA), 
                                  Aggrecan.RNA.sd = sd(cy.RNA),
                                  Aggrecan.RNA.CV = sd(cy.RNA)/mean(cy.RNA),
                                  N = length(Colony.Number))

# if (systemType=="Mac"){
#   fish_spots <- read.csv('~/Dropbox/Cow FISHing/all_csvs/G69.csv')
#   fish_spots2 <- read.csv('~/Dropbox/Cow FISHing/all_csvs/G70.csv')
#   division_timing <- readWorksheetFromFile('~/Dropbox/Cow FISHing/all_csvs/division_timing_left_well.xlsx', sheet = 1)
#   timelapse_FISH_registration <- readWorksheetFromFile('~/Dropbox/Cow FISHing/all_csvs/timelapse_fish_cell_matchup.xlsx', sheet = 1)
# }
# 
# fish_spots <- cbind(fish_spots, 'Well' = 'left')
# fish_spots2 <- cbind(fish_spots2, 'Well' = 'right')
# fish_spots <- rbind(fish_spots, fish_spots2)
# rm(fish_spots2)
# dataFile <- str_split(fish_spots$ImageFile, pattern = 'cy|\\.')
# dataFile <- sapply(dataFile, '[', 2)
# fish_spots <- cbind(fish_spots, dataFile)
# fish_spots$dataFile <- as.numeric(as.character(fish_spots$dataFile))

# division_timing <- getURL('https://docs.google.com/spreadsheets/d/15tuKGYHwqy6AbON8L3ZsJ7MWedGUCJJweyQ4-uywgQI/pubhtml?gid=0&single=true&export?format=csv')
# division_timing <- read.csv(textConnection(division_timing), header = TRUE, sep = ',', na.strings = 'NA')
# division_timing <- division_timing[c('Parent_Cell', 'Division_Scan', 'Daughter1', 'Daughter2', 'Well')]
# division_timing <- melt(division_timing, id.vars = c('Parent_Cell', 'Division_Scan'))
# 
# timelapse_FISH_registration_left <- subset(timelapse_FISH_registration, Well == 'left')
# timelapse_FISH_registration <- subset(timelapse_FISH_registration, Well == 'right')
# timelapse_FISH_registration_left <- transform(timelapse_FISH_registration_left, DataFileNumber = SegmentFileNumber + 33 + round(SegmentFileNumber/32))
# timelapse_FISH_registration <- rbind(timelapse_FISH_registration, timelapse_FISH_registration_left)
# rm(timelapse_FISH_registration_left)
# 
# fish_spots_with_timelapse_number <- merge(x = fish_spots, y = timelapse_FISH_registration, by.y = c('Well' , 'DataFileNumber', 'objNum'), by.x = c('Well', 'dataFile', 'objNum'))
# fish_spots_with_div_number <- merge(x = fish_spots_with_timelapse_number, y = division_timing, by.x = 'Cell.Number', by.y = 'value')
# 
# all_div <- fish_spots_with_div_number[c('Cell.Number', 'cy.RNA', 'GAPDH.RNA', 'tmr.RNA', 'alexa.RNA', 'variable', 'Division_Scan', 'Parent_Cell', 'Well')]
# all_div <- rename(all_div, c('variable' = 'daughter_variable'))
# all_div <- melt(all_div, id.vars = c('Cell.Number', 'daughter_variable', 'Division_Scan', 'Parent_Cell', 'Well'))
# all_div <- dcast(all_div, 'Division_Scan+Parent_Cell+Well ~ daughter_variable+variable', value.var = 'value')
# all_div_left <- subset(all_div, Well == 'left')
# all_div_right <- subset(all_div, Well == 'right')
# all_div_left <- transform(all_div_left, TimeSinceDivision = ((131 - (Division_Scan))*30)/60)
# all_div_right <- transform(all_div_right, TimeSinceDivision = ((138 - (Division_Scan))*30)/60)
# all_div <- rbind(all_div_left, all_div_right)
# rm(all_div_left, all_div_right)
# all_div$maxDaughter_cy.RNA <- pmax(all_div$Daughter1_cy.RNA, all_div$Daughter2_cy.RNA)
# all_div$minDaughter_cy.RNA <- pmin(all_div$Daughter1_cy.RNA, all_div$Daughter2_cy.RNA)
# all_div <- transform(all_div, Daughter1_aggpergaph = (Daughter1_cy.RNA/Daughter1_GAPDH.RNA), 
#                      Daughter2_aggpergapdh = (Daughter2_cy.RNA/Daughter2_GAPDH.RNA))
# all_div$maxDaughter_aggpergapdh <- pmax(all_div$Daughter1_aggpergaph, all_div$Daughter2_aggpergapdh)
# all_div$minDaughter_aggpergapdh <- pmin(all_div$Daughter1_aggpergaph, all_div$Daughter2_aggpergapdh)
# all_div$maxDaughter_tmr.RNA <- pmax(all_div$Daughter1_tmr.RNA, all_div$Daughter2_tmr.RNA)
# all_div$minDaughter_tmr.RNA <- pmin(all_div$Daughter1_tmr.RNA, all_div$Daughter2_tmr.RNA)
# all_div$maxDaughter_GAPDH.RNA <- pmax(all_div$Daughter1_GAPDH.RNA, all_div$Daughter2_GAPDH.RNA)
# all_div$minDaughter_GAPDH.RNA <- pmin(all_div$Daughter1_GAPDH.RNA, all_div$Daughter2_GAPDH.RNA)
# all_div <- transform(all_div, cy_distance = abs(minDaughter_cy.RNA - maxDaughter_cy.RNA)/sqrt(2))
# all_div <- transform(all_div, cy_difference = abs(minDaughter_cy.RNA - maxDaughter_cy.RNA))
# all_div <- transform(all_div, aggpergapdh_difference = abs(minDaughter_aggpergapdh - maxDaughter_aggpergapdh))
# all_div <- transform(all_div, tmr_distance = abs(minDaughter_tmr.RNA - maxDaughter_tmr.RNA)/sqrt(2))
# all_div <- transform(all_div, GAPDH_distance = abs(minDaughter_GAPDH.RNA - maxDaughter_GAPDH.RNA)/sqrt(2), 
#                      aggpergapdh_distance = abs(minDaughter_aggpergapdh - maxDaughter_aggpergapdh)/sqrt(2),
#                      GAPDH_difference = abs(minDaughter_GAPDH.RNA - maxDaughter_GAPDH.RNA))
# all_div <- transform(all_div, mean_cy = (Daughter1_cy.RNA + Daughter2_cy.RNA)/2)
# all_div <- subset(all_div, Daughter2_cy.RNA != 'NA' & Daughter1_cy.RNA != 'NA' & Daughter1_GAPDH.RNA >= 10 & Daughter2_GAPDH.RNA >= 10)
# # all_div$TimeCategory <- ifelse(all_div$TimeSinceDivision > 12, '> 12', '0-12')
# all_div$TimeCategory<-ifelse(all_div$TimeSinceDivision < 12, '< 12',
#                       ifelse(all_div$TimeSinceDivision >= 12 & all_div$TimeSinceDivision < 24, '12-23',
#                       ifelse(all_div$TimeSinceDivision >= 24 & all_div$TimeSinceDivision < 36, '24-35',
#                       ifelse(all_div$TimeSinceDivision >= 36 & all_div$TimeSinceDivision <= 47, '36-47',
#                       ifelse(all_div$TimeSinceDivision > 47, '> 47', 'NA'
#                               )))))
# all_div$TimeCategory <- factor(all_div$TimeCategory, levels=c('< 12', '12-23', '24-35', '36-47', '> 47'))


## graphing -------
smallcolony_aggrecan_dotplot <- 
  ggplot(colonyspotcounts, aes(x = factor(Colony.Number))) + 
  geom_dotplot(binaxis = 'y',binwidth = 1, stackdir = 'center', dotsize = 4, color = '#333333', fill = 'chartreuse4', aes(y = cy.RNA)) + 
  xlab('Colony') + ylab('Aggrecan RNA') + 
  theme(legend.position = "none") + theme(text = element_text(size = 30) ) + 
  geom_boxplot(aes(y = cy.RNA), alpha = 0.8, width = 0.5, fill = 'white',color = 'chartreuse4')

smallcolony_GAPDH_dotplot <-
  ggplot(colonyspotcounts, aes(x = factor(Colony.Number))) + 
  geom_dotplot(binaxis = 'y',binwidth = 1, stackdir = 'center', dotsize = 30, color = '#333333', fill = 'chartreuse4', aes(y = GAPDH.RNA)) + 
  xlab('Colony') + ylab('GAPDH RNA') + 
  theme(legend.position = "none") + theme(text = element_text(size = 30) ) + 
  geom_boxplot(aes(y = GAPDH.RNA), alpha = 0.8, width = 0.5, fill = 'white',color = 'chartreuse4') +
  ylim(c(0, 1900))

sistercell_fluctuations_AGG_GAPDH <-
  ggplot(all_div, aes(x = Daughter1_cy.RNA-Daughter2_cy.RNA, y = Daughter1_GAPDH.RNA-Daughter2_GAPDH.RNA, color = TimeSinceDivision)) + geom_point(size = 4) + 
  scale_color_continuous(low = 'sienna1',  high = 'cornflowerblue', name = 'Time Since Divisions\n(hours)') +
  geom_segment(aes(x = -600, y = -600, xend = 600, yend = 600), color = 'grey', linetype = 'dashed', alpha = 0.5) +
  theme(legend.position = c(1,0), legend.justification = c(1,0)) +
  xlim(-600, 600) + ylim(-600, 600)+
  xlab('Aggrecan Difference (Sister #1-Sister #2)') + ylab('GAPDH Difference (Sister #1-Sister #2)')  



sistercell_fluctuations_AGG_OPN <-
  ggplot(all_div, aes(x = Daughter1_cy.RNA-Daughter2_cy.RNA, y = Daughter1_tmr.RNA-Daughter2_tmr.RNA, color = TimeSinceDivision)) + geom_point(size = 4) + 
  scale_color_continuous(low = 'sienna1',  high = 'cornflowerblue', name = 'Time Since Divisions\n(hours)') +
  geom_segment(aes(x = -600, y = -600, xend = 600, yend = 600), color = 'grey', linetype = 'dashed', alpha = 0.5) +
  theme(legend.position = c(1,0), legend.justification = c(1,0)) +
  xlim(-600, 600) + ylim(-600, 600)+
  xlab('Aggrecan Difference (Sister #1-Sister #2)') + ylab('Osteopontin Difference (Sister #1-Sister #2)')  





sistercell_aggrecan <-
  ggplot(all_div, aes(x = minDaughter_cy.RNA, y = maxDaughter_cy.RNA, color = TimeSinceDivision)) + geom_point(size = 4) + 
  scale_color_continuous(low = 'sienna1',  high = 'cornflowerblue', name = 'Time Since Divisions\n(hours)') + xlab('Sister Cell # 1 Aggrecan RNA') +
  ylab('Sister Cell # 2 Aggrecan RNA') + 
  geom_segment(aes(x = 0, y = 0, xend = 475, yend = 475), color = 'grey', linetype = 'dashed', alpha = 0.5) +
  theme(legend.position = c(1,0), legend.justification = c(1,0)) +
  xlim(0, 475) + ylim(0, 475)

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


aggpergapdh_difference_vs_minaggpergapdh<-
  ggplot(all_div, aes(y =aggpergapdh_difference, x = minDaughter_aggpergapdh, color = TimeSinceDivision)) + geom_point(size = 4) + 
  scale_color_continuous(low = 'sienna1',  high = 'cornflowerblue', name = 'Time Since Divisions\n(hours)') + xlab('Sister Cell # 1 GAPDH RNA') +
  ylab('Difference in AGG/GAPDH between sister cells') + 
  geom_smooth(method=lm) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = 'grey', linetype = 'dashed', alpha = 0.5)

maxAggvsTime<-
  ggplot(all_div, aes(x = as.factor(TimeCategory), y = maxDaughter_cy.RNA)) + 
  geom_jitter(stackdir = 'center', dotsize = 4, color = 'cornflowerblue', position = position_jitter(width = .1)) + 
  ylab('Max Number Agg') + xlab('Time Since Division (in hours)') + 
  theme(legend.position = "none") + theme(text = element_text(size = 20) ) + 
  geom_boxplot(alpha = 0.8, width = 0.5, fill = 'white',color = 'royalblue4')


maxGAPDHvsTime<-
  ggplot(all_div, aes(x = as.factor(TimeCategory), y = maxDaughter_GAPDH.RNA)) + 
  geom_jitter(stackdir = 'center', dotsize = 4, color = 'cornflowerblue', position = position_jitter(width = .1)) + 
  ylab('Max Number GAPDH') + xlab('Time Since Division (in hours)') + 
  theme(legend.position = "none") + theme(text = element_text(size = 20) ) + 
  geom_boxplot(alpha = 0.8, width = 0.5, fill = 'white',color = 'royalblue4')



# ggplot(all_div, aes(x = cy_distance)) + geom_histogram() + facet_grid(Well ~ .)


# distance_with_category <- 
# ggplot(all_div, aes(x = as.factor(TimeCategory), y = cy_distance)) + 
#   geom_dotplot(binaxis = 'y',binwidth = 1, stackdir = 'center', dotsize = 4, color = '#333333', fill = 'mediumseagreen') + 
#   ylab('Distance from y = x ') + xlab('Time Since Division (in hours)') + 
#   theme(legend.position = "none") + theme(text = element_text(size = 30) ) + 
#   geom_boxplot(alpha = 0.8, width = 0.5, fill = 'white',color = 'mediumseagreen')

agg_difference_with_category <- 
  ggplot(all_div, aes(x = as.factor(TimeCategory), y = cy_difference)) + 
  geom_jitter(stackdir = 'center', dotsize = 4, color = 'cornflowerblue', position = position_jitter(width = .1)) + 
  ylab('Difference in \nAggrecan RNA Number \nBetween Sister Cells') + xlab('Time Since Division (in hours)') + 
  theme(legend.position = "none") + theme(text = element_text(size = 20) ) + 
  geom_boxplot(alpha = 0.8, width = 0.5, fill = 'white',color = 'royalblue4')+
  ylim(0, 300)

aggpergapdh_difference_with_category <- 
  ggplot(all_div, aes(x = as.factor(TimeCategory), y = aggpergapdh_difference)) + 
  geom_jitter(stackdir = 'center', dotsize = 4, color = 'cornflowerblue', position = position_jitter(width = .1)) + 
  ylab('Difference in \nAggrecan/GAPDH \nBetween Sister Cells') + xlab('Time Since Division (in hours)') + 
  theme(legend.position = "none") + theme(text = element_text(size = 20) ) + 
  geom_boxplot(alpha = 0.8, width = 0.5, fill = 'white',color = 'royalblue4')


fractional_aggpergapdh_difference_with_category <- 
  ggplot(all_div, aes(x = as.factor(TimeCategory), y = aggpergapdh_difference/maxDaughter_aggpergapdh)) + 
  geom_jitter(stackdir = 'center', dotsize = 4, color = 'cornflowerblue', position = position_jitter(width = .1)) + 
  ylab('Fractional Difference in \nAggrecan/GAPDH \nBetween Sister Cells') + xlab('Time Since Division (in hours)') + 
  theme(legend.position = "none") + theme(text = element_text(size = 20) ) + 
  geom_boxplot(alpha = 0.8, width = 0.5, fill = 'white',color = 'royalblue4')

agg_distance_with_category <- 
  ggplot(all_div, aes(x = as.factor(TimeCategory), y = cy_distance)) + 
  geom_jitter(stackdir = 'center', dotsize = 4, color = 'cornflowerblue', position = position_jitter(width = .1)) + 
  ylab('Distance from Equal\nAggrecan Amounts') + xlab('Time Since Division (in hours)') + 
  theme(legend.position = "none") + theme(text = element_text(size = 20) ) + 
  geom_boxplot(alpha = 0.8, width = 0.5, fill = 'white',color = 'royalblue4')

gapdh_distance_with_category <-
  ggplot(all_div, aes(x = as.factor(TimeCategory), y = GAPDH_distance)) + 
  geom_jitter(stackdir = 'center', dotsize = 4, color = 'cornflowerblue', position = position_jitter(width = .1)) + 
  ylab('Distance from Equal\nGAPDH Amounts') + xlab('Time Since Division (in hours)') + 
  theme(legend.position = "none") + theme(text = element_text(size = 20) ) + 
  geom_boxplot(alpha = 0.8, width = 0.5, fill = 'white',color = 'royalblue4')

gapdh_difference_with_category <-
  ggplot(all_div, aes(x = as.factor(TimeCategory), y = GAPDH_difference)) + 
  geom_jitter(stackdir = 'center', dotsize = 4, color = 'cornflowerblue', position = position_jitter(width = .1)) + 
  ylab('Difference in GAPDH Count \nBetween Sister Cells') + xlab('Time Since Division (in hours)') + 
  theme(legend.position = "none") + theme(text = element_text(size = 20) ) + 
  geom_boxplot(alpha = 0.8, width = 0.5, fill = 'white',color = 'royalblue4')


# sistercell_GAPDH <- 
#   ggplot(all_div, aes(x = minDaughter_GAPDH.RNA, y = maxDaughter_GAPDH.RNA, color = TimeSinceDivision)) + geom_point(size = 4) + 
#   scale_color_continuous(low = 'sienna1',  high = 'cyan3', name = 'Time Since Divisions\n(hours)') + xlab('Sister Cell # 1 GAPDH RNA') +
#   ylab('Sister Cell # 2 GAPDH RNA') +
#   geom_segment(aes(x = 0, y = 0, xend = 1100, yend = 1100), color = 'grey', linetype = 'dashed', alpha = 0.5) +
#   theme(legend.position = c(1,0), legend.justification = c(1,0))
# aggrecan_distance <-
# ggplot(all_div, aes(x = TimeSinceDivision, y = cy_distance)) + geom_point(color = 'chartreuse4', size = 3) +
#   xlab('Time Since Division (hours)') + ylab('Distance From Equal (for aggrecan RNA)')
# GAPDH_distance <- 
# ggplot(all_div, aes(x = TimeSinceDivision, y = GAPDH_distance)) + geom_point(color = 'chartreuse4', size = 3) +
#   xlab('Time Since Division (hours)') + ylab('Distance From Equal (for GAPDH RNA)')




##printing plots
plots <- c('smallcolony_aggrecan_dotplot', 'smallcolony_GAPDH_dotplot',
           'sistercell_aggrecan', 'sistercell_GAPDH', 'sistercell_aggperGAPDH',
           'agg_distance_with_category', 'gapdh_distance_with_category',
           'agg_difference_with_category', 'gapdh_difference_with_category',
           'aggpergapdh_difference_with_category'         
           )



if(systemType=="Mac"){
  
  directory <- paste('./Figure4/graphs/', Sys.Date(),'/', sep = '')
  
  if(directory %in% dir() == FALSE) dir.create(directory) 
  for (i in seq(along = plots)){
    savefile <- paste(Sys.Date(),'_',as.name(plots[i]),'.pdf',sep = '')
    ggsave(filename = savefile, plot = get(plots[i]), device = pdf,
           path = directory, scale = 1, 
           width = 20, height = 20, units = "cm", dpi = 300)
  }
}


if(systemType=="Windows"){
  setwd("Figure4\\")
  directory <- paste('graphs\\', Sys.Date(),'\\', sep = "")
  if(directory %in% dir() == FALSE) dir.create(directory)
  for (i in seq(along = plots)){
    savefile <- paste(directory, Sys.Date(),'_',as.name(plots[i]),'.pdf',sep = '')
    ggsave(filename = savefile, plot = get(plots[i]), device = pdf,
           scale = 1, 
           width = 20, height = 20, units = "cm", dpi = 300)
  }
  setwd("..\\")
}


# directory <- paste('./Figure2/graphs/', Sys.Date(),'/', sep = '')
# if(directory %in% dir() == FALSE) dir.create(directory) 
# for (i in seq(along = plots)){
#   savefile <- paste(Sys.Date(),'_',as.name(plots[i]),'.pdf',sep = '')
#   ggsave(filename = savefile, plot = get(plots[i]), device = pdf,
#          path = directory, scale = 1, 
#          width = 20, height = 20, units = "cm", dpi = 300)
# }
#rm(colonypositions, colonyspotcounts, colonyspotcounts.summary, plots, savefile, smallcolony_aggrecan_dotplot, i, directory)

