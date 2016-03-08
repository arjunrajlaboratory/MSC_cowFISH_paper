if(!exists('forscatters')) source('./data_import/import_alldata_core.R')
forscattersGAPDH10 <- subset(forscatters, GAPDH.RNA >=10)
systemType <- 'Mac'
#calculating means and others statistics
calculate_all_stats <- function(df) ddply(df, .(experimentid), summarise, 
                                          N    = length(Cell.Number),
                                          Aggrecan.RNA.mean = mean(cy.RNA),
                                          Aggrecan.RNA.sd   = sd(cy.RNA),
                                          Aggrecan.RNA.se   = sd(cy.RNA) / sqrt(length(Cell.Number)),
                                          Aggrecan.RNA.CV  = sd(cy.RNA)/mean(cy.RNA),
                                          GAPDH.RNA.mean = mean(GAPDH.RNA),
                                          GAPDH.RNA.sd   = sd(GAPDH.RNA),
                                          GAPDH.RNA.se   = sd(GAPDH.RNA) / sqrt(length(Cell.Number)),
                                          GAPDH.RNA.CV  = sd(GAPDH.RNA)/mean(GAPDH.RNA),
                                          Osteopontin.RNA.mean = mean(tmr.RNA),
                                          Ostepontin.RNA.sd   = sd(tmr.RNA),
                                          Osteopontin.RNA.se   = sd(tmr.RNA) / sqrt(length(Cell.Number)),
                                          Osteopontin.RNA.CV  = sd(tmr.RNA)/mean(tmr.RNA),
                                          LPL.RNA.mean = mean(alexa.RNA),
                                          LPL.RNA.sd   = sd(alexa.RNA),
                                          LPL.RNA.se   = sd(alexa.RNA) / sqrt(length(Cell.Number)),
                                          LPL.RNA.CV  = sd(alexa.RNA)/mean(alexa.RNA),
                                          AggperGAPDH.RNA.mean = mean(cy.RNA/GAPDH.RNA),
                                          AggperGAPDH.RNA.sd   = sd(cy.RNA/GAPDH.RNA),
                                          AggperGAPDH.RNA.se   = sd(cy.RNA/GAPDH.RNA) / sqrt(length(Cell.Number)),
                                          AggperGAPDH.RNA.CV  = sd(cy.RNA/GAPDH.RNA)/mean(cy.RNA/GAPDH.RNA),
                                          AggGAPDH.RNA.cor = cor(GAPDH.RNA,cy.RNA),
                                          AGGperOPN.mean = mean(cy.RNA/(tmr.RNA+1)),
                                          AGGperOPN.median=median(cy.RNA/(tmr.RNA+1)),
                                          GAPDH.RNA.range = max(GAPDH.RNA) - min(GAPDH.RNA),
                                          cy.RNA.range = max(cy.RNA) - min(cy.RNA),
                                          AggGAPDH.RNA.corsquared = (cor(GAPDH.RNA,cy.RNA))^2,
                                          AggGAPDH.RNA.cor.p = cor(GAPDH.RNA,cy.RNA, method = c("pearson")),
                                          AggGAPDH.RNA.cor.s = cor(GAPDH.RNA,cy.RNA, method = c("spearman")),
                                          AggGAPDH.RNA.cor.pZ = fisherz(cor(GAPDH.RNA,cy.RNA, method = c("pearson"))),                                           AggGAPDH.RNA.cor.sZ = fisherz(cor(GAPDH.RNA,cy.RNA, method = c("spearman"))),
                                          AggMEANperGAPDHMEAN = mean(cy.RNA)/mean(GAPDH.RNA),
                                          AggOPN.RNA.corsquared = (cor(cy.RNA,tmr.RNA))^2,
                                          AggLPL.RNA.corsquared = (cor(cy.RNA,alexa.RNA))^2,
                                          GAPDHOPN.RNA.corsquared = (cor(GAPDH.RNA,tmr.RNA))^2,
                                          GAPDHLPL.RNA.corsquared = (cor(GAPDH.RNA,alexa.RNA))^2,
                                          OPNLPL.RNA.corsquared = (cor(tmr.RNA,alexa.RNA))^2,
                                          OPNGAPDH.RNA.cor.p = cor(GAPDH.RNA,tmr.RNA, method = c("pearson")),
                                          AggOPN.RNA.cor.p = cor(cy.RNA,tmr.RNA, method = c("pearson")),
                                          OPNGAPDH.RNA.cor.s = cor(GAPDH.RNA,tmr.RNA, method = c("spearman")),
                                          LPLGAPDH.RNA.cor.p = cor(GAPDH.RNA,alexa.RNA, method = c("pearson")),
                                          LPLGAPDH.RNA.cor.s = cor(GAPDH.RNA,alexa.RNA, method = c("spearman"))
                                          )

# calculate_lm <- function(df) dlply(df, .(experimentid), lm, formula = GAPDH.RNA ~ cy.RNA)
# lm_stats <- calculate_lm(forscattersGAPDH10)
# lm.intercept <- ldply(lm_stats, coef)[2]
# names(lm.intercept) <- c('lm.intercept')
# lm.slope <- ldply(lm_stats, coef)[3]
# names(lm.slope) <- ldply('lm.slope')
# allstats <- cbind(calculate_all_stats(forscattersGAPDH10), lm.intercept, lm.slope)
allstats <- calculate_all_stats(forscattersGAPDH10)

noStain_forscattersGAPDH10<-subset(forscattersGAPDH10, Pre.Fix.Stain=="" & Post.Fix.Stain == "")
noStain_stats<- calculate_all_stats(noStain_forscattersGAPDH10)
noStain_stats <- merge(x = noStain_stats, y = gsheetdata, by.x = 'experimentid', by.y = 'Experiment.ID')
noStain_means <- subset(noStain_stats, Purpose == 'MSCs het differentiation in gels core' & Days.CM %in% c(1,4,7,14,21))
noStain_means <- ddply(noStain_means, .(Days.CM, Media), 
                      summarize, 
                      Aggrecan.RNA.mean.mean = mean(Aggrecan.RNA.mean), 
                      AggPerGAPDH.mean.mean=mean(AggperGAPDH.RNA.mean),
                      GAPDH.RNA.mean.mean = mean(GAPDH.RNA.mean),
                      LPL.RNA.mean.mean = mean(LPL.RNA.mean),
                      Osteopontin.RNA.mean.mean = mean(Osteopontin.RNA.mean),
                      AggGAPDH.RNA.corsq.mean = mean(AggGAPDH.RNA.corsquared), 
                      AggGAPDH.RNA.corsq.sd = sd(AggGAPDH.RNA.corsquared), 
                      AggperGAPDH.RNA.mean.mean = mean(AggperGAPDH.RNA.mean),
                      AggMEANperGAPDHMEAN.mean = mean(AggMEANperGAPDHMEAN),
                      AggMEANperGAPDHMEAN.se=sd(AggMEANperGAPDHMEAN)/sqrt(3), # n=4 donors
                      Aggrecan.RNA.CV.mean = mean(Aggrecan.RNA.CV),
                      GAPDH.RNA.CV.mean = mean(GAPDH.RNA.CV),
                      Aggrecan.RNA.mean.sd = sd(Aggrecan.RNA.mean),
                      GAPDH.RNA.mean.sd = sd(GAPDH.RNA.mean),
                      AggGAPDH.RNA.cor.p.mean = mean(AggGAPDH.RNA.cor.p),
                      AggGAPDH.RNA.cor.p.se = sd(AggGAPDH.RNA.cor.p)/sqrt(3),
                      OPNGAPDH.RNA.cor.p.mean = mean(OPNGAPDH.RNA.cor.p),
                      OPNGAPDH.RNA.cor.p.se = sd(OPNGAPDH.RNA.cor.p)/sqrt(3),
                      AggOPN.RNA.cor.p.mean = mean(AggOPN.RNA.cor.p),
                      LPLGAPDH.RNA.cor.p.mean = mean(LPLGAPDH.RNA.cor.p),
                      LPLGAPDH.RNA.cor.p.se = sd(LPLGAPDH.RNA.cor.p)/sqrt(3),
                      AGGperOPN.median.mean=mean(AGGperOPN.median),
                      N = length(experimentid))


# allstats <- calculate_all_stats(forscattersGAPDH10)
allstats <- merge(x = allstats, y = gsheetdata, by.x = 'experimentid', by.y = 'Experiment.ID')
melt_allstats <- melt(allstats, id.vars = c('experimentid', 'Purpose', 
                      'Cell.Type', 'Media', 'Donor', 'Days.CM', 'Passage', 
                      'Pre.Fix.Stain', 'Post.Fix.Stain', 'Substrate'))
melt_allstats$value <- as.numeric(as.character(melt_allstats$value)) 
  
# write.table(allstats, file = '~/code/cowfish/data/allstats.csv')

# this part is for plotting the means across donors
means_dediff <- subset(allstats, Purpose == 'chondrocyte dedifferentiation core' & Passage %in% c(0,1,3,5,7,9))
means_dediff <- ddply(means_dediff, .(Passage), 
                      summarize, 
                      Aggrecan.RNA.mean.mean = mean(Aggrecan.RNA.mean), 
                      AggPerGAPDH.mean.mean=mean(AggperGAPDH.RNA.mean),
                      GAPDH.RNA.mean.mean = mean(GAPDH.RNA.mean),
                      LPL.RNA.mean.mean = mean(LPL.RNA.mean),
                      Osteopontin.RNA.mean.mean = mean(Osteopontin.RNA.mean),
                      AggGAPDH.RNA.corsq.mean = mean(AggGAPDH.RNA.corsquared), 
                      AggGAPDH.RNA.corsq.sd = sd(AggGAPDH.RNA.corsquared), 
                      AggperGAPDH.RNA.mean.mean = mean(AggperGAPDH.RNA.mean),
                      AggMEANperGAPDHMEAN.mean = mean(AggMEANperGAPDHMEAN),
                      AggMEANperGAPDHMEAN.se=sd(AggMEANperGAPDHMEAN)/sqrt(3), # n=4 donors
                      Aggrecan.RNA.CV.mean = mean(Aggrecan.RNA.CV),
                      GAPDH.RNA.CV.mean = mean(GAPDH.RNA.CV),
                      Aggrecan.RNA.mean.sd = sd(Aggrecan.RNA.mean),
                      GAPDH.RNA.mean.sd = sd(GAPDH.RNA.mean),
                      AggGAPDH.RNA.cor.p.mean = mean(AggGAPDH.RNA.cor.p),
                      AggGAPDH.RNA.cor.p.se = sd(AggGAPDH.RNA.cor.p)/sqrt(3),
                      OPNGAPDH.RNA.cor.p.mean = mean(OPNGAPDH.RNA.cor.p),
                      OPNGAPDH.RNA.cor.p.se = sd(OPNGAPDH.RNA.cor.p)/sqrt(3),
                      AggOPN.RNA.cor.p.mean = mean(AggOPN.RNA.cor.p),
                      LPLGAPDH.RNA.cor.p.mean = mean(LPLGAPDH.RNA.cor.p),
                      LPLGAPDH.RNA.cor.p.se = sd(LPLGAPDH.RNA.cor.p)/sqrt(3),
                      AGGperOPN.median.mean=mean(AGGperOPN.median),
                      N = length(experimentid))

mean_rediff <- subset(allstats, Purpose == 'chondrocyte redifferentiation core' & Substrate == 'agarose')
mean_rediff <- ddply(mean_rediff, .(Passage, Days.CM), 
                     summarize, 
                     Aggrecan.RNA.mean.mean = mean(Aggrecan.RNA.mean), 
                     AggPerGAPDH.mean.mean=mean(AggperGAPDH.RNA.mean),
                     GAPDH.RNA.mean.mean = mean(GAPDH.RNA.mean),
                     LPL.RNA.mean.mean = mean(LPL.RNA.mean),
                     Osteopontin.RNA.mean.mean = mean(Osteopontin.RNA.mean),
                     AggGAPDH.RNA.corsq.mean = mean(AggGAPDH.RNA.corsquared), 
                     AggGAPDH.RNA.corsq.sd = sd(AggGAPDH.RNA.corsquared), 
                     AggperGAPDH.RNA.mean.mean = mean(AggperGAPDH.RNA.mean),
                     AggMEANperGAPDHMEAN.mean = mean(AggMEANperGAPDHMEAN),
                     AggMEANperGAPDHMEAN.se=sd(AggMEANperGAPDHMEAN)/sqrt(3), # n=4 donors
                     Aggrecan.RNA.CV.mean = mean(Aggrecan.RNA.CV),
                     GAPDH.RNA.CV.mean = mean(GAPDH.RNA.CV),
                     Aggrecan.RNA.mean.sd = sd(Aggrecan.RNA.mean),
                     GAPDH.RNA.mean.sd = sd(GAPDH.RNA.mean),
                     AggGAPDH.RNA.cor.p.mean = mean(AggGAPDH.RNA.cor.p),
                     AggGAPDH.RNA.cor.p.se = sd(AggGAPDH.RNA.cor.p)/sqrt(3),
                     OPNGAPDH.RNA.cor.p.mean = mean(OPNGAPDH.RNA.cor.p),
                     OPNGAPDH.RNA.cor.p.se = sd(OPNGAPDH.RNA.cor.p)/sqrt(3),
                     AggOPN.RNA.cor.p.mean = mean(AggOPN.RNA.cor.p),
                     LPLGAPDH.RNA.cor.p.mean = mean(LPLGAPDH.RNA.cor.p),
                     LPLGAPDH.RNA.cor.p.se = sd(LPLGAPDH.RNA.cor.p)/sqrt(3),
                     AGGperOPN.median.mean=mean(AGGperOPN.median),
                     N = length(experimentid))

mean_mscdiff <- subset(allstats, Purpose == 'MSCs het differentiation in gels core' & Days.CM %in% c(1, 4, 7, 14, 21))
mean_mscdiff <- ddply(mean_mscdiff, .(Media, Days.CM), 
                      summarize, 
                      Aggrecan.RNA.mean.mean = mean(Aggrecan.RNA.mean), 
                      AggPerGAPDH.mean.mean=mean(AggperGAPDH.RNA.mean),
                      GAPDH.RNA.mean.mean = mean(GAPDH.RNA.mean),
                      LPL.RNA.mean.mean = mean(LPL.RNA.mean),
                      Osteopontin.RNA.mean.mean = mean(Osteopontin.RNA.mean),
                      AggGAPDH.RNA.corsq.mean = mean(AggGAPDH.RNA.corsquared), 
                      AggGAPDH.RNA.corsq.sd = sd(AggGAPDH.RNA.corsquared), 
                      AggperGAPDH.RNA.mean.mean = mean(AggperGAPDH.RNA.mean),
                      AggMEANperGAPDHMEAN.mean = mean(AggMEANperGAPDHMEAN),
                      AggMEANperGAPDHMEAN.se=sd(AggMEANperGAPDHMEAN)/sqrt(3), # n=4 donors
                      Aggrecan.RNA.CV.mean = mean(Aggrecan.RNA.CV),
                      GAPDH.RNA.CV.mean = mean(GAPDH.RNA.CV),
                      Aggrecan.RNA.mean.sd = sd(Aggrecan.RNA.mean),
                      GAPDH.RNA.mean.sd = sd(GAPDH.RNA.mean),
                      AggGAPDH.RNA.cor.p.mean = mean(AggGAPDH.RNA.cor.p),
                      AggGAPDH.RNA.cor.p.se = sd(AggGAPDH.RNA.cor.p)/sqrt(3),
                      OPNGAPDH.RNA.cor.p.mean = mean(OPNGAPDH.RNA.cor.p),
                      OPNGAPDH.RNA.cor.p.se = sd(OPNGAPDH.RNA.cor.p)/sqrt(3),
                      AggOPN.RNA.cor.p.mean = mean(AggOPN.RNA.cor.p),
                      LPLGAPDH.RNA.cor.p.mean = mean(LPLGAPDH.RNA.cor.p),
                      LPLGAPDH.RNA.cor.p.se = sd(LPLGAPDH.RNA.cor.p)/sqrt(3),
                      AGGperOPN.median.mean=mean(AGGperOPN.median),
                      N = length(experimentid))

mean_mscdiff_nostain <- subset(allstats, Purpose == 'MSCs het differentiation in gels core' & Days.CM %in% c(1, 4, 7, 14, 21) & Pre.Fix.Stain %in% c(''))
mean_mscdiff_nostain <- ddply(mean_mscdiff_nostain, .(Media, Days.CM), 
                              summarize, 
                              Aggrecan.RNA.mean.mean = mean(Aggrecan.RNA.mean), 
                              AggPerGAPDH.mean.mean=mean(AggperGAPDH.RNA.mean),
                              GAPDH.RNA.mean.mean = mean(GAPDH.RNA.mean),
                              LPL.RNA.mean.mean = mean(LPL.RNA.mean),
                              Osteopontin.RNA.mean.mean = mean(Osteopontin.RNA.mean),
                              AggGAPDH.RNA.corsq.mean = mean(AggGAPDH.RNA.corsquared), 
                              AggGAPDH.RNA.corsq.sd = sd(AggGAPDH.RNA.corsquared), 
                              AggperGAPDH.RNA.mean.mean = mean(AggperGAPDH.RNA.mean),
                              AggMEANperGAPDHMEAN.mean = mean(AggMEANperGAPDHMEAN),
                              AggMEANperGAPDHMEAN.se=sd(AggMEANperGAPDHMEAN)/sqrt(3), # n=4 donors
                              Aggrecan.RNA.CV.mean = mean(Aggrecan.RNA.CV),
                              GAPDH.RNA.CV.mean = mean(GAPDH.RNA.CV),
                              Aggrecan.RNA.mean.sd = sd(Aggrecan.RNA.mean),
                              GAPDH.RNA.mean.sd = sd(GAPDH.RNA.mean),
                              AggGAPDH.RNA.cor.p.mean = mean(AggGAPDH.RNA.cor.p),
                              AggGAPDH.RNA.cor.p.se = sd(AggGAPDH.RNA.cor.p)/sqrt(3),
                              OPNGAPDH.RNA.cor.p.mean = mean(OPNGAPDH.RNA.cor.p),
                              OPNGAPDH.RNA.cor.p.se = sd(OPNGAPDH.RNA.cor.p)/sqrt(3),
                              AggOPN.RNA.cor.p.mean = mean(AggOPN.RNA.cor.p),
                              LPLGAPDH.RNA.cor.p.mean = mean(LPLGAPDH.RNA.cor.p),
                              LPLGAPDH.RNA.cor.p.se = sd(LPLGAPDH.RNA.cor.p)/sqrt(3),
                              AGGperOPN.median.mean=mean(AGGperOPN.median),
                              N = length(experimentid))
melt_mean_mscdiff <- melt(mean_mscdiff_nostain, id.vars = c('Media', 'Days.CM'))
melt_mean_mscdiff$value <- as.numeric(as.character(melt_mean_mscdiff$value)) 


# rm(lm.intercept, lm.slope, alldata, lm_stats, calculate_all_stats, calculate_lm)

#single cell scoring of aggrecan IF signal
singleCellScore$Donor<-as.factor(singleCellScore$Donor)
singleCellScore<- subset(singleCellScore, isGood == TRUE & GAPDH.RNA >= 10 & Score!='NA')
singleCellScore$plotScore<- ifelse(singleCellScore$Score > 3, 'Extracellular Staining', ' No Extracellular Staining')

singleCellScore_stats<- ddply(singleCellScore, .(Donor, plotScore), summarise,
                                           N = length(GAPDH.RNA),
                                           AggGAPDH.RNA.cor.p = cor(GAPDH.RNA,cy.RNA, method = c("pearson")),
                                           AggGAPDH.RNA.cor.s = cor(GAPDH.RNA,cy.RNA, method = c("spearman")),
                                           OPNGAPDH.RNA.cor.p = cor(GAPDH.RNA,tmr.RNA, method = c("pearson")),
                                           OPNGAPDH.RNA.cor.s = cor(GAPDH.RNA,tmr.RNA, method = c("spearman")),
                                           LPLGAPDH.RNA.cor.p = cor(GAPDH.RNA,alexa.RNA, method = c("pearson")),
                                           LPLGAPDH.RNA.cor.s = cor(GAPDH.RNA,alexa.RNA, method = c("spearman")),
                                           Aggrecan.RNA.mean = mean(cy.RNA),
                                           Aggrecan.RNA.sd   = sd(cy.RNA),
                                           Aggrecan.RNA.se   = sd(cy.RNA) / sqrt(length(cy.RNA)),
                                           Aggrecan.RNA.CV  = sd(cy.RNA)/mean(cy.RNA),
                                           GAPDH.RNA.mean = mean(GAPDH.RNA),
                                           GAPDH.RNA.sd   = sd(GAPDH.RNA),
                                           GAPDH.RNA.se   = sd(GAPDH.RNA) / sqrt(length(GAPDH.RNA)),
                                           GAPDH.RNA.CV  = sd(GAPDH.RNA)/mean(GAPDH.RNA),
                                           Osteopontin.RNA.mean = mean(tmr.RNA),
                                           Ostepontin.RNA.sd   = sd(tmr.RNA),
                                           Osteopontin.RNA.se   = sd(tmr.RNA) / sqrt(length(tmr.RNA)),
                                           Osteopontin.RNA.CV  = sd(tmr.RNA)/mean(tmr.RNA),
                                           LPL.RNA.mean = mean(alexa.RNA),
                                           LPL.RNA.sd   = sd(alexa.RNA),
                                           LPL.RNA.se   = sd(alexa.RNA) / sqrt(length(alexa.RNA)),
                                           LPL.RNA.CV  = sd(alexa.RNA)/mean(alexa.RNA),
                                           AggperGAPDH.RNA.mean = mean(cy.RNA/GAPDH.RNA),
                                           AggperGAPDH.RNA.sd   = sd(cy.RNA/GAPDH.RNA),
                                           AggperGAPDH.RNA.se   = sd(cy.RNA/GAPDH.RNA) / sqrt(length(cy.RNA)),
                                           AggperGAPDH.RNA.CV  = sd(cy.RNA/GAPDH.RNA)/mean(cy.RNA/GAPDH.RNA),
                                           AGGperOPN = mean(cy.RNA/(tmr.RNA+1)),
                                           AGGperOPNinf = mean(cy.RNA/(tmr.RNA)),
                                           AGGperOPN.se=sd(cy.RNA/(tmr.RNA+1))/sqrt(length(cy.RNA)),
                                           AGGperOPN.median=median(cy.RNA/(tmr.RNA+1)),
                                           OPNperGAPDH=mean(tmr.RNA/GAPDH.RNA)
)
singleCellScore_stats_AGGperOPN_removeinf <- ddply(subset(singleCellScore, tmr.RNA > 0), .(plotScore), summarize, 
                                                   AGGperOPN = mean(cy.RNA/tmr.RNA))


means_singleCellScore <-ddply(singleCellScore_stats, .(plotScore), summarise, 
                              AggGAPDH.cor.p.mean=mean(AggGAPDH.RNA.cor.p),
                              OPNGAPDH.cor.p.mean=mean(OPNGAPDH.RNA.cor.p),
                              LPLGAPDH.cor.p.mean=mean(LPLGAPDH.RNA.cor.p),
                              AggGAPDH.cor.s.mean=mean(AggGAPDH.RNA.cor.s),
                              OPNGAPDH.cor.s.mean=mean(OPNGAPDH.RNA.cor.s),
                              LPLGAPDH.cor.s.mean=mean(LPLGAPDH.RNA.cor.s),
                              Aggrecan.RNA.mean.mean = mean(Aggrecan.RNA.mean), 
                              AggPerGAPDH.mean.mean=mean(AggperGAPDH.RNA.mean),
                              GAPDH.RNA.mean.mean = mean(GAPDH.RNA.mean),
                              AggperGAPDH.RNA.mean.mean = mean(AggperGAPDH.RNA.mean),
                              Aggrecan.RNA.CV.mean = mean(Aggrecan.RNA.CV),
                              GAPDH.RNA.CV.mean = mean(GAPDH.RNA.CV),
                              Osteopontin.RNA.mean.mean = mean(Osteopontin.RNA.mean),
                              LPL.RNA.mean.mean = mean(LPL.RNA.mean),
                              AGGperOPN.mean=mean(AGGperOPN),
                              AGGperOPN.median.mean=mean(AGGperOPN.median),
                              OPNperGAPDH.mean=mean(OPNperGAPDH)                              
                              
) 



# sister cell tracking

fish_spots <- cbind(fish_spots, 'Well' = 'left')
fish_spots2 <- cbind(fish_spots2, 'Well' = 'right')
fish_spots <- rbind(fish_spots, fish_spots2)
#rm(fish_spots2)
dataFile <- str_split(fish_spots$ImageFile, pattern = 'cy|\\.')
dataFile <- sapply(dataFile, '[', 2)
fish_spots <- cbind(fish_spots, dataFile)
fish_spots$dataFile <- as.numeric(as.character(fish_spots$dataFile))

division_timing <- division_timing[c('Parent_Cell', 'Division_Scan', 'Daughter1', 'Daughter2', 'Well')]
division_timing <- melt(division_timing, id.vars = c('Parent_Cell', 'Division_Scan'))

timelapse_FISH_registration_left <- subset(timelapse_FISH_registration, Well == 'left')
timelapse_FISH_registration <- subset(timelapse_FISH_registration, Well == 'right')
timelapse_FISH_registration_left <- transform(timelapse_FISH_registration_left, DataFileNumber = SegmentFileNumber + 33 + round(SegmentFileNumber/32))
timelapse_FISH_registration <- rbind(timelapse_FISH_registration, timelapse_FISH_registration_left)
rm(timelapse_FISH_registration_left)

fish_spots_with_timelapse_number <- merge(x = fish_spots, y = timelapse_FISH_registration, by.y = c('Well' , 'DataFileNumber', 'objNum'), by.x = c('Well', 'dataFile', 'objNum'))
fish_spots_with_div_number <- merge(x = fish_spots_with_timelapse_number, y = division_timing, by.x = 'Cell.Number', by.y = 'value')

all_div <- fish_spots_with_div_number[c('Cell.Number', 'cy.RNA', 'GAPDH.RNA', 'tmr.RNA', 'alexa.RNA', 'variable', 'Division_Scan', 'Parent_Cell', 'Well')]
all_div <- rename(all_div, daughter_variable = variable)
all_div <- melt(all_div, id.vars = c('Cell.Number', 'daughter_variable', 'Division_Scan', 'Parent_Cell', 'Well'))
all_div <- dcast(all_div, 'Division_Scan+Parent_Cell+Well ~ daughter_variable+variable', value.var = 'value')
all_div_left <- subset(all_div, Well == 'left')
all_div_right <- subset(all_div, Well == 'right')
all_div_left <- transform(all_div_left, TimeSinceDivision = ((131 - (Division_Scan))*30)/60)
all_div_right <- transform(all_div_right, TimeSinceDivision = ((138 - (Division_Scan))*30)/60)
all_div <- rbind(all_div_left, all_div_right)
rm(all_div_left, all_div_right)
all_div$maxDaughter_cy.RNA <- pmax(all_div$Daughter1_cy.RNA, all_div$Daughter2_cy.RNA)
all_div$minDaughter_cy.RNA <- pmin(all_div$Daughter1_cy.RNA, all_div$Daughter2_cy.RNA)
all_div <- transform(all_div, Daughter1_aggpergaph = (Daughter1_cy.RNA/Daughter1_GAPDH.RNA), 
                     Daughter2_aggpergapdh = (Daughter2_cy.RNA/Daughter2_GAPDH.RNA))
all_div$maxDaughter_aggpergapdh <- pmax(all_div$Daughter1_aggpergaph, all_div$Daughter2_aggpergapdh)
all_div$minDaughter_aggpergapdh <- pmin(all_div$Daughter1_aggpergaph, all_div$Daughter2_aggpergapdh)
all_div$maxDaughter_tmr.RNA <- pmax(all_div$Daughter1_tmr.RNA, all_div$Daughter2_tmr.RNA)
all_div$minDaughter_tmr.RNA <- pmin(all_div$Daughter1_tmr.RNA, all_div$Daughter2_tmr.RNA)
all_div$maxDaughter_GAPDH.RNA <- pmax(all_div$Daughter1_GAPDH.RNA, all_div$Daughter2_GAPDH.RNA)
all_div$minDaughter_GAPDH.RNA <- pmin(all_div$Daughter1_GAPDH.RNA, all_div$Daughter2_GAPDH.RNA)
all_div <- transform(all_div, cy_distance = abs(minDaughter_cy.RNA - maxDaughter_cy.RNA)/sqrt(2))
all_div <- transform(all_div, cy_difference = abs(minDaughter_cy.RNA - maxDaughter_cy.RNA))
all_div <- transform(all_div, aggpergapdh_difference = abs(minDaughter_aggpergapdh - maxDaughter_aggpergapdh))

all_div <- transform(all_div, tmr_distance = abs(minDaughter_tmr.RNA - maxDaughter_tmr.RNA)/sqrt(2))
all_div <- transform(all_div, GAPDH_distance = abs(minDaughter_GAPDH.RNA - maxDaughter_GAPDH.RNA)/sqrt(2), 
                     aggpergapdh_distance = abs(minDaughter_aggpergapdh - maxDaughter_aggpergapdh)/sqrt(2),
                     GAPDH_difference = abs(minDaughter_GAPDH.RNA - maxDaughter_GAPDH.RNA))
all_div <- transform(all_div, mean_cy = (Daughter1_cy.RNA + Daughter2_cy.RNA)/2)
all_div <- subset(all_div, Daughter2_cy.RNA != 'NA' & Daughter1_cy.RNA != 'NA' & Daughter1_GAPDH.RNA >= 10 & Daughter2_GAPDH.RNA >= 10)
# all_div$TimeCategory <- ifelse(all_div$TimeSinceDivision > 12, '> 12', '0-12')
all_div$TimeCategory<-ifelse(all_div$TimeSinceDivision < 12, '< 12',
                             ifelse(all_div$TimeSinceDivision >= 12 & all_div$TimeSinceDivision < 24, '12-23',
                                    ifelse(all_div$TimeSinceDivision >= 24 & all_div$TimeSinceDivision < 36, '24-35',
                                           ifelse(all_div$TimeSinceDivision >= 36 & all_div$TimeSinceDivision <= 47, '36-47',
                                                  ifelse(all_div$TimeSinceDivision > 47, '> 47', 'NA'
                                                  )))))
all_div$TimeCategory <- factor(all_div$TimeCategory, levels=c('< 12', '12-23', '24-35', '36-47', '> 47'))
