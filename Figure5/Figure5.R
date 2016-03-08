#Figures for paper
theme_set(theme_bw(base_size = 32))
if(!exists('forscattersGAPDH10')) source('./data_import/calculate_all_stats_forpaper.R')

#set color assignments for overlaid bars (color by passage)
p0f<-"#0099cc"
p1f<-"#CCCCCC"
p3f<-p1f
p5f<-"#660099"
p7f<-p1f
p9f<-p1f

p0c<-"#000000" #"#005C7A"
p1c<-"#000000"
p3c<-p1c
p5c<-"#000000" #"#993D3D"
p7c<-p1c
p9c<-p1c

donA<-"#D6D6D6"
donB<-"#999999"
donC<-"#373737"
donD<-"#212121"

# AGG means
dediff_aggrecan_means<-
ggplot(subset(allstats, Purpose == 'chondrocyte dedifferentiation core' & Passage %in% c(0,1,3,5,7,9)), 
       aes(x = Passage, y = Aggrecan.RNA.mean, fill = Donor)) + 
  scale_fill_manual(values=c(donA,donB,donC,donD)) +
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge(), color = '#000000') + 
  geom_errorbar(aes(ymin = Aggrecan.RNA.mean - .001, ymax = Aggrecan.RNA.mean + Aggrecan.RNA.se), 
                width = 0.3, position=position_dodge(.75)) + 
  geom_bar(data = means_dediff, aes(x = as.factor(Passage), y = Aggrecan.RNA.mean.mean), 
           fill = c(p0f,p1f,p3f,p5f,p7f,p9f), color = c(p0c,p1c,p3c,p5c,p7c,p9c) , alpha = 0.6, stat = 'identity', color = '#000000') +
  theme(legend.position = 'none') +
  ylab('Mean Aggrecan RNA per Cell') + xlab('Passage')


# AGGperOPN median
  ggplot(subset(allstats, Purpose == 'chondrocyte dedifferentiation core' & Passage %in% c(0,1,3,5,7,9)), 
         aes(x = Passage, y = AGGperOPN.median, fill = Donor)) + 
  scale_fill_manual(values=c(donA,donB,donC,donD)) +
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge(), color = '#000000')+
  geom_bar(data = means_dediff, aes(x = as.factor(Passage), y = AGGperOPN.median.mean), 
           fill = c(p0f,p1f,p3f,p5f,p7f,p9f), color = c(p0c,p1c,p3c,p5c,p7c,p9c) , alpha = 0.6, stat = 'identity', color = '#000000') 


#GAPDH means
dediff_gapdh_means<-
ggplot(subset(allstats, Purpose == 'chondrocyte dedifferentiation core' & Passage %in% c(0,1,3,5,7,9)), 
       aes(x = Passage, y = GAPDH.RNA.mean, fill = Donor)) + 
  scale_fill_manual(values=c(donA,donB,donC,donD)) +
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge(), color = '#000000') + 
  geom_errorbar(aes(ymin = GAPDH.RNA.mean - .001, ymax = GAPDH.RNA.mean + GAPDH.RNA.se), 
                width = 0.3, position=position_dodge(.75)) + 
  geom_bar(data = means_dediff, aes(x = as.factor(Passage), y = GAPDH.RNA.mean.mean), 
           fill = c(p0f,p1f,p3f,p5f,p7f,p9f), color = c(p0c,p1c,p3c,p5c,p7c,p9c) , alpha = 0.6, stat = 'identity', color = '#000000') +
  theme(legend.position = 'none') +
  ylab('Mean GAPDH RNA per Cell') + xlab('Passage')

#Aggrecan per GAPDH
dediff_aggrecanPerGAPDH_pooled<- #faux PCR 
ggplot( data = means_dediff, aes(x = as.factor(Passage), y = AggMEANperGAPDHMEAN.mean)) +
  scale_fill_manual(values=c(donA,donB,donC,donD)) +
  geom_bar(data=subset(allstats, Purpose == 'chondrocyte dedifferentiation core' & Passage %in% c(0,1,3,5,7,9)), 
           aes(x = Passage, y = AggMEANperGAPDHMEAN, fill=Donor), width = 0.75, stat = 'identity', position = position_dodge()) + 
  geom_bar( 
           fill = c(p0f,p1f,p3f,p5f,p7f,p9f), color = c(p0c,p1c,p3c,p5c,p7c,p9c), alpha = 0.6, stat = 'identity', size=1.1) +
  geom_errorbar( aes(x = as.factor(Passage), ymin = AggMEANperGAPDHMEAN.mean, ymax = AggMEANperGAPDHMEAN.mean +AggMEANperGAPDHMEAN.se), 
                width = 0.3, position=position_dodge(.75), size = 1.1) + 
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  ylab('Pooled Aggrecan per GAPDH') + xlab('Passage')


dediff_aggrecanPerGAPDH_pooled_nolegend<- #faux PCR 
  ggplot( data = means_dediff, aes(x = as.factor(Passage), y = AggMEANperGAPDHMEAN.mean)) +
  scale_fill_manual(values=c(donA,donB,donC,donD)) +
  geom_bar(data=subset(allstats, Purpose == 'chondrocyte dedifferentiation core' & Passage %in% c(0,1,3,5,7,9)), 
           aes(x = Passage, y = AggMEANperGAPDHMEAN, fill=Donor), width = 0.75, stat = 'identity', position = position_dodge()) + 
  geom_bar( 
    fill = c(p0f,p1f,p3f,p5f,p7f,p9f), color = c(p0c,p1c,p3c,p5c,p7c,p9c), alpha = 0.6, stat = 'identity', size=1.1) +
  geom_errorbar( aes(x = as.factor(Passage), ymin = AggMEANperGAPDHMEAN.mean, ymax = AggMEANperGAPDHMEAN.mean +AggMEANperGAPDHMEAN.se), 
                 width = 0.3, position=position_dodge(.75), size = 1.1) + 
  theme(legend.position = 'none') +
  ylab('Pooled Aggrecan per GAPDH') + xlab('Passage')

dediff_aggrecanPerGAPDH_perCell<-
ggplot(subset(allstats, Purpose == 'chondrocyte dedifferentiation core' & Passage %in% c(0,1,3,5,7,9)), 
       aes(x = Passage, y = AggperGAPDH.RNA.mean, fill=Donor) ) +
  scale_fill_manual(values=c(donA,donB,donC,donD)) +
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge(), color ="#000000") + 
  geom_errorbar(aes(ymin = AggperGAPDH.RNA.mean, ymax = AggperGAPDH.RNA.mean + AggperGAPDH.RNA.se), 
                width = 0.3, position=position_dodge(.75)) + 
  geom_bar(data = means_dediff, aes(x = as.factor(Passage), y = AggperGAPDH.RNA.mean.mean), 
           fill = c(p0f,p1f,p3f,p5f,p7f,p9f), color = c(p0c,p1c,p3c,p5c,p7c,p9c), alpha = 0.6, stat = 'identity') +
  ylab('Mean AGG/GAPDH per Cell') + xlab('Passage')



#chondrocyte diameter
diam_all<-diam
diam<-subset(diam, passage %in% c("0","5") & donor %in% c("1","3"))
diam$passage<-as.factor(diam$passage)
diam$donor<-as.factor(diam$donor)
diam$volume<- (4/3*3.14*diam$radiusMicrons^3) # assumes sphere

diam_means <- ddply(diam, c("passage", "donor"), summarise,
                    rad_mean=mean(radiusMicrons), rad_sem=sd(radiusMicrons)/sqrt(length(radiusMicrons)),
                    volume_mean=mean(volume), volume_sem=sd(volume)/sqrt(length(volume)))
diam_means.mean <- ddply(diam_means, c("passage"), summarise,
                         rad_mean.mean=mean(rad_mean),
                         volume_mean.mean=mean(volume_mean))

chond_volume<-
ggplot(data=diam_means, aes(x = passage, y = volume_mean, fill = donor))+
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge()) + 
  geom_errorbar(data=diam_means, aes(ymin = volume_mean - rad_sem, ymax = volume_mean + rad_sem), 
                width = 0.3, position=position_dodge(.75)) + 
  scale_fill_manual(values=c(donA, donC)) + 
  geom_bar(data = diam_means.mean, aes(x = passage, y = volume_mean.mean),
           fill=c(p0f, p5f), color=c(p0c, p5c), alpha = 0.6, stat = 'identity') +
  ylab('Suspended Cell Volume (um^3)') + xlab('Passage')


## figure 4 - chondrocyte redifferentiation --------
chond_rediff_meanGAPDH <- 
  ggplot(subset(allstats, Purpose == 'chondrocyte redifferentiation core' & Substrate=="agarose"), 
       aes(x = interaction(Days.CM, Passage), y = GAPDH.RNA.mean, fill = Donor)) + 
  scale_x_discrete(limits = c('1.0', '14.0', '1.5', '14.5'), labels = c('Day 1\nP0', 'Day 14\nP0', 'Day 1\nP5', 'Day 14\nP5')) +
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE) +
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge(), color = '#000000', show_guide = FALSE, alpha = 0) + 
#   facet_grid(Passage ~ Days.CM) + 
  scale_fill_manual(breaks = c('1B.0', '2.0', '3.0', '1B.5', '2.5', '3.5'), 
                    values=c(donA, donB, donC)) +
  geom_errorbar(aes(ymin = GAPDH.RNA.mean - 0.001, ymax = GAPDH.RNA.mean + GAPDH.RNA.se), 
                width = 0.3, position=position_dodge(.75)) + 
  geom_bar(data = mean_rediff, aes(x = interaction(Days.CM, Passage), y = GAPDH.RNA.mean.mean), 
           fill = c(p0f, p0f, p5f, p5f), alpha = 0.6, stat = 'identity', color = '#000000') +
  ylab('Mean GAPDH RNA per Cell') + xlab(NULL)

#AGGperOPNmedian
  ggplot(subset(allstats, Purpose == 'chondrocyte redifferentiation core' & Substrate=="agarose"), 
         aes(x = interaction(Days.CM, Passage), y = AGGperOPN.median, fill = Donor)) + 
  scale_x_discrete(limits = c('1.0', '14.0', '1.5', '14.5'), labels = c('Day 1\nP0', 'Day 14\nP0', 'Day 1\nP5', 'Day 14\nP5')) +
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge()) +
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge(), color = '#000000', show_guide = FALSE, alpha = 0) + 
  #   facet_grid(Passage ~ Days.CM) + 
  scale_fill_manual(breaks = c('1B.0', '2.0', '3.0', '1B.5', '2.5', '3.5'), 
                    values=c(donA, donB, donC)) + 
  geom_bar(data = mean_rediff, aes(x = interaction(Days.CM, Passage), y = AGGperOPN.median.mean), 
           fill = c(p0f, p0f, p5f, p5f), alpha = 0.6, stat = 'identity', color = '#000000') 



#Aggrecan
chond_rediff_meanAgg <- 
  ggplot(subset(allstats, Purpose == 'chondrocyte redifferentiation core' & Substrate=="agarose"), 
       aes(x = interaction(Days.CM, Passage), y = Aggrecan.RNA.mean, fill = Donor)) + 
  scale_x_discrete(limits = c('1.0', '14.0', '1.5', '14.5'), labels = c('Day 1\nP0', 'Day 14\nP0', 'Day 1\nP5', 'Day 14\nP5')) +
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge()) +
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge(), color = '#000000', show_guide = FALSE, alpha = 0) + 
  #   facet_grid(Passage ~ Days.CM) + 
  scale_fill_manual(breaks = c('1B.0', '2.0', '3.0', '1B.5', '2.5', '3.5'), 
                    values=c(donA, donB, donC)) +
  geom_errorbar(aes(ymin = Aggrecan.RNA.mean - 0.001, ymax = Aggrecan.RNA.mean + Aggrecan.RNA.se), 
                width = 0.3, position=position_dodge(.75)) + 
  geom_bar(data = mean_rediff, aes(x = interaction(Days.CM, Passage), y = Aggrecan.RNA.mean.mean), 
           fill = c(p0f, p0f, p5f, p5f), alpha = 0.6, stat = 'identity', color = '#000000') +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  ylab('Mean Aggrecan RNA per Cell') + xlab(NULL)

chond_rediff_AggPerGaPDHperCell <- 
  ggplot(subset(allstats, Purpose == 'chondrocyte redifferentiation core' & Substrate=="agarose"), 
         aes(x = interaction(Days.CM, Passage), y = AggperGAPDH.RNA.mean, fill = Donor)) + 
  scale_x_discrete(limits = c('1.0', '14.0', '1.5', '14.5'), labels = c('Day 1\nP0', 'Day 14\nP0', 'Day 1\nP5', 'Day 14\nP5')) +
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge()) +
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge(), color = '#000000', show_guide = FALSE, alpha = 0) + 
  #   facet_grid(Passage ~ Days.CM) + 
  scale_fill_manual(breaks = c('1B.0', '2.0', '3.0', '1B.5', '2.5', '3.5'), 
                    values=c(donA, donB, donC)) +
  geom_errorbar(aes(ymin = AggperGAPDH.RNA.mean - 0.001, ymax = AggperGAPDH.RNA.mean + AggperGAPDH.RNA.se), 
                width = 0.3, position=position_dodge(.75)) + 
  geom_bar(data = mean_rediff, aes(x = interaction(Days.CM, Passage), y = AggperGAPDH.RNA.mean.mean), 
           fill = c(p0f, p0f, p5f, p5f), alpha = 0.6, stat = 'identity', color = '#000000') +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  ylab('Aggrecan/GAPDH per Cell') + xlab(NULL)

chond_rediff_AggPerGaPDHpooled <- 
  ggplot(subset(allstats, Purpose == 'chondrocyte redifferentiation core' & Substrate=="agarose"), 
         aes(x = interaction(Days.CM, Passage))) + 
  scale_x_discrete(limits = c('1.0', '14.0', '1.5', '14.5'), labels = c('Day 1\nP0', 'Day 14\nP0', 'Day 1\nP5', 'Day 14\nP5')) +
  geom_bar(aes(y= AggMEANperGAPDHMEAN, fill = Donor), width = 0.75, stat = 'identity', position = position_dodge()) +
  geom_bar(aes(y= AggMEANperGAPDHMEAN, fill = Donor), width = 0.75, stat = 'identity', position = position_dodge(), color = '#000000', show_guide = FALSE, alpha = 0) + 
  scale_fill_manual(breaks = c('1B.0', '2.0', '3.0', '1B.5', '2.5', '3.5'), 
                    values=c(donA, donB, donC)) +
  geom_errorbar(data = mean_rediff, aes(ymin = AggMEANperGAPDHMEAN.mean - 0.001, ymax = AggMEANperGAPDHMEAN.mean + AggMEANperGAPDHMEAN.se), 
                width = 0.3, position=position_dodge(.75)) + 
  geom_bar(data = mean_rediff, aes(x = interaction(Days.CM, Passage), y =  AggMEANperGAPDHMEAN.mean), 
           fill = c(p0f, p0f, p5f, p5f), alpha = 0.6, stat = 'identity', color = '#000000') +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  ylab('Aggrecan/GAPDH pooled') + xlab(NULL)



#osteopontin
ggplot(subset(allstats, Purpose == 'chondrocyte redifferentiation core' & Substrate=="agarose"), 
       aes(x = interaction(Days.CM, Passage), y = Osteopontin.RNA.mean, fill = Donor)) + 
  scale_x_discrete(limits = c('1.0', '14.0', '1.5', '14.5'), labels = c('Day 1\nP0', 'Day 14\nP0', 'Day 1\nP5', 'Day 14\nP5')) +
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge()) +
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge(), color = '#000000', show_guide = FALSE, alpha = 0) + 
  #   facet_grid(Passage ~ Days.CM) + 
  scale_fill_manual(breaks = c('1B.0', '2.0', '3.0', '1B.5', '2.5', '3.5'), 
                    values=c(donA, donB, donC)) +
  geom_errorbar(aes(ymin = Osteopontin.RNA.mean - 0.001, ymax = Osteopontin.RNA.mean + Osteopontin.RNA.se), 
                width = 0.3, position=position_dodge(.75)) + 
  geom_bar(data = mean_rediff, aes(x = interaction(Days.CM, Passage), y = Osteopontin.RNA.mean.mean), 
           fill = c(p0f, p0f, p5f, p5f), alpha = 0.6, stat = 'identity', color = '#000000') +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  ylab('Mean Osteopontin RNA per Cell') + xlab(NULL)


# chond_rediff_correlation <- 
#   ggplot(subset(boots, Purpose == 'chondrocyte redifferentiation core' & Substrate == 'agarose'), 
#                                    aes(x = interaction(Days.CM, Passage), y = actual.rho, fill = Donor)) +
#   scale_x_discrete(limits = c('1.0', '14.0', '1.5', '14.5'), labels = c('Day 1\nP0', 'Day 14\nP0', 'Day 1\nP5', 'Day 14\nP5')) +
#   scale_fill_manual(breaks = c('1B', '2', '3'), 
#                     values = c(donA, donB, donC)) +
#   geom_errorbar(aes(ymin = actual.rho, ymax = upper.rho), 
#                 width = 0.3, position=position_dodge(.75)) + 
#   geom_bar(width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE) +
#   geom_bar(width = 0.75, stat = 'identity', position = position_dodge(), color = '#000000', show_guide = FALSE, alpha = 0) + 
#   geom_bar(data = mean_rediff, aes(x = interaction(Days.CM, Passage), y = AggGAPDH.RNA.cor.mean), fill = c(p0f, p0f, p5f, p5f), 
#            alpha = 0.6, stat = 'identity', color = '#000000') +
#   scale_y_continuous(limits=c(0, 1)) +
#   ylab(expression(paste(rho, ' for GAPDH vs Aggrecan', sep = "\ "))) + xlab(NULL)

#chondrocyte gel scatter plots, Donor 1B only
chondGel<-subset(forscattersGAPDH10, Purpose == 'chondrocyte redifferentiation core' & Substrate=="agarose" & Donor == '1B')
chondMeans <- ddply(chondGel, c("Passage", "Days.CM"), summarise, mean.AGG=mean(cy.RNA), sem.AGG=sd(cy.RNA)/sqrt(length(cy.RNA)), mean.GAP=mean(GAPDH.RNA), sem.GAP=sd(GAPDH.RNA)/sqrt(length(GAPDH.RNA)) )
chondMeans <- transform(chondMeans, lower.AGG=mean.AGG-sem.AGG, upper.AGG=mean.AGG+sem.AGG, lower.GAP=mean.GAP-sem.GAP, upper.GAP=mean.GAP+sem.GAP)

chond_rediff_scatter<-  
  ggplot() + 
  geom_point(data=chondGel, aes(x = GAPDH.RNA, y = cy.RNA, color=Passage, fill=Passage), alpha = 0.25, pch=21) +
  facet_grid(.~Days.CM) + 
  geom_smooth(data=chondGel, aes(x = GAPDH.RNA, y = cy.RNA, color=Passage), method = 'lm', se = FALSE, size = 2) +
  geom_errorbar(data=chondMeans, aes(x=mean.GAP, ymin=lower.AGG, ymax=upper.AGG), size=1) +
  geom_errorbarh(data=chondMeans, aes(x=mean.GAP, y=mean.AGG, xmin=lower.GAP, xmax=upper.GAP), size=1.5) +
  geom_point(data=chondMeans, aes(x=mean.GAP, y=mean.AGG, fill=Passage),pch=21, color="black", size=2) +
  ylab('Aggrecan RNA') + xlab('GAPDH RNA') +
  scale_fill_manual(values=c( p0f, p5f))+
  scale_color_manual(values=c( p0f, p5f)) +
  theme(legend.position = 'none')

# chondrocyte viability in gels

chondViability$Replicate.. <- as.factor(chondViability$Donor)
chondViability$Donor <- as.factor(chondViability$Donor)
chondViability$Day <- as.factor(chondViability$Day)
chondViability$Passage <- as.factor(chondViability$Passage)
chondViability$X.Live...Filter <- as.numeric(chondViability$X.Live...Filter)

chondViability<-subset(chondViability, Day %in% c("1", "14") )

chondViabilityMeans <- ddply(chondViability, c("Passage", "Day", "Donor"), summarise,
                        mean.X.Live=mean(X.Live...Filter), sem.X.Live=sd(X.Live...Filter)/sqrt(length(X.Live...Filter)))
chondViabilityMeans.mean  <- ddply(chondViabilityMeans, c("Passage", "Day"), summarise,
                                   mean.mean.X.Live=mean(mean.X.Live))

chond_viability<-
  ggplot(data=chondViabilityMeans, aes(x = interaction(Passage, Day), y = mean.X.Live, fill = Donor))+
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge()) + 
  geom_errorbar(data=chondViabilityMeans, aes(ymin = mean.X.Live - sem.X.Live, ymax = mean.X.Live + sem.X.Live), 
                width = 0.3, position=position_dodge(.75)) + 
  scale_fill_manual(values=c(donA, donC)) + 
  geom_bar(data = chondViabilityMeans.mean, aes(x = interaction(Passage, Day), y =  mean.mean.X.Live),
           fill=c(p0f, p5f, p0f, p5f), color=c(p0c,  p5c, p0c,  p5c), alpha = 0.6, stat = 'identity') +
  ylab('Fraction Live Cells ') + xlab('Passage*Day')


##printing plots
plots <- c(#'chond_rediff_correlation',
           'dediff_aggrecan_means',
           'dediff_gapdh_means',
           'dediff_aggrecanPerGAPDH_pooled',
           'dediff_aggrecanPerGAPDH_perCell',
           'chond_volume',
           'chond_viability',  
           'chond_rediff_meanGAPDH',
           'chond_rediff_meanAgg',
           #'chond_rediff_correlation',
           'chond_rediff_scatter',
           'dediff_aggrecanPerGAPDH_pooled_nolegend',
           'chond_rediff_AggPerGaPDHpooled',
           'chond_rediff_AggPerGaPDHperCell')

plotWidth <- c(#17.4, #'chond_rediff_correlation',
               17.4, #'dediff_aggrecan_means',
               17.4, #'dediff_gapdh_means',
               17.4, #'dediff_aggrecanPerGAPDH_pooled',
               17.4, #'dediff_aggrecanPerGAPDH_perCell',
               17.4, #'chond_volume',
               17.4, #'chond_viability',
               17.4, #'chond_rediff_meanGAPDH',
               17.4, #'chond_rediff_meanAgg',
              # 17.4, #'chond_rediff_correlation',
               36, #'chond_rediff_scatter',
               17.4,#'dediff_aggrecanPerGAPDH_pooled_nolegend'
               17.4,
               17.4) 
systemType <- 'Mac'

if(systemType=="Mac"){
  
  directory <- paste('./Figure5/graphs/', Sys.Date(),'/', sep = '')
  
  if(directory %in% dir() == FALSE) dir.create(directory) 
  for (i in seq(along = plots)){
    savefile <- paste(Sys.Date(),'_',as.name(plots[i]),'.pdf',sep = '')
    ggsave(filename = savefile, plot = get(plots[i]), device = pdf,
           path = directory, scale = 1, 
           width = 20, height = 10, units = "cm", dpi = 300)
  }
}


if(systemType=="Windows"){
  setwd("Figure5\\")
  directory <- paste('graphs\\', Sys.Date(),'\\', sep = "")
  if(directory %in% dir() == FALSE) dir.create(directory)
  for (i in seq(along = plots)){
    savefile <- paste(directory, Sys.Date(),'_',as.name(plots[i]),'.pdf',sep = '')
    ggsave(filename = savefile, plot = get(plots[i]), device = pdf,
           scale = 1, 
           width = plotWidth[i], height = 14, units = "cm", dpi = 300)
  }
  setwd("..\\")
}

