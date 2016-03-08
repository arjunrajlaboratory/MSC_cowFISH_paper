#Figures for paper
systemType <- 'Mac'
if(!exists('forscattersGAPDH10')) source('./data_import/calculate_all_stats_forpaper.R')
theme_set(theme_bw(base_size = 18))

#correlation
singleCellScore$plotBin<-singleCellScore$plotScore=="Extracellular Staining"
# cor(singleCellScore$cy.RNA, singleCellScore$plotBin)

#MSC single cell stuff
donor_colors_gray <- c('grey77', 'grey51', 'grey36')
MSCfillColor<-'mediumseagreen'

singleCellScore_density<-subset(singleCellScore, Donor %in% c('9'))

#N
sc_N<-ggplot()+
  geom_bar(data=singleCellScore_stats, aes(y=N, x=plotScore, fill=Donor), width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE)+
   scale_fill_manual(values = donor_colors_gray, labels = c('H', 'I', 'J'))+
  ylab('Number of Cells')

# N from randomly sampled images
sampledSingleCell<-subset(singleCellScore, objArrayNum<14)

unscreenedStainingFractions<-ggplot()+
  geom_bar(data = sampledSingleCell,
         aes(x =Donor, fill = plotScore),
         position = "fill")




#AGG/GAPDH
sc_AGGperGAPDH<-ggplot()+
  geom_bar(data=singleCellScore_stats, aes(y=AggperGAPDH.RNA.mean, x=plotScore, fill=Donor), width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE)+
  geom_errorbar(data=singleCellScore_stats, 
                aes(x=plotScore, ymin = AggperGAPDH.RNA.mean - .001, ymax=AggperGAPDH.RNA.mean + AggperGAPDH.RNA.se, fill=Donor), 
                width = 0.3, position=position_dodge(.75)) + 
  scale_fill_manual(values = donor_colors_gray, labels = c('H', 'I', 'J'))+
  geom_bar(data=means_singleCellScore, aes(y=AggPerGAPDH.mean.mean, x=plotScore),  fill=MSCfillColor, width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE, alpha=.5, color='black')+
  xlab('Staining')+
  ylab('Aggrecan/GAPDH')

#AGG
sc_AGG<-ggplot()+
  geom_bar(data=singleCellScore_stats, aes(y=Aggrecan.RNA.mean, x=plotScore, fill=Donor), width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE)+
  geom_errorbar(data=singleCellScore_stats, 
                aes(x=plotScore, ymin = Aggrecan.RNA.mean - .001, ymax=Aggrecan.RNA.mean + Aggrecan.RNA.se, fill=Donor), 
                width = 0.3, position=position_dodge(.75)) + 
  scale_fill_manual(values = donor_colors_gray, labels = c('H', 'I', 'J'))+
  geom_bar(data=means_singleCellScore, aes(y=Aggrecan.RNA.mean.mean, x=plotScore),  fill=MSCfillColor, width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE, alpha=.5, color='black')+
  xlab('Staining')+
  ylab('Aggrecan')


sc_AGG_density<-ggplot()+
  geom_density(data=singleCellScore_density, aes(x=cy.RNA, fill=plotScore, color=plotScore), alpha=.5)+
  xlab('Aggrecan RNA')

sc_AGGperGAPDH_density<-ggplot()+
  geom_density(data=singleCellScore_density, aes(x=cy.RNA/GAPDH.RNA, fill=plotScore, color=plotScore), alpha=.5)+
  xlab('Aggracan/GAPDH (per cell)')

sc_AGGperOPN_density<-ggplot()+
  geom_density(data=singleCellScore_density, aes(x=cy.RNA/tmr.RNA, fill=plotScore, color=plotScore), alpha=.5)+
  xlab('Aggracan/Osteopontin (per cell)')+
  scale_x_continuous(limits=c(0,20))


#GAPDH
sc_GAPDH<-ggplot()+
  geom_bar(data=singleCellScore_stats, aes(y=GAPDH.RNA.mean, x=plotScore, fill=Donor), width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE)+
  geom_errorbar(data=singleCellScore_stats, 
                aes(x=plotScore, ymin = GAPDH.RNA.mean - .001, ymax=GAPDH.RNA.mean + GAPDH.RNA.se, fill=Donor), 
                width = 0.3, position=position_dodge(.75)) +   
  scale_fill_manual(values = donor_colors_gray, labels = c('H', 'I', 'J'))+
  geom_bar(data=means_singleCellScore, aes(y=GAPDH.RNA.mean.mean, x=plotScore),  fill=MSCfillColor, width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE, alpha=.5, color='black')+
  xlab('Staining')+
    ylab('GAPDH')

#OPN
sc_OPN<-ggplot()+
  geom_bar(data=singleCellScore_stats, aes(y=Osteopontin.RNA.mean, x=plotScore, fill=Donor), width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE)+
  geom_errorbar(data=singleCellScore_stats, 
                aes(x=plotScore, ymin = Osteopontin.RNA.mean - .001, ymax=Osteopontin.RNA.mean + Osteopontin.RNA.se, fill=Donor), 
                width = 0.3, position=position_dodge(.75)) +   
  scale_fill_manual(values = donor_colors_gray, labels = c('H', 'I', 'J'))+
  geom_bar(data=means_singleCellScore, aes(y=Osteopontin.RNA.mean.mean, x=plotScore),  fill=MSCfillColor, width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE, alpha=.5, color='black')+
  xlab('Staining') +
  ylab('Osteopontin')

#AGG/OPN
sc_AGGperOPN<-ggplot()+
  geom_bar(data=singleCellScore_stats, aes(y=AGGperOPN, x=plotScore, fill=Donor), width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE)+
  geom_errorbar(data=singleCellScore_stats, 
                aes(x=plotScore, ymin = AGGperOPN - .001, ymax=AGGperOPN + AGGperOPN.se, fill=Donor), 
                width = 0.3, position=position_dodge(.75)) +     
  scale_fill_manual(values = donor_colors_gray, labels = c('H', 'I', 'J'))+
  geom_bar(data=means_singleCellScore, aes(y=AGGperOPN.mean, x=plotScore),  fill=MSCfillColor, width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE, alpha=.5, color='black')+
  xlab('Staining')+
  ylab('Aggrecan/Osteopontin')
  
ggplot()+
  geom_bar(data=singleCellScore_stats, aes(y=AGGperOPN.median, x=plotScore, fill=Donor), width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE)+
  scale_fill_manual(values = donor_colors_gray, labels = c('H', 'I', 'J'))+
  geom_bar(data=means_singleCellScore, aes(y=AGGperOPN.median.mean, x=plotScore),  fill=MSCfillColor, width = 0.75, stat = 'identity', position = position_dodge(), show_guide = FALSE, alpha=.5, color='black')+
  xlab('Staining')+
  ggtitle("AGG/OPN (median): MSC Day 7 CM+ Single Cell - Donors 8 9 10")


sc_scattersAGG_OPN<-ggplot(data=singleCellScore, aes(x=tmr.RNA, y=cy.RNA, color=plotScore))+
  geom_point()+
  geom_smooth(method=lm, se=TRUE, )+
  facet_grid(.~Donor)

#roc plots
library(pROC)  

#for MSCs
forROC_MSC<-subset(singleCellScore, Donor %in% c('8','9','10'))
#forROC_MSC<-subset(singleCellScore, Donor %in% c('10'))
forROC_MSC$AGGperOPN=forROC_MSC$cy.RNA/forROC_MSC$tmr.RNA

forROC_MSC$binary<-ifelse(forROC_MSC$Score>3, 1, 0)

MSCoutcomeAGG<-roc(response=forROC_MSC$binary, predictor=forROC_MSC$cy.RNA )
MSCoutcomeGAPDH<-roc(response=forROC_MSC$binary, predictor=forROC_MSC$GAPDH.RNA )
MSCoutcomeOPN<-roc(response=forROC_MSC$binary, predictor=forROC_MSC$tmr.RNA )
MSCoutcomeAGGperOPN<-roc(response=forROC_MSC$binary, predictor=forROC_MSC$AGGperOPN )
MSCoutcomeAGGperGAPDH<-roc(response=forROC_MSC$binary, predictor=(forROC_MSC$cy.RNA/forROC_MSC$GAPDH.RNA) )
MSCoutcomeLPL<-roc(response=forROC_MSC$binary, predictor=forROC_MSC$alexa.RNA )

#all MSCs
directory <- paste('./Figure2/graphs/', Sys.Date(),'/', sep = '')
if(directory %in% dir() == FALSE) dir.create(directory) 
pdf(paste(directory,'ROCcurves.pdf', sep = ''))
t1<-plot.roc(MSCoutcomeAGG, col="black")
t2<-lines.roc(MSCoutcomeOPN, col="maroon2")
t3<-lines.roc(MSCoutcomeGAPDH, col="mediumslateblue")
#t4<-lines.roc(MSCoutcomeAGGperOPN, col="#7570b3")
#t5<-lines.roc(MSCoutcomeAGGperGAPDH, col="#e7298a")
#t6<-lines.roc(MSCoutcomeLPL, col="limegreen")

t1<-plot.roc(MSCoutcomeAGG, col="black")
#t2<-lines.roc(MSCoutcomeOPN, col="#e6ab02")
#t3<-lines.roc(MSCoutcomeGAPDH, col="#66a61e")
t4<-plot.roc(MSCoutcomeAGGperOPN, col="tomato", add = TRUE)
t5<-plot.roc(MSCoutcomeAGGperGAPDH, col="turquoise3", add = TRUE)
#t6<-lines.roc(MSCoutcomeLPL, col="#d95f02")
dev.off()

forROC_MSCobj=roc.test(t5,t1, method="d")
forROC_MSCobj$p.value

#find the threshold for aggrecan RNA that maximizes Youden J statistic = 405.5 RNA
thresh<-coords(MSCoutcomeAGG, "b", ret="t", best.method="youden") # default

sc_AGG_density_threshold<-ggplot()+
  geom_histogram(data=singleCellScore_density, aes(x=cy.RNA, fill=plotScore, color=plotScore), position=, alpha=.5)+
  xlab('Aggrecan RNA')+
  facet_grid(.~plotScore)

singleCellScoreThresh=singleCellScore
singleCellScoreThresh$thresholdAgg405<-singleCellScoreThresh$cy.RNA>thresh

truePositives<-subset(singleCellScoreThresh, thresholdAgg405 & plotScore == "Extracellular Staining")
falsePositives<-subset(singleCellScoreThresh, thresholdAgg405 & plotScore != "Extracellular Staining")
trueNegatives<-subset(singleCellScoreThresh, !thresholdAgg405 & plotScore != "Extracellular Staining")
falseNegatives<-subset(singleCellScoreThresh, !thresholdAgg405 & plotScore == "Extracellular Staining")

summary(truePositives$Donor)
summary(falsePositives$Donor)
summary(trueNegatives$Donor)
summary(falseNegatives$Donor)


posECM<-subset(singleCellScoreThresh, plotScore=="Extracellular Staining")
negECM<-subset(singleCellScoreThresh, plotScore==" No Extracellular Staining")

fauxSort1<-
  ggplot()+
  geom_histogram(data=posECM, aes(x='_Unsorted', y= -..count.., fill=plotScore) )+
  geom_histogram(data=negECM, aes(x='_Unsorted', y= ..count.., fill=plotScore) )+
  geom_histogram(data=posECM, aes(x=thresholdAgg405, y= -..count.., fill=plotScore) )+
  geom_histogram(data=negECM, aes( x=thresholdAgg405, y= ..count.., fill=plotScore) )+
  scale_fill_manual(values = c("darkblue", "blue")) +
  scale_y_reverse() +
  scale_x_discrete(limits = c('TRUE', 'FALSE', '_Unsorted')) + 
  coord_flip() +
  #facet_grid(.~Donor)+
  xlab(thresh)

fauxSortPercentile<-ggplot() + 
  geom_bar(data = singleCellScoreThresh,
           aes(x = factor(thresholdAgg405),fill = interaction(thresholdAgg405,plotScore)),
           position = "fill")+
  geom_bar(data = singleCellScoreThresh,
           aes(x = "_unsorted",fill = plotScore),
           position = "fill")+
  facet_grid(.~Donor)

##printing plots
plots <- c(#'sc_N',
           #'sc_AGG',
           #'sc_GAPDH',
           #'sc_OPN',
           #'sc_AGGperGAPDH',
           #'sc_AGGperOPN',
           #'sc_scattersAGG_OPN',
           'fauxSort1',
           'fauxSortPercentile',
           'sc_AGG_density',
           #'sc_AGGperGAPDH_density',
           #'sc_AGGperOPN_density',
           'unscreenedStainingFractions'
           )

plotWidth <- c(17.4,
               17.4, 
               17.4, 
               17.4, 
               17.4, 
               17.4,
               34, #scatters
               34,
               34,
               34,
               34,
               34, 
               17
               ) 


if(systemType=="Mac"){
  
  directory <- paste('./Figure2/graphs/', Sys.Date(),'/', sep = '')
  
  if(directory %in% dir() == FALSE) dir.create(directory) 
  for (i in seq(along = plots)){
    savefile <- paste(Sys.Date(),'_',as.name(plots[i]),'.pdf',sep = '')
    ggsave(filename = savefile, plot = get(plots[i]), device = pdf,
           path = directory, scale = 1, 
           width = 20, height = 10, units = "cm", dpi = 300)
  }
}


if(systemType=="Windows"){
  setwd("Figure2\\")
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
