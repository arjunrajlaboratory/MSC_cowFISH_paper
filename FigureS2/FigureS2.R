#Figures for paper
if(!exists('forscattersGAPDH10')) source('./data_import/calculate_all_stats_forpaper.R')
theme_set(theme_bw(base_size = 18))
systemType <- 'Mac'

#correlation
singleCellScore$plotBin<-singleCellScore$plotScore=="Extracellular Staining"
cor(singleCellScore$cy.RNA, singleCellScore$plotBin)

#MSC single cell stuff
donor_colors_gray <- c('grey77', 'grey51', 'grey36')
MSCfillColor<-'mediumseagreen'
singleCellScore_density<-subset(singleCellScore, Donor %in% c('9', '8', '10'))

# N from randomly sampled images
sampledSingleCell<-subset(singleCellScore, objArrayNum<14)
N <- ddply(sampledSingleCell, .(plotScore, Donor), summarize, N = length(CellNumber))

fill_colors <- c('darkslategray', 'cyan3')

sc_AGG_density <- ggplot()+
  geom_density(data=singleCellScore_density, aes(x=cy.RNA, fill=plotScore, color=plotScore), alpha=.7)+
  xlab('Aggrecan RNA') + facet_grid(Donor~.) + 
  geom_vline(data = singleCellScore_stats, aes(xintercept=Aggrecan.RNA.mean, color = plotScore), 
             linetype="dashed", size=1) +
  scale_fill_manual(values = fill_colors, name = '', labels = c('Low Staining', 'High Staining')) + 
  scale_color_manual(values = fill_colors, guide = FALSE) +
  theme(legend.position = c(1,1), legend.justification = c(1,1))

sc_GAPDH_density <- ggplot()+
  geom_density(data=singleCellScore_density, aes(x=GAPDH.RNA, fill=plotScore, color=plotScore), alpha=.7)+
  xlab('GAPDH RNA') + facet_grid(Donor~.) +
  geom_vline(data = singleCellScore_stats, aes(xintercept=GAPDH.RNA.mean, color = plotScore), 
             linetype="dashed", size=1) +
  scale_fill_manual(values = fill_colors, name = '', labels = c('Low Staining', 'High Staining')) + 
  scale_color_manual(values = fill_colors, guide = FALSE) +
  theme(legend.position = 'none')

sc_OPN_density <- ggplot()+
  geom_density(data=singleCellScore_density, aes(x=tmr.RNA, fill=plotScore, color=plotScore), alpha=.7)+
  xlab('Osteopontin RNA') + facet_grid(Donor~.) + 
  geom_vline(data = singleCellScore_stats, aes(xintercept=Osteopontin.RNA.mean, color = plotScore), 
             linetype="dashed", size=1) +
  scale_fill_manual(values = fill_colors, name = '', labels = c('Low Staining', 'High Staining')) + 
  scale_color_manual(values = fill_colors, guide = FALSE) +
  theme(legend.position = 'none')

sc_AGGperGAPDH_density<-
  ggplot()+
  geom_density(data=singleCellScore_density, aes(x=cy.RNA/GAPDH.RNA, fill=plotScore, color=plotScore), alpha=.5)+
  xlab('Aggrecan/GAPDH (per cell)') + facet_grid(Donor~.) + 
  geom_vline(data = singleCellScore_stats, aes(xintercept=AggperGAPDH.RNA.mean, color = plotScore), 
             linetype="dashed", size=1) +
  scale_fill_manual(values = fill_colors, name = '', labels = c('Low Staining', 'High Staining')) + 
  scale_color_manual(values = fill_colors, guide = FALSE) +
  theme(legend.position = 'none')

sc_AGGperOPN_density<-ggplot()+
  geom_density(data=singleCellScore_density, aes(x=cy.RNA/tmr.RNA, fill=plotScore, color=plotScore), alpha=.5)+
  xlab('Aggrecan/Osteopontin (per cell)')+
#   scale_x_continuous(limits=c(0,20)) + 
  facet_grid(Donor~.) + 
  scale_fill_manual(values = fill_colors, name = '', labels = c('Low Staining', 'High Staining')) + 
  scale_color_manual(values = fill_colors, guide = FALSE) +
  geom_vline(data = singleCellScore_stats_AGGperOPN_removeinf, aes(xintercept=AGGperOPN, color = plotScore), 
             linetype="dashed", size=1) +
  theme(legend.position = 'none')


#Figure 2 R code renames columns in newIHC.FISH - if this isn't working, run the import script fresh
library(plyr)
newFISH.summarystats<-ddply(newIHC.FISH, .(Donor, Binary.Low.High), summarize,
Sox9.mean = mean(tmr.RNA),
COMP.mean = mean(alexa.RNA),
GAPDH.mean = mean(GAPDH.RNA),
AGG.mean = mean(cy.RNA)
)

Sox9_density<-ggplot()+
  geom_density(data=newIHC.FISH, aes(x=Sox9, fill=Binary.Low.High, color=Binary.Low.High), alpha=.5)+
  xlab('Sox9 (per cell)') + facet_grid(Donor~.) + 
  geom_vline(data = newFISH.summarystats, aes(xintercept=Sox9.mean, color = Binary.Low.High), 
             linetype="dashed", size=1) +
  scale_fill_manual(values = fill_colors, name = '', labels = c('Low Staining', 'High Staining')) + 
  scale_color_manual(values = fill_colors, guide = FALSE) +
  theme(legend.position = 'none')


COMP_density<-ggplot()+
  geom_density(data=newIHC.FISH, aes(x=COMP, fill=Binary.Low.High, color=Binary.Low.High), alpha=.5)+
  xlab('COMP (per cell)') + facet_grid(Donor~.) + 
  geom_vline(data = newFISH.summarystats, aes(xintercept=COMP.mean, color = Binary.Low.High), 
             linetype="dashed", size=1) +
  scale_fill_manual(values = fill_colors, name = '', labels = c('Low Staining', 'High Staining')) + 
  scale_color_manual(values = fill_colors, guide = FALSE) +
  theme(legend.position = 'none')

GAPDH2_density<-ggplot()+
  geom_density(data=newIHC.FISH, aes(x=GAPDH.RNA, fill=Binary.Low.High, color=Binary.Low.High), alpha=.5)+
  xlab('GAPDH (per cell)') + facet_grid(Donor~.) + 
  geom_vline(data = newFISH.summarystats, aes(xintercept=GAPDH.mean, color = Binary.Low.High), 
             linetype="dashed", size=1) +
  scale_fill_manual(values = fill_colors, name = '', labels = c('Low Staining', 'High Staining')) + 
  scale_color_manual(values = fill_colors, guide = FALSE) +
  theme(legend.position = 'none')

AGG2_density<-ggplot()+
  geom_density(data=newIHC.FISH, aes(x=cy.RNA, fill=Binary.Low.High, color=Binary.Low.High), alpha=.5)+
  xlab('ACAN (per cell)') + facet_grid(Donor~.) + 
  geom_vline(data = newFISH.summarystats, aes(xintercept=AGG.mean, color = Binary.Low.High), 
             linetype="dashed", size=1) +
  scale_fill_manual(values = fill_colors, name = '', labels = c('Low Staining', 'High Staining')) + 
  scale_color_manual(values = fill_colors, guide = FALSE) +
  theme(legend.position = 'none')

#roc plots
library(pROC)  

#for MSCs
forROC_MSC<-subset(singleCellScore, Donor %in% c('8','9','10'))
#forROC_MSC<-subset(singleCellScore, Donor %in% c('10'))
forROC_MSC$AGGperOPN=forROC_MSC$cy.RNA/forROC_MSC$tmr.RNA

forROC_MSC$binary<-ifelse(forROC_MSC$Score>3, 1, 0)

MSCoutcomeAGG_don8<-roc(response=subset(forROC_MSC, Donor == 8)$binary, predictor=subset(forROC_MSC, Donor == 8)$cy.RNA )
MSCoutcomeGAPDH_don8<-roc(response=subset(forROC_MSC, Donor == 8)$binary, predictor=subset(forROC_MSC, Donor == 8)$GAPDH.RNA )
MSCoutcomeOPN_don8<-roc(response=subset(forROC_MSC, Donor == 8)$binary, predictor=subset(forROC_MSC, Donor == 8)$tmr.RNA )
MSCoutcomeAGGperOPN_don8<-roc(response=subset(forROC_MSC, Donor == 8)$binary, predictor=subset(forROC_MSC, Donor == 8)$AGGperOPN )
MSCoutcomeAGGperGAPDH_don8<-roc(response=subset(forROC_MSC, Donor == 8)$binary, predictor=(subset(forROC_MSC, Donor == 8)$cy.RNA/subset(forROC_MSC, Donor == 8)$GAPDH.RNA) )

MSCoutcomeAGG_don9<-roc(response=subset(forROC_MSC, Donor == 9)$binary, predictor=subset(forROC_MSC, Donor == 9)$cy.RNA )
MSCoutcomeGAPDH_don9<-roc(response=subset(forROC_MSC, Donor == 9)$binary, predictor=subset(forROC_MSC, Donor == 9)$GAPDH.RNA )
MSCoutcomeOPN_don9<-roc(response=subset(forROC_MSC, Donor == 9)$binary, predictor=subset(forROC_MSC, Donor == 9)$tmr.RNA )
MSCoutcomeAGGperOPN_don9<-roc(response=subset(forROC_MSC, Donor == 9)$binary, predictor=subset(forROC_MSC, Donor == 9)$AGGperOPN )
MSCoutcomeAGGperGAPDH_don9<-roc(response=subset(forROC_MSC, Donor == 9)$binary, predictor=(subset(forROC_MSC, Donor == 9)$cy.RNA/subset(forROC_MSC, Donor == 9)$GAPDH.RNA) )

MSCoutcomeAGG_don10<-roc(response=subset(forROC_MSC, Donor == 10)$binary, predictor=subset(forROC_MSC, Donor == 10)$cy.RNA )
MSCoutcomeGAPDH_don10<-roc(response=subset(forROC_MSC, Donor == 10)$binary, predictor=subset(forROC_MSC, Donor == 10)$GAPDH.RNA )
MSCoutcomeOPN_don10<-roc(response=subset(forROC_MSC, Donor == 10)$binary, predictor=subset(forROC_MSC, Donor == 10)$tmr.RNA )
MSCoutcomeAGGperOPN_don10<-roc(response=subset(forROC_MSC, Donor == 10)$binary, predictor=subset(forROC_MSC, Donor == 10)$AGGperOPN )
MSCoutcomeAGGperGAPDH_don10<-roc(response=subset(forROC_MSC, Donor == 10)$binary, predictor=(subset(forROC_MSC, Donor == 10)$cy.RNA/subset(forROC_MSC, Donor == 10)$GAPDH.RNA) )


#all MSCs
directory <- paste('./FigureS2/graphs/', Sys.Date(),'/', sep = '')
if(directory %in% dir() == FALSE) dir.create(directory) 
pdf(paste(directory,'ROCcurves_byDonor.pdf', sep = ''))

plot.roc(MSCoutcomeAGG_don8, col="black")
lines.roc(MSCoutcomeOPN_don8, col="maroon2")
lines.roc(MSCoutcomeGAPDH_don8, col="mediumslateblue")
plot.roc(MSCoutcomeAGG_don8, col="black")
lines.roc(MSCoutcomeAGGperOPN_don8, col="tomato")
lines.roc(MSCoutcomeAGGperGAPDH_don8, col="turquoise3")

plot.roc(MSCoutcomeAGG_don9, col="black")
lines.roc(MSCoutcomeOPN_don9, col="maroon2")
lines.roc(MSCoutcomeGAPDH_don9, col="mediumslateblue")
plot.roc(MSCoutcomeAGG_don9, col="black")
lines.roc(MSCoutcomeAGGperOPN_don9, col="tomato")
lines.roc(MSCoutcomeAGGperGAPDH_don9, col="turquoise3")

plot.roc(MSCoutcomeAGG_don10, col="black")
lines.roc(MSCoutcomeOPN_don10, col="maroon2")
lines.roc(MSCoutcomeGAPDH_don10, col="mediumslateblue")
plot.roc(MSCoutcomeAGG_don10, col="black")
lines.roc(MSCoutcomeAGGperOPN_don10, col="tomato")
lines.roc(MSCoutcomeAGGperGAPDH_don10, col="turquoise3")
dev.off()


#find the threshold for aggrecan RNA that maximizes Youden J statistic = 405.5 RNA

##printing plots
plots <- c('sc_AGG_density',
           'sc_AGGperGAPDH_density',
           'sc_AGGperOPN_density',
           'sc_GAPDH_density',
           'sc_OPN_density',
           'Sox9_density',
         'COMP_density',
         'GAPDH2_density'
           )

if(systemType=="Mac"){
  
  directory <- paste('./FigureS2/graphs/', Sys.Date(),'/', sep = '')
  
  if(directory %in% dir() == FALSE) dir.create(directory) 
  for (i in seq(along = plots)){
    savefile <- paste(Sys.Date(),'_',as.name(plots[i]),'.pdf',sep = '')
    ggsave(filename = savefile, plot = get(plots[i]), device = pdf,
           path = directory, scale = 1, 
           width = 20, height = 10, units = "cm", dpi = 300)
  }
}


if(systemType=="Windows"){
  setwd("FigureS2\\")
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
