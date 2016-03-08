library(stringr)
library(ggplot2)
library(reshape2)
library(RCurl)
library(dplyr)
library(plyr)
library(stats)
library(psych)
library(pROC)
library(XLConnect)
library(data.table)
library(dplyr)
library(tidyr)
library(DESeq)
library(gplots)
library(RColorBrewer)
library(ggfortify)

systemType <- "Mac"
gsheetdata<-read.csv(file = './all_csvs/metadata.csv', 
                     header=TRUE,sep=",",na.strings='N/A')

files<-list.files(pattern="^[SGAP]..?.?\\.csv$",recursive=TRUE,full.names=TRUE) 
filecodes <- str_split( files, pattern = '/|\\.')
alldata <- list()
for(i in seq(along = files) ) {
  
  x<-read.table( files[i] ,header=TRUE,sep=",")
  x <- subset(x, isGood != 'false')
  x <- subset(x, isGood != 0)
  x<-cbind(x,"Cell.Number"=c(1:nrow(x)))
  experimentid<- filecodes[[i]][length(filecodes[[i]])-1]
  #CHECK TO SEE IF IT ALREADY EXISTS, IF IT DOES, ASSIGN IT AS 
  x <- cbind( x, experimentid )
  
  
  rm(experimentid)
  alldata[[i]] <- x
}
# the following are various manipulations of the alldata list to get it into a graphing-friendly format
alldata <- melt(alldata, id.vars=c('objNum','isGood','experimentid',
                                   'Cell.Number'))

alldata1 <- merge(x=alldata,y=gsheetdata,by.x="experimentid",by.y="Experiment.ID")

forscatters<- dcast(alldata1, 'experimentid+Cell.Number+Pre.Fix.Stain+Post.Fix.Stain+
                    Cell.Type+Substrate+Loading+Media+Donor+Days.CM+Passage+Purpose~ variable',
                    value.var = 'value')
forscatters[['GAPDH.RNA']]<- as.numeric(forscatters[['GAPDH.RNA']])
forscatters[['cy.RNA']]<- as.numeric(forscatters[['cy.RNA']])
forscatters[['tmr.RNA']]<- as.numeric(forscatters[['tmr.RNA']])
forscatters[['alexa.RNA']]<- as.numeric(forscatters[['alexa.RNA']])


theme_set(theme_gray(base_size = 18))
rm(x,filecodes,files,i)

#chondrocyte diameter and area data
diam<- read.csv("./all_csvs/ChondrocyteDiameterData.csv")
area<- read.csv("./all_csvs/chondrocyteSpreadAreas.csv")

#viability data
chondViability <- read.csv("./all_csvs/ChondroLiveDead.csv")


#bootstrapped correlation coefficient statistics
boots<-read.csv("./all_csvs/boots.csv")

#single cell data from donors 8, 9, 10
A140<-read.csv("./all_csvs/A140_withScore.csv")
A144<-read.csv("./all_csvs/A144_withScore.csv")
A148<-read.csv("./all_csvs/A148_withScore.csv")

singleCellScore <- rbind(A140, A144, A148)

A156 <- cbind(read.csv("./all_csvs/A156.csv"), "Donor"=c("1C"))
A157 <- cbind(read.csv("./all_csvs/A157.csv"), "Donor"=c("2C"))
A158 <- cbind(read.csv("./all_csvs/A158.csv"), "Donor"=c("3C"))

newIHC <- read.csv("./all_csvs/2015_09_18AggrecanIHCScoring.csv")
FISH<-rbind(A156,A157, A158)

newIHC.FISH<-cbind(newIHC, FISH)
newIHC.FISH<-subset(newIHC.FISH, GAPDH.RNA>10 & isGood=='true')
newIHC.FISH<-subset(newIHC.FISH, !is.na(Binary.Low.High))


#colony data
colonypositions<-read.csv('./all_csvs/ColonyPositions.csv')
colonypositions$Colony.Number <- as.factor(colonypositions$Colony.Number)
levels(colonypositions$Colony.Number) <- c('A','B','C','D','E')

colonyspotcounts<-read.table("./all_csvs/G3.csv",header=TRUE,sep=",")
colonyspotcounts<-subset(colonyspotcounts, isGood==1)
colonyspotcounts<-merge(colonyspotcounts,colonypositions,'dataFile')
colonyspotcounts<-subset(colonyspotcounts, Colony.Number != 'E')
colonyspotcounts.summary <- ddply(colonyspotcounts, .(Colony.Number), summarize, 
                                  Aggrecan.RNA.mean = mean(cy.RNA), 
                                  Aggrecan.RNA.sd = sd(cy.RNA),
                                  Aggrecan.RNA.CV = sd(cy.RNA)/mean(cy.RNA),
                                  GAPDH.RNA.CV = sd(GAPDH.RNA)/mean(GAPDH.RNA),
                                  N = length(Colony.Number))

#cell tracking data
fish_spots <- read.csv('./all_csvs/G69.csv')
fish_spots2 <- read.csv('./all_csvs/G70.csv')
timelapse_FISH_registration <- read.csv('./all_csvs/timelapse_fish_cell_matchup.csv')
division_timing <- read.csv('./all_csvs/division_timing_left_well.csv')

#live dead for MSCs
MSC_LD<-read.csv('./all_csvs/LiveDead_MSC.csv')

#chondocyte dediff PCR
chond_PCR<-read.csv('./all_csvs/ChondroDediffPCR.csv')
