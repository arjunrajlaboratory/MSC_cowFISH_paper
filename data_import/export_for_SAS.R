##### Export data and statistics to cvs for use in SAS
#export data and stats for use in SAS

forscattersGAPDH10$AGGperGAPDH<-forscattersGAPDH10$cy.RNA/forscattersGAPDH10$GAPDH.RNA
forscattersGAPDH10$log.cy.RNA<-log1p(forscattersGAPDH10$cy.RNA)
forscattersGAPDH10$log.GAPDH.RNA<-log1p(forscattersGAPDH10$GAPDH.RNA)
forscattersGAPDH10$log.AGGperGAPDH<-log1p(forscattersGAPDH10$AGGperGAPDH)

mscGel<-subset(forscattersGAPDH10, Purpose == 'MSCs het differentiation in gels core' & Days.CM>0)
write.csv(mscGel, "Figure3\\forSAS\\MSCDiffCore.csv")


chondGel<-subset(forscattersGAPDH10, Purpose == 'chondrocyte redifferentiation core' & Substrate=="agarose")
write.csv(chondGel, "Figure4\\forSAS\\ChondRediffCore.csv")


chondGelGlass<-subset(forscattersGAPDH10, Purpose == 'chondrocyte redifferentiation core' & Days.CM<2)
write.csv(chondGelGlass, "Figure4\\forSAS\\ChondRediffCore_gelAndGlass.csv")

chondMono_rediff<-subset(forscattersGAPDH10,  Purpose=='chondrocyte redifferentiation core' & Substrate == "glass" & Passage %in% c(0,1,3,5,7,9))
write.csv(chondMono_rediff, "Figure4\\forSAS\\ChondMonolayer_Rediff.csv")

chondMono_dediff<-subset(forscattersGAPDH10,  Purpose=='chondrocyte dedifferentiation core' & Substrate == "glass" & Passage %in% c(0,1,3,5,7,9))
write.csv(chondMono_dediff, "Figure4\\forSAS\\ChondMonolayer_Dediff.csv")

chondMono_all<- rbind(chondMono_rediff, chondMono_dediff)
write.csv(chondMono_all, "Figure4\\forSAS\\ChondMonolayer_all.csv")


mscGel_allstats <- subset(allstats, Purpose == 'MSCs het differentiation in gels core' & Days.CM>0)
mscGel_allstats$label<-paste( mscGel_allstats$Days.CM, mscGel_allstats$Media, sep="_")
write.csv(mscGel_allstats, "Figure3\\forSAS\\MSCDiffCore_allstats.csv")

chondGel_allstats <- subset(allstats, Purpose == 'chondrocyte redifferentiation core' & Substrate=="agarose")
chondGel_allstats$label<-paste(chondGel_allstats$Days.CM, chondGel_allstats$Passage, sep="_")
write.csv(chondGel_allstats, "Figure4\\forSAS\\ChondRediffCore_allstats.csv")

chondMono_rediff_allstats <- subset(allstats,  Purpose=='chondrocyte redifferentiation core' & Substrate == "glass" & Passage %in% c(0,1,3,5,7,9))
write.csv(chondMono_rediff_allstats, "Figure4\\forSAS\\ChondMonolayer_Rediff_allstats.csv")

chondMono_dediff_allstats <- subset(allstats,  Purpose=='chondrocyte dedifferentiation core' & Substrate == "glass" & Passage %in% c(0,1,3,5,7,9))
write.csv(chondMono_dediff_allstats, "Figure4\\forSAS\\ChondMonolayer_Dediff_allstats.csv")


chondMono_all_allstats <- rbind(chondMono_rediff_allstats, chondMono_dediff_allstats)
write.csv(chondMono_all_allstats, "Figure4\\forSAS\\ChondMonolayer_all_allstats.csv")
