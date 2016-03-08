#Figures for paper
theme_set(theme_bw(base_size = 32))
if(!exists('forscattersGAPDH10')) source('./data_import/calculate_all_stats_forpaper.R')
systemType <- 'Mac'

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

# PCR

chond_PCR$AGGperGAP <- as.numeric(chond_PCR$AGGperGAP)
chond_PCR$OPNperGAP <- as.numeric(chond_PCR$OPNperGAP)

melted<- melt(chond_PCR, id=c( "Donor", "Passage"))

#Aggrecan

melted.AGG <- subset(melted, variable == "AGGperGAP")
means <- ddply(melted.AGG, c( "Passage", "variable"), summarise, mean=mean(value))

melted.AGG <- subset(melted.AGG, variable == "AGGperGAP")
melted.AGG$value<-as.numeric(melted.AGG$value)
means.AGG <- ddply(melted.AGG, c("Passage"), summarise, mean=mean(value))

means.AGGsem <- ddply(melted.AGG, c("Passage"), summarise, mean=mean(value), sem=sd(value)/sqrt(length(value)))
means.AGGsem <- transform(means.AGGsem, lower=mean-sem, upper=mean+sem)

P1.AGG.mean<-means.AGG$mean[2]

norm.AGG<- melted.AGG
norm.AGG$value <- norm.AGG$value / P1.AGG.mean
norm.AGG.means <- ddply(norm.AGG, c("Passage"), summarise, mean=mean(value))

norm.means.AGGsem <- ddply(norm.AGG, c("Passage"), summarise, mean=mean(value), sem=sd(value)/sqrt(length(value)))
norm.means.AGGsem <- transform(norm.means.AGGsem, lower=mean-sem, upper=mean+sem)


#=========== LINE GRAPH
AGG_PCR <- ggplot() +
  geom_point(data=melted.AGG, aes(x=Passage, y=value,  color=Donor),  size =2) +
  geom_point(data=means.AGG, aes(x=Passage, y=mean), shape=21, size = 3.7) +
  geom_line(data=means.AGG, aes(x=Passage, y=mean)) +
  ylab("ACAN per GAPDH (PCR)")


#Osteopontin

melted.OPN <- subset(melted, variable == "OPNperGAP")
means <- ddply(melted.OPN, c( "Passage", "variable"), summarise, mean=mean(value))

melted.OPN <- subset(melted.OPN, variable == "OPNperGAP")
melted.OPN$value<-as.numeric(melted.OPN$value)
means.OPN <- ddply(melted.OPN, c("Passage"), summarise, mean=mean(value))

means.OPNsem <- ddply(melted.OPN, c("Passage"), summarise, mean=mean(value), sem=sd(value)/sqrt(length(value)))
means.OPNsem <- transform(means.OPNsem, lower=mean-sem, upper=mean+sem)

P1.OPN.mean<-means.OPN$mean[2]

norm.OPN<- melted.OPN
norm.OPN$value <- norm.OPN$value / P1.OPN.mean
norm.OPN.means <- ddply(norm.OPN, c("Passage"), summarise, mean=mean(value))

norm.means.OPNsem <- ddply(norm.OPN, c("Passage"), summarise, mean=mean(value), sem=sd(value)/sqrt(length(value)))
norm.means.OPNsem <- transform(norm.means.OPNsem, lower=mean-sem, upper=mean+sem)


#=========== LINE GRAPH
OPN_PCR <- ggplot() +
  geom_point(data=melted.OPN, aes(x=Passage, y=value,  color=Donor),  size =3) +
  geom_point(data=means.OPN, aes(x=Passage, y=mean), shape=21, size = 3.7) +
  geom_line(data=means.OPN, aes(x=Passage, y=mean)) +
  ylab("OPN per GAPDH (PCR)")

# area quantification
area$Passage<-as.factor(area$Passage)
area$Donor<-as.factor(area$Donor)

area_means <- ddply(area, c("Passage", "Donor"), summarise,
                    area_mean=mean(AreaMicrons), 
                    area_sem=sd(AreaMicrons)/sqrt(length(AreaMicrons)))
area_means.mean <- ddply(area_means, c("Passage"), summarise,
                         area_mean.mean=mean(area_mean))

chond_area<-
  ggplot(data=area_means, aes(x = Passage, y = area_mean, fill = Donor))+
  geom_errorbar(data=area_means, aes(ymin = area_mean - area_sem, ymax = area_mean + area_sem), 
                width = 0.3, position=position_dodge(.75)) + 
  geom_bar(width = 0.75, stat = 'identity', position = position_dodge()) + 
  scale_fill_manual(values=c(donB, donD)) + 
  geom_bar(data = area_means.mean, aes(x = Passage, y = area_mean.mean),
           fill=c(p0f, p1f, p5f), color=c(p0c, p1c, p5c), alpha = 0.6, stat = 'identity') +
  ylab('Chondrocyte Spread Area (um^2)') + xlab('Passage')+
  ylim(c(0,1750))


##printing plots
plots_square <- c()
plots_rectangle <- c(
  'AGG_PCR',
  'OPN_PCR',
  'chond_area'
  )


for (i in seq(along = plots_square)){
  savefile <- paste(Sys.Date(),'_',as.name(plots_square[i]),'.pdf',sep = '')
  ggsave(filename = savefile, plot = get(plots_square[i]), device = pdf,
         path = directory, scale = 1, 
         width = 10, height = 10, units = "cm", dpi = 300)
}


if(systemType=="Mac"){
  
  directory <- paste('./FigureS4/graphs/', Sys.Date(),'/', sep = '')
  
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
  setwd("FigureS4\\")
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