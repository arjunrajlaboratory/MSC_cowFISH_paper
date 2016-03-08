if(!exists('forscatters')) source('~/code/cowfish/acote_import.R')
forscattersGAPDH10 <- subset(forscatters, GAPDH.RNA >=10)

#calculate confidence intervals for rho squareds
library(boot)

boot.rsq <- function(mydata,i){
  d <- mydata[i,]
  GAPDH.RNA <- d$GAPDH.RNA
  cy.RNA <- d$cy.RNA
  rsq <- (cor(GAPDH.RNA,cy.RNA))
  return(rsq)
}

boots <- ddply(forscattersGAPDH10, .(experimentid), function(df) {
  results <- boot(data=df, statistic=boot.rsq, R=5000)
  orderedResults <- results$t[order(results$t)]
  actual <- results$t0
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  out <- c(actual,lower,upper)
  names(out) <- c('actual.rho','lower.rho','upper.rho');
  out
})

boots_MSC <- ddply(subset(forscattersGAPDH10, 
                          Purpose == 'MSCs het differentiation in gels core' & Days.CM %in% c(1, 4, 7, 14, 21)), 
                   .(Media, Days.CM, Donor), function(df) {
  results <- boot(data=df, statistic=boot.rsq, R=5000)
  orderedResults <- results$t[order(results$t)]
  actual <- results$t0
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  out <- c(actual,lower,upper)
  names(out) <- c('actual.rho','lower.rho','upper.rho');
  out
})

boots <- merge(x = boots, y = gsheetdata, by.x = 'experimentid', by.y = 'Experiment.ID')
write.csv(boots, file = '~/Dropbox/Cow FISHing/all_csvs/boots.csv')

write.csv(boots_MSC, file = '~/Dropbox/Cow FISHing/all_csvs/boots_MSC.csv')

