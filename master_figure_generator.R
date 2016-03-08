#NOTE: running this script will generate a great deal of warnings, which can be ignored. They have to do with loading specific packages and 
# introduction of NAs into data structures that don't have certain factor levels. 
#  It requires a folder in Dropbox called "cow fishing" and two subfolders called "cowFISH_figures" (which contains all the scripts and where the graphs
#  are generated) and "all_csvs" (which contains all of the data). 

#required packages:
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


#generate figure 1
source('./Figure1/Figure1.R')

#generate figure 2
source('./Figure2/Figure2.R')

#generate figure 3
source('./Figure3/Figure3.R')

#generate figure 4
source('./Figure4/Figure4.R')

#generate figure 5
source('./Figure4/Figure5.R')

#generate supplementary figure 1
source('./FigureS1/FigureS1.R')

#generate supplementary figure 2
source('./FigureS2/FigureS2.R')

#generate supplementary figure 3
source('./FigureS3/FigureS3.R')

#generate supplementary figure 4
source('./FigureS4/FigureS4.R')

rm(list = ls())
