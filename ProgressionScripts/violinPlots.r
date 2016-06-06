# Violin Plots
library(R.matlab)
library(vioplot)
library(ggplot2)

drugs <-c('Doxorubicin','Vincristine','Paclitaxel','Cisplatin')
Dox <-readMat('/Users/sdalin/Dropbox\ (MIT)/Biology\ PhD/2016/Hemann\ Lab/CR.CS/Creating\ Resistant\ Cells/Round\ 3/Mondays/DoxPIDataFiltered.mat')

Vin <-readMat('/Users/sdalin/Dropbox\ (MIT)/Biology\ PhD/2016/Hemann\ Lab/CR.CS/Creating\ Resistant\ Cells/Round\ 3/Mondays/VinPIDataFiltered.mat')

Pac <-readMat('/Users/sdalin/Dropbox\ (MIT)/Biology\ PhD/2016/Hemann\ Lab/CR.CS/Creating\ Resistant\ Cells/Round\ 3/Mondays/PacPIDataFiltered.mat')

Cis <-readMat('/Users/sdalin/Dropbox\ (MIT)/Biology\ PhD/2016/Hemann\ Lab/CR.CS/Creating\ Resistant\ Cells/Round\ 3/Mondays/CisPIDataFiltered.mat')
CisData <- matrix(Cis, ncol = 8, nrow = 96)

x1 <- mtcars$mpg[mtcars$cyl==4]
x2 <- mtcars$mpg[mtcars$cyl==6]
x3 <- mtcars$mpg[mtcars$cyl==8]
vioplot(x1, x2, x3, names=c("4 cyl", "6 cyl", "8 cyl"),col="gold")
title("Violin Plots of Miles Per Gallon")

frame()
vioplot(CisData[,1],CisData[,2])

frame()
p <- ggplot(CisData, aes(factor(cyl), mpg))
p + geom_violin()
