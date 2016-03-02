# Violin Plots
library(R.matlab)
library(vioplot)

data <-readMat('/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/Creating Resistant Cells/Round 3/Mondays/CSV_Files/PIDataFiltered')
drugs <-c('Doxorubicin','Vincristine','Paclitaxel','Cisplatin')
(PIdata <- data.frame(data, drugs))

x1 <- mtcars$mpg[mtcars$cyl==4]
x2 <- mtcars$mpg[mtcars$cyl==6]
x3 <- mtcars$mpg[mtcars$cyl==8]
vioplot(x1, x2, x3, names=c("4 cyl", "6 cyl", "8 cyl"),
        col="gold")
title("Violin Plots of Miles Per Gallon")