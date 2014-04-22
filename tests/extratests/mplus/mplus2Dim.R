#
# Code estimating item parameters and factor scores for the data2Dim dataset in 
# Mplus
#

library(plyr)
library(stringr)
source('inst/extratests/mplus/mplusFunctions.R')

if(!is.null(data2Dim$theta1)){
  data2Dim <- data2Dim[, -grep('theta', names(data2Dim))]
}
dataTmp = data2Dim
dataTmp[is.na(dataTmp)] <- 99
write.table(dataTmp, 'mplus2Dim.csv', sep=',', row.names=F, col.names=F)
rm(dataTmp)

mplusSyntax <- '
DATA:
FILE IS "@DIR/mplus2Dim.csv";

VARIABLE:
NAMES ARE id g1 g2 i1-i12 j1-j40;
MISSING IS i1-j40 (99);
USEVARIABLES ARE id i1-i12 j1-j40;
CATEGORICAL ARE i1-i12 j1-j40;
IDVARIABLE IS id;

MODEL:
theta1 BY i1-i12*;
theta2 BY j1-j40*;
theta1@1;
theta2@1;

ANALYSIS:
ESTIMATOR = MLR;

OUTPUT:
STANDARDIZED;
TECH4;
TECH8;

SAVEDATA:
FILE IS "@DIR/mplus2DimOut.csv";
SAVE IS FSCORES;
'
mplusSyntax <- gsub('@DIR', getwd(), mplusSyntax)
writeLines(mplusSyntax, 'mplus2Dim.inp')
system('mplus mplus2Dim.inp')

fscores <- readMplusOutputData('mplus2DimOut.csv')
fscores <- fscores[, (ncol(fscores) - 4):ncol(fscores)]
names(fscores) <- c('id', 'theta1', 'theta1.se', 'theta2', 'theta2.se')
data2Dim <- join(data2Dim, fscores, by='id', type='left', match='first')
param2Dim <- readMplusOutput2Dim('mplus2dim.out')
for(i in 3:4){
  param2Dim[, i] = as.numeric(as.character(param2Dim[, i]))
}
save('data2Dim', 'param2Dim', file='data/data2Dim.RData')
unlink(c('mplus2Dim.inp', 'mplus2dim.out', 'mplus2Dim.csv'))
