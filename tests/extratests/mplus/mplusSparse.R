#
# Code estimating item parameters and factor scores for the dataSparse dataset
# in Mplus
#

library(plyr)
library(stringr)
source('inst/extratests/mplus/mplusFunctions.R')

if(!is.null(dataSparse$theta)){
  dataSparse <- dataSparse[, -grep('theta', names(dataSparse))]
}
dataTmp = dataSparse
dataTmp[is.na()] <- 99
write.table(dataTmp, 'mplusSparse.csv', sep=',', row.names=F, col.names=F)
rm(dataTmp)

mplusSyntax <- '
DATA:
FILE IS "@DIR/mplusSparse.csv";

VARIABLE:
NAMES ARE id g1 i1-i142;
MISSING IS i1-i142 (99);
USEVARIABLES ARE id i1-i142;
CATEGORICAL ARE i1-i142;
IDVARIABLE IS id;

MODEL:
theta BY i1-i142*;
theta@1;

ANALYSIS:
ESTIMATOR = MLR;

OUTPUT:
STANDARDIZED;
TECH4;
TECH8;

SAVEDATA:
FILE IS "@DIR/mplusSparseOut.csv";
SAVE IS FSCORES;
'
mplusSyntax <- gsub('@DIR', getwd(), mplusSyntax)
writeLines(mplusSyntax, 'mplusSparse.inp')
system('mplus mplusSparse.inp')

fscores <- readMplusOutputData('mplusSparseOut.csv')
fscores <- fscores[, (ncol(fscores) - 2):ncol(fscores)]
names(fscores) <- c('id', 'theta', 'theta.se')
dataSparse <- join(dataSparse, fscores, by='id', type='left', match='first')
paramSparse <- readMplusOutput('mplussparse.out')
for(i in 3:4){
  paramSparse[, i] = as.numeric(as.character(paramSparse[, i]))
}
save('dataSparse', 'paramSparse', file='data/dataSparse.RData')
unlink(c('mplusSparse.inp', 'mplussparse.out', 'mplusSparseOut.csv'))

