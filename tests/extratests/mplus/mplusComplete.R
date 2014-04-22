#
# Code estimating item parameters and factor scores for the dataComplete 
# dataset in Mplus
#

library(plyr)
library(stringr)
source('inst/extratests/mplus/mplusFunctions.R')

if(!is.null(dataComplete$theta)){
  dataComplete <- dataComplete[, -grep('theta', names(dataComplete))]
}
dataTmp = dataComplete
dataTmp[is.na(dataTmp)] <- 99
write.table(dataTmp, 'mplusComplete.csv', sep=',', row.names=F, col.names=F)
rm(dataTmp)

mplusSyntax <- '
DATA:
  FILE IS "@DIR/mplusComplete.csv";

VARIABLE:
  NAMES ARE id g1 g2 g3 i1-i27;
  MISSING IS i1-i27 (99);
  USEVARIABLES ARE id i1-i27;
  CATEGORICAL ARE i1-i27;
  IDVARIABLE IS id;

MODEL:
  theta BY i1-i27*;
  theta@1;

ANALYSIS:
  ESTIMATOR = MLR;

OUTPUT:
  STANDARDIZED;
  TECH4;
  TECH8;

SAVEDATA:
  FILE IS "@DIR/mplusCompleteOut.csv";
  SAVE IS FSCORES;
'
mplusSyntax <- gsub('@DIR', getwd(), mplusSyntax)
writeLines(mplusSyntax, 'mplusComplete.inp')
system('mplus mplusComplete.inp')

fscores <- readMplusOutputData('mplusCompleteOut.csv')
fscores <- fscores[, (ncol(fscores) - 2):ncol(fscores)]
names(fscores) <- c('id', 'theta', 'theta.se')
dataComplete <- join(dataComplete, fscores, by='id', type='left', match='first')
paramComplete <- readMplusOutput('mpluscomplete.out')
for(i in 3:4){
  paramComplete[, i] = as.numeric(as.character(paramComplete[, i]))
}
save('dataComplete', 'paramComplete', file='data/dataComplete.RData')
unlink(c('mplusComplete.inp', 'mpluscomplete.out', 'mplusCompleteOut.csv'))

