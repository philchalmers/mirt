setwd('C:\\Users\\Phil\\Desktop\\mirt\\pkg\\mirt')
mirt <- getwd()

#initial convertion to roxygen
library(Rd2roxygen)
Rd2roxygen(mirt)

library(devtools)

load_all(mirt)
document(mirt)