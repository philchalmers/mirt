#
# Set of simple functions to run Mplus from R and obtain results
#

#' Reads factor scores file produced by Mplus into a data.frame
#' @param file path to Mplus factor scores file 
#' @return file read as data.frame
readMplusOutputData <- function(file){
  content <- readLines(file)
  content <- gsub('(^ +)|( +$)', '', content)
  content <- gsub(' +', '\t', content)
  writeLines(content, file)
  content <- read.table(file, sep='\t', header=F, stringsAsFactors=F)
  return(content)
}

#' Reads item parameters from a Mplus output file
#' Valid only for 1-dimensional models
#' @param file path to Mplus output file 
#' @return data.frame symilar to one returned by mirt mod2param()
readMplusOutput <- function(file){
  require(stringr)
  
  content <- readLines(file)
  tmp1 <- which(grepl("STDYX Standardization", content))
  tmp2 <- which(grepl("STDY Standardization", content))
  content <- content[(tmp1 + 6):(tmp2 - 6)]
  content <- gsub('(^ +)|( +$)', '', content)
  content <- gsub(' +', '\t', content)
  tmp <- which(grepl("Thresholds", content))
  a <- str_split(content[1:(tmp - 2)], '\t')
  a <- as.data.frame(matrix(unlist(a), nrow=length(a), ncol=5, byrow=T), stringsAsFactors=F)
  a <- a[, 1:3]
  names(a) <- c('item', 'value', 'se')
  a$name <- 'a1'
  a$item <- tolower(a$item)
  b <- str_split(content[(tmp + 1):length(content)], '\t')
  b <- as.data.frame(matrix(unlist(b), nrow=length(b), ncol=5, byrow=T), stringsAsFactors=F)
  b <- b[, 1:3]
  names(b) <- c('item', 'value', 'se')
  b$name <- sub('^.*[$]', 'd', b$item)
  b$item <- tolower(sub('[$].*$', '', b$item))
  param <- rbind(a, b)
  return(param[, c('item', 'name', 'value', 'se')])
}

#' Reads item parameters from a Mplus output file
#' Valid only for 2-dimensional models
#' @param file path to Mplus output file 
#' @return data.frame symilar to one returned by mirt mod2param()
readMplusOutput2Dim <- function(file){
  require(stringr)
  
  content <- readLines(file)
  tmp1 <- which(grepl("STDYX Standardization", content))
  tmp2 <- which(grepl("STDY Standardization", content))
  content <- content[(tmp1 + 6):(tmp2 - 7)]
  content <- gsub('(^ +)|( +$)', '', content)
  content <- gsub(' +', '\t', content)
  tmp1 <- which(grepl("THETA2\tBY", content))
  a1 <- str_split(content[1:(tmp1 - 2)], '\t')
  a1 <- as.data.frame(matrix(unlist(a1), nrow=length(a1), ncol=5, byrow=T), stringsAsFactors=F)
  a1 <- a1[, 1:3]
  names(a1) <- c('item', 'value', 'se')
  a1$name <- 'a1'
  a1$item <- tolower(a1$item)
  tmp2 <- which(grepl("Thresholds", content))
  a2 <- str_split(content[(tmp1 + 1):(tmp2 - 5)], '\t')
  a2 <- as.data.frame(matrix(unlist(a2), nrow=length(a2), ncol=5, byrow=T), stringsAsFactors=F)
  a2 <- a2[, 1:3]
  names(a2) <- c('item', 'value', 'se')
  a2$name <- 'a2'
  a2$item <- tolower(a2$item)
  cor <- c(
    'GROUP',
    unlist(str_split(content[tmp2 - 2], '\t')[1:3])[2:3],
    'COV_21'
  )
  a2[nrow(a2), ] <- cor
  b <- str_split(content[(tmp2 + 1):length(content)], '\t')
  b <- as.data.frame(matrix(unlist(b), nrow=length(b), ncol=5, byrow=T), stringsAsFactors=F)
  b <- b[, 1:3]
  names(b) <- c('item', 'value', 'se')
  b$name <- sub('^.*[$]', 'd', b$item)
  b$item <- tolower(sub('[$].*$', '', b$item))
  param <- rbind(a1, a2, b)
  return(param[, c('item', 'name', 'value', 'se')])
}
