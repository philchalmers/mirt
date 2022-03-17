#' Generic item summary statistics
#'
#' Function to compute generic item summary statistics that do not require
#' prior fitting of IRT models. Contains information about coefficient alpha
#' (and alpha if an item is deleted), mean/SD and frequency of total scores,
#' reduced item-total correlations, response frequencies, and conditional
#' mean/sd information given the unweighted sum scores.
#'
#' @param data An object of class \code{data.frame} or \code{matrix}
#'   with the response patterns
#' @param group optional grouping variable to condition on when computing
#'   summary information
#' @param proportions logical; include response proportion information for
#'   each item?
#' @param use_ts logical; include information that is conditional on a
#'   meaningful total score?
#' @param ts.tables logical; include mean/sd summary information
#'   pertaining to the unweighted total score?
#' @return Returns a list containing the summary statistics
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords data
#' @export
#' @examples
#'
#' # dichotomous data example
#' LSAT7full <- expand.table(LSAT7)
#' head(LSAT7full)
#' itemstats(LSAT7full)
#'
#' # behaviour with missing data
#' LSAT7full[1:5,1] <- NA
#' itemstats(LSAT7full)
#'
#' # data with no meaningful total score
#' head(SAT12)
#' itemstats(SAT12, use_ts=FALSE)
#'
#' # extra total scores tables
#' dat <- key2binary(SAT12,
#'                    key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,
#'                            5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' itemstats(dat, ts.tables=TRUE)
#'
#' # grouping information
#' group <- gl(2, 300, labels=c('G1', 'G2'))
#' itemstats(dat, group=group)
#'
#'
#' #####
#' # polytomous data example
#' itemstats(Science)
#'
#' # polytomous data with missing
#' newScience <- Science
#' newScience[1:5,1] <- NA
#' itemstats(newScience)
#'
#' # unequal categories
#' newScience[,1] <- ifelse(Science[,1] == 1, NA, Science[,1])
#' itemstats(newScience)
#'
#' merged <- data.frame(LSAT7full[1:392,], Science)
#' itemstats(merged)
#'
itemstats <- function(data, group,
                      use_ts=TRUE,
                      proportions=TRUE,
                      ts.tables=FALSE){
    data <- as.matrix(data)
    if(!missing(group) && !is.null(group)){
        groups <- unique(group)
        out <- lapply(groups, function(g){
            itemstats(data=data[group == g, , drop=FALSE], group=NULL,
                      use_ts=use_ts, proportions=proportions,
                      ts.tables=ts.tables)
        })
        names(out) <- groups
        return(out)
    }
    TS <- rowSums(data)
    if(use_ts){
        itemcor <- apply(data, 2, function(x, drop){
            tsx <- if(drop) TS-x else TS
            cor(x, tsx, use = 'complete.obs')
        }, drop=TRUE)
        itemalpha <- sapply(1:ncol(data), function(x){
            tmpdat <- na.omit(data[,-x, drop=FALSE])
            CA(tmpdat)
        })
        overall <- data.frame(N.complete=sum(!is.na(TS)), N=nrow(data),
                              mean_total.score=mean(TS, na.rm=TRUE),
                              sd_total.score=sd(TS, na.rm=TRUE),
                              alpha = CA(na.omit(data)))
        rownames(overall) <- ""
        df <- data.frame(N=apply(data, 2, function(x) sum(!is.na(x))),
                         mean=colMeans(data, na.rm = TRUE),
                         sd=apply(data, 2, sd, na.rm = TRUE),
                         item.total_cor=itemcor,
                         alpha_if_deleted=itemalpha)
    } else {
        overall <- data.frame(N.complete=sum(!is.na(TS)), N=nrow(data))
        df <- data.frame(N=apply(data, 2, function(x) sum(!is.na(x))),
                         mean=colMeans(data, na.rm = TRUE),
                         sd=apply(data, 2, sd, na.rm = TRUE))
    }

    if(overall$N == overall$N.complete) overall$N.complete <- NULL

    ret <- list(overall=as.mirt_df(overall),
                itemstats=as.mirt_df(df))
    if(proportions){
        useNA <- ifelse(any(is.na(data)), 'always', 'ifany')
        props <- apply(data, 2, function(x){
            out <- prop.table(table(x, useNA=useNA))
            if(useNA == 'always' && out[length(out)] == 0)
                out[length(out)] <- NA
            out
        })
        if(is.list(props)){
            vals <- sort(unique(as.numeric(data)), na.last=TRUE)
            proportions <- matrix(NA, ncol(data), length(vals))
            colnames(proportions) <- vals
            rownames(proportions) <- colnames(data)
            for(i in 1:ncol(data)){
                pick <- names(props[[i]])
                proportions[i,colnames(proportions) %in% pick] <- props[[i]]
            }
            ret$proportions <- as.mirt_df(as.data.frame(proportions))
        } else {
            ret$proportions <- as.mirt_df(as.data.frame(t(props)))
        }
    }
    if(ts.tables){
        ret$total.score_frequency <- as.data.frame(t(as.matrix(table(TS))))
        rownames(ret$total.score_frequency) <- "Freq"
        ret$total.score_means <- t(apply(data, 2, function(x){
            tapply(TS, x, mean, na.rm=TRUE)
        }))
        ret$total.score_sds <- t(apply(data, 2, function(x){
            tapply(TS, x, sd, na.rm=TRUE)
        }))
    }
    ret
}

