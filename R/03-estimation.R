ESTIMATION <- function(data, model, group, itemtype = NULL, guess = 0, upper = 1,
                       invariance = '', pars = NULL, constrain = NULL, key = NULL,
                       parprior = NULL, mixed.design = NULL, customItems = NULL,
                       customItemsData = NULL, customGroup = NULL,
                       GenRandomPars = FALSE, large = FALSE,
                       survey.weights = NULL, latent.regression = NULL,
                       gpcm_mats=list(), spline_args=list(), monopoly.k=1,
                       control = list(), ...)
{
    start.time <- proc.time()[3L]
    dots <- list(...)
    if(!is.null(itemtype))
        itemtype <- ifelse(itemtype == 'grsm', 'grsmIRT', itemtype)
    if(missing(data)) missingMsg('data')
    if(length(unique(colnames(data))) != ncol(data))
        stop('items must have unique names in data input', call.=FALSE)
    if(missing(model)) missingMsg('model')
    if(missing(group)) missingMsg('group')
    if(!(is.factor(group) || is.character(group)) || length(group) != nrow(data))
        stop('group input provided is not valid', call.=FALSE)
    if(is.character(large) && large == 'return'){
        Data <- opts <- list()
        opts$zeroExtreme <- FALSE
        opts$dentype <- 'default'
        if(!is.null(dots$technical$zeroExtreme)) opts$zeroExtreme <- dots$technical$zeroExtreme
        data <- as.matrix(data)
        if(is.numeric(data))
            data <- matrix(as.integer(data), nrow(data), dimnames=list(rownames(data), colnames(data)))
        rownames(data) <- 1L:nrow(data)
        if(is.null(colnames(data)))
            colnames(data) <- paste0('Item.', 1L:ncol(data))
        Data$data <- data
        Data$group <- factor(group)
        Data$groupNames <- levels(Data$group)
        Data$ngroups <- length(Data$groupNames)
        Data$nitems <- ncol(data)
        K <- apply(Data$data, 2L, function(x) length(unique(na.omit(x))))
        tmp <- list(...)
        if(!is.null(tmp$technical$customK)) K <- tmp$technical$customK
        itemloc <- c(1L, cumsum(K) + 1L)
        Names <- NULL
        for(i in 1L:length(K))
            Names <- c(Names, paste("Item.",i,"_",1L:K[i],sep=""))
        PrepListFull <- list(K=K, itemloc=itemloc,
                             Names=Names, itemnames=colnames(Data$data))
    } else {
        if(missing(data) || is.null(nrow(data)))
            stop('data argument is required', call.=FALSE)
        if(missing(model) || !is(model, 'numeric') && !is(model, 'mirt.model') &&
           !is(model, 'character'))
            stop('model argument (numeric, character, or from mirt.model() function) is required',
                 call.=FALSE)
        if(!(is.factor(group) || is.character(group)) || length(group) != nrow(data))
            stop('group input provided is not valid', call.=FALSE)
        if(is.matrix(itemtype)){
            itemtypefull <- itemtype
            itemtype <- itemtypefull[1L,]
        } else itemtypefull <- NULL
        if(!is.null(itemtype))
            stopifnot(is(itemtype, 'character'))
        if(!is.null(constrain))
            stopifnot(is(constrain, 'list'))
        if(!is.null(parprior))
            stopifnot(is(parprior, 'list'))
        if(!is.null(customItems))
            stopifnot(is(customItems, 'list'))
        if(!is.null(customItemsData))
            stopifnot(is(customItemsData, 'list') || length(customItemsData) == ncol(data))
        if(!is.null(customGroup))
            stopifnot(is(customGroup, 'GroupPars') || is.list(customGroup))
        stopifnot(is(invariance, 'character'))
        stopifnot(is(GenRandomPars, 'logical'))
        stopifnot(is(large, 'logical') || is(large, 'list') || is(large, 'character'))
        opts <- makeopts(GenRandomPars=GenRandomPars,
                         hasCustomGroup=!is.null(customGroup), ...)
        if(opts$Moptim == 'NR'){
            if(is.null(control$tol)) control$tol <- opts$TOL/1000
            if(is.null(control$maxit)) control$maxit <- 50L
        }
        if(opts$Moptim == 'nlminb')
            if(is.null(control$rel.tol)) control$rel.tol <- opts$TOL/100

        if(opts$dentype == "discrete" && is.null(customGroup)){
            den <- Theta_discrete_den
            par <- if(!is.null(dots$structure)){
                if(is.null(opts$technical$customTheta))
                    stop('customTheta must be defined when using structure input', call.=FALSE)
                opts$structure <- model.matrix(dots$structure, as.data.frame(opts$technical$customTheta))
                tmpnfact <- ncol(opts$technical$customTheta)
                numeric(ncol(opts$structure))
            } else if(is.null(opts$technical$customTheta)){
                tmpnfact <- model
                numeric(model-1L)
            } else {
                tmpnfact <- ncol(opts$technical$customTheta)
                numeric(nrow(opts$technical$customTheta) - 1L)
            }
            if(length(par)){
                names(par) <- paste0('c', 1:length(par))
                est <- rep(TRUE, length(par))
            } else {
                par <- c(c = 0)
                est <- FALSE
            }
            customGroup <- createGroup(par=par, est=est, den=den, nfact=tmpnfact,
                                       gen=function(object) rnorm(length(object@par), 0, 1/2))
            customGroup@itemclass <- -1L
            rm(tmpnfact)
        }
        if(!is.null(survey.weights)){
            stopifnot(opts$method %in% c('EM', 'QMCEM', 'MCEM'))
            stopifnot(length(survey.weights) == nrow(data))
        }
        if(any(is.na(group))){
            if(opts$message)
                message('NA values in group removed, along with associated rows in data')
            data <- data[!is.na(group), , drop=FALSE]
            group <- group[!is.na(group)]
        }
        if(!is.null(customItems) || any(itemtype %in% Experimental_itemtypes()))
            opts$calcNull <- FALSE
        opts$times <- list(start.time=start.time)
        # on exit, reset the seed to override internal
        if(opts$method == 'MHRM' || opts$method == 'MIXED' || opts$method == 'SEM' &&
           opts$plausible.draws == 0L)
            on.exit(set.seed((as.numeric(Sys.time()) - floor(as.numeric(Sys.time()))) * 1e8))
        #change itemtypes if NULL.MODEL
        if(opts$NULL.MODEL){
            constrain <- NULL
            if(!is.null(itemtype)){
                itemtype[itemtype == 'grsm'] <- 'graded'
                itemtype[itemtype == 'rsm'] <- 'gpcm'
                itemtype[itemtype == '3PL' | itemtype == '3PLu' | itemtype == '4PL'] <- '2PL'
                itemtype[itemtype == '3PLNRM' | itemtype == '3PLuNRM' | itemtype == '4PLNRM'] <- '2PLNRM'
                itemtype[itemtype == 'spline'] <- '2PL'
            }
        }
        if(!is.null(itemtype)){
            if(any(itemtype == 'spline') && !(opts$method %in% c('EM', 'QMCEM', 'MCEM')))
                stop('spline itemtype only supported for EM algorithm', call.=FALSE)
        }
        if(length(group) != nrow(data))
            stop('length of group not equal to number of rows in data.', call.=FALSE)
        if(any(is.na(group)))
            stop('Unknown group memberships cannot be estimated. Please remove the NA values in group
                    and the associated rows in the data input.', call.=FALSE)

        #####
        ### uncomment when building html examples
        # opts$verbose <- FALSE
        ####
        opts$times$start.time.Data <- proc.time()[3L]
        Data <- list()
        data <- as.matrix(data)
        if(!all(apply(data, 2L, class) %in% c('integer', 'numeric')))
            stop('Input data must be integer or numeric values only', call.=FALSE)
        if(is.numeric(data))
            data <- matrix(as.integer(data), nrow(data),
                           dimnames=list(rownames(data), colnames(data)))
        rownames(data) <- 1L:nrow(data)
        if(is.null(colnames(data)))
            colnames(data) <- paste0('Item.', 1L:ncol(data))
        if(nrow(data) > 1L && is.null(opts$technical$customK)){
            if(!is.null(key)){
                key <- sapply(1L:length(key), function(ind, data, key){
                    if(is.na(key[ind])) return(NA)
                    s <- sort(unique(data[,ind]))
                    se <- min(s, na.rm = TRUE):(length(s) + min(s) - 1L)
                    ret <- se[s == key[ind]]
                    if(!length(ret)) ret <- NA
                    ret
                }, data=data, key=key)
            }
            data <- remap.distance(data, message = opts$message)
        }
        Data$rowID <- 1L:nrow(data)
        Data$completely_missing <- which(rowSums(is.na(data)) == ncol(data))
        if(length(Data$completely_missing)){
            if(opts$warn)
                warning('data contains response patterns with only NAs', call.=FALSE)
            pick <- rowSums(is.na(data)) != ncol(data)
            Data$rowID <- unname(which(pick))
            data <- subset(data, pick)
            group <- subset(group, pick)
            if(!is.null(latent.regression) || !is.null(mixed.design))
                covdata <- covdata[-Data$completely_missing, , drop=FALSE]
        }
        Data$data <- data

        if(is.null(opts$grsm.block)) Data$grsm.block <- rep(1L, ncol(data))
        else Data$grsm.block <- opts$grsm.block
        if(is.null(opts$rsm.block)) Data$rsm.block <- rep(1L, ncol(data))
        else Data$rsm.block <- opts$rsm.block
        Data$group <- factor(group)
        Data$groupNames <- levels(Data$group)
        if(any(grepl('-', Data$groupNames)))
            stop('Group names cannot contain a dash (-) character', call.=FALSE)
        Data$ngroups <- if(!is.null(opts$ngroups)) opts$ngroups else length(Data$groupNames)
        if(!is.null(itemtypefull))
            if(Data$ngroups != nrow(itemtypefull))
                stop("number of groups does not match number of rows in itemtype", call.=FALSE)
        Data$nitems <- ncol(Data$data)
        Data$N <- nrow(Data$data)
        Data$mins <- suppressWarnings(apply(data, 2L, min, na.rm=TRUE))
        Data$mins[!is.finite(Data$mins)] <- 0L
        if(!is.null(opts$technical$mins)) Data$mins <- opts$technical$mins
        if(is.character(model)){
            tmp <- any(sapply(colnames(data), grepl, x=model))
            model <- mirt.model(model, itemnames = if(tmp) colnames(data) else NULL)
        }
        oldmodel <- model
        PrepList <- vector('list', Data$ngroups)
        tmp <- 1L:Data$ngroups
        if(opts$dentype == "mixture")
            Data$groupNames <- paste0("MIXTURE_", seq_len(Data$ngroups))
        names(PrepList) <- Data$groupNames
        model <- buildModelSyntax(model, J=Data$nitems, groupNames=Data$groupNames,
                                  itemtype=itemtype)
        Data$model <- model
        if(!is.null(customGroup)){
            if(Data$ngroups == 1L) customGroup <- list(customGroup)
            else {
                if(Data$ngroups != length(customGroup) && length(customGroup) == 1L){
                    tmplst <- vector('list', Data$ngroups)
                    for(g in seq_len(Data$ngroups))
                        tmplst[[g]] <- customGroup
                    customGroup <- tmplst
                    rm(tmplst)
                }
            }
        }
        if(!is.null(dots$PrepList)) {
            PrepListFull <- PrepList[[1L]] <- dots$PrepList
        } else {
            PrepListFull <- PrepList[[1L]] <-
                PrepData(data=Data$data, model=Data$model, itemtype=itemtype, guess=guess,
                         upper=upper, parprior=parprior, verbose=opts$verbose,
                         technical=opts$technical, parnumber=1L, BFACTOR=opts$dentype == 'bfactor',
                         grsm.block=Data$grsm.block, rsm.block=Data$rsm.block,
                         mixed.design=mixed.design, customItems=customItems,
                         customItemsData=customItemsData,
                         customGroup=customGroup[[1L]], spline_args=spline_args, monopoly.k=monopoly.k,
                         fulldata=opts$PrepList[[1L]]$fulldata, key=key, opts=opts,
                         gpcm_mats=gpcm_mats, internal_constraints=opts$internal_constraints,
                         dcIRT_nphi=opts$dcIRT_nphi, dentype=opts$dentype, item.Q=opts$item.Q)
            if(!is.null(dots$Return_PrepList)) return(PrepListFull)
            if(!is.null(itemtypefull)){
                for(g in 2L:nrow(itemtypefull)){
                    PrepList[[g]] <-
                        PrepData(data=Data$data, model=Data$model, itemtype=itemtypefull[g,], guess=guess,
                                 upper=upper, parprior=parprior, verbose=opts$verbose,
                                 technical=opts$technical, parnumber=1L, BFACTOR=opts$dentype == 'bfactor',
                                 grsm.block=Data$grsm.block, rsm.block=Data$rsm.block,
                                 mixed.design=mixed.design, customItems=customItems,
                                 customItemsData=customItemsData,
                                 customGroup=customGroup[[1L]], spline_args=spline_args, monopoly.k=monopoly.k,
                                 fulldata=opts$PrepList[[1L]]$fulldata, key=key, opts=opts,
                                 gpcm_mats=gpcm_mats, internal_constraints=opts$internal_constraints,
                                 dcIRT_nphi=opts$dcIRT_nphi, dentype=opts$dentype, item.Q=opts$item.Q)
                }
            }
        }
        if(!is.null(opts$structure)){
            PrepList[[1L]]$pars[[Data$nitems+1L]]@structure <- opts$structure
            PrepList[[1L]]$pars[[Data$nitems+1L]]@parnames <-
                names(PrepList[[1L]]$pars[[Data$nitems+1L]]@est) <- colnames(opts$structure)
            PrepList[[1L]]$pars[[Data$nitems+1L]]@itemclass <- -999L
        }
        parnumber <- max(PrepList[[1L]]$pars[[Data$nitems+1L]]@parnum) + 1L
        attr(PrepListFull$pars, 'nclasspars') <- attr(PrepList[[1L]]$pars, 'nclasspars') <-
            matrix(sapply(PrepListFull$pars, function(y) length(y@parnum)), nrow=1L)
        for(g in seq_len(Data$ngroups)){
            if(g != 1L){
                if(is.null(itemtypefull))
                    PrepList[[g]] <- list(pars=PrepList[[1L]]$pars)
                else attr(PrepList[[g]]$pars, 'nclasspars') <-
                        sapply(PrepList[[g]]$pars, function(y) length(y@parnum))
                for(i in 1L:length(PrepList[[g]]$pars)){
                    PrepList[[g]]$pars[[i]]@parnum <- parnumber:(parnumber +
                                                    length(PrepList[[g]]$pars[[i]]@parnum) - 1L)
                    parnumber <- max(PrepList[[g]]$pars[[i]]@parnum) + 1L
                }
                if(!is.null(customGroup)){
                    customGroup[[g]]@parnum <- PrepList[[g]]$pars[[Data$nitems + 1L]]@parnum
                    PrepList[[g]]$pars[[Data$nitems + 1L]] <- customGroup[[g]]
                }
            }
        }
        if(length(mixed.design$random) > 0L){
            for(i in seq_len(length(mixed.design$random))){
                mixed.design$random[[i]]@parnum <- parnumber:(parnumber - 1L +
                                                                  length(mixed.design$random[[i]]@par))
                parnumber <- max(mixed.design$random[[i]]@parnum) + 1L
            }
        }
        if(!is.null(latent.regression)){
            if(length(PrepListFull$prodlist))
                stop('Polynomial combinations currently not supported when latent regression effects are used', call.=FALSE)
            lrPars <- make.lrdesign(df=latent.regression$df, formula=latent.regression$formula,
                                    factorNames=PrepListFull$factorNames, EM=latent.regression$EM,
                                    TOL=opts$TOL)
            lrPars@parnum <- parnumber:(parnumber - 1L + length(lrPars@par))
            parnumber <- max(lrPars@parnum) + 1L
            if(opts$dentype == 'discrete'){
                tmp <- matrix(1L:length(lrPars@beta), nrow(lrPars@beta), ncol(lrPars@beta))
                tmp2 <- tmp[1, ]
                lrPars@est[tmp2[-length(tmp2)]] <- TRUE
                lrPars@est[tmp[,ncol(tmp)]] <- FALSE
                PrepList$all$pars[[ncol(data) + 1L]]@est <- FALSE
            }
        } else lrPars <- list()
        if(length(latent.regression$lr.random) > 0L){
            for(i in seq_len(length(latent.regression$lr.random))){
                latent.regression$lr.random[[i]]@parnum <- parnumber:(parnumber - 1L +
                                                                  length(latent.regression$lr.random[[i]]@par))
                parnumber <- max(latent.regression$lr.random[[i]]@parnum) + 1L
            }
        }
    }
    if(!is.null(opts$PrepList)){
        PrepList <- opts$PrepList
    } else {
        if(!is.list(large)){
            tmpdata <- Data$data
            if(opts$zeroExtreme){
                n_is_na <- rowSums(is.na(tmpdata))
                sums <- colSums(t(tmpdata) - Data$mins, na.rm=TRUE)
                maxs <- apply(tmpdata, 2L, function(x) length(unique(na.omit(x)))) - 1L
                if(is.null(survey.weights)) survey.weights <- rep(1, nrow(tmpdata))
                survey.weights[sums == 0L | (sums + n_is_na) == sum(maxs)] <- 0
            }
            tmptabdata <- if(is.logical(large) && large){
                maketabDataLarge(tmpdata=tmpdata, group=Data$group,
                                 groupNames=if(opts$dentype != 'mixture') Data$groupNames else 'full',
                                 nitem=Data$nitems, K=PrepListFull$K, itemloc=PrepListFull$itemloc,
                                 Names=PrepListFull$Names, itemnames=PrepListFull$itemnames,
                                 survey.weights=survey.weights)
            } else {
                maketabData(tmpdata=tmpdata, group=Data$group,
                            groupNames=if(opts$dentype != 'mixture') Data$groupNames else 'full',
                            nitem=Data$nitems, K=PrepListFull$K, itemloc=PrepListFull$itemloc,
                            Names=PrepListFull$Names, itemnames=PrepListFull$itemnames,
                            survey.weights=survey.weights)
            }
            if(is.character(large) && large == 'return'){
                return(tmptabdata)
            } else large <- tmptabdata
            rm(tmpdata, tmptabdata)
        }
        Data$tabdatalong <- large$tabdata
        Data$tabdata <- large$tabdata2
        for(g in seq_len(Data$ngroups)){
            if(opts$dentype == 'mixture'){
                Data$fulldata[[g]] <- PrepListFull$fulldata
                Data$Freq[[g]] <- large$Freq[[1L]]
            } else {
                Data$fulldata[[g]] <- PrepListFull$fulldata[Data$group == Data$groupNames[g],
                                                            , drop=FALSE]
                Data$Freq[[g]] <- large$Freq[[g]]
            }
        }
    }
    if(opts$returnPrepList) return(PrepList)
    if(opts$dentype == 'bfactor'){
        #better start values
        if((PrepList[[1L]]$nfact - attr(model, 'nspec')) == 1L){
            nfact <- PrepListFull$nfact
            for(g in seq_len(Data$ngroups)){
                for(i in seq_len(Data$nitems)){
                    tmp <- PrepList[[g]]$pars[[i]]@par[1L:nfact]
                    tmp2 <- tmp[1L]
                    tmp[PrepList[[g]]$pars[[i]]@est[1L:nfact]] <-
                        tmp[PrepList[[g]]$pars[[i]]@est[1L:nfact]]/2
                    tmp[1L] <- tmp2
                    PrepList[[g]]$pars[[i]]@par[1L:nfact] <- tmp
                }
            }
        }
    }
    PrepList <- UpdateParameters(PrepList, model, groupNames=Data$groupNames)
    if(GenRandomPars){
        for(g in seq_len(Data$ngroups))
            for(i in seq_len(length(PrepList[[g]]$pars)))
                PrepList[[g]]$pars[[i]] <- GenRandomPars(PrepList[[g]]$pars[[i]])

    }
    if(opts$dentype == "discrete"){
        PrepList[[1L]]$exploratory <- FALSE
        if(is.null(opts$technical$customTheta))
            opts$technical$customTheta <- diag(PrepList[[1L]]$nfact)
    }
    RETURNVALUES <- FALSE
    if(!is.null(pars)){
        if(is(pars, 'data.frame')){
            PrepList <- UpdatePrepList(PrepList, pars, random=mixed.design$random,
                                       lrPars=lrPars, lr.random=latent.regression$lr.random,
                                       MG = TRUE)
            mixed.design$random <- attr(PrepList, 'random')
            if(any(pars$class == 'lrPars')) lrPars <- update.lrPars(pars, lrPars)
            attr(PrepList, 'random') <- NULL
        }
        if(!is.null(attr(pars, 'values')) || (is.character(pars) && pars == 'values'))
            RETURNVALUES <- TRUE
    }
    pars <- vector('list', Data$ngroups)
    for(g in seq_len(Data$ngroups))
        pars[[g]] <- PrepList[[g]]$pars
    nitems <- Data$nitems
    Data$K <- PrepList[[1L]]$K
    nfact <- PrepList[[1L]]$pars[[nitems+1L]]@nfact
    if(nfact != 1L && any(c('Rasch') %in% itemtype ) && PrepList[[1L]]$exploratory)
       stop('Rasch itemtype is for confirmatory models only.', call.=FALSE)
    nLambdas <- PrepList[[1L]]$pars[[1L]]@nfact
    if(is.null(constrain)) constrain <- list()
    nspec <- ifelse(!is.null(attr(model, 'nspec')), attr(model, 'nspec'), 1L)
    #default MG uses configural model (independent groups but each identified)
    if('free_means' %in% invariance ){ #Free factor means (means 0 for ref)
        if(opts$dentype == 'bfactor'){
            for(g in 2L:Data$ngroups)
                pars[[g]][[nitems + 1L]]@est[1L:(nfact-nspec)] <- TRUE
        } else {
            for(g in 2L:Data$ngroups)
                pars[[g]][[nitems + 1L]]@est[1L:nfact] <- TRUE
        }
    }
    dummymat <- matrix(FALSE, pars[[1L]][[nitems + 1L]]@nfact, pars[[1L]][[nitems + 1L]]@nfact)
    if(any('free_var' %in% invariance)){ #Free factor vars (vars 1 for ref)
        if(opts$dentype == 'bfactor'){
            tmp <- dummymat[1L:(nfact-nspec),1L:(nfact-nspec), drop=FALSE]
            diag(tmp) <- TRUE
            dummymat[1L:(nfact-nspec),1L:(nfact-nspec)] <- tmp
        } else {
            diag(dummymat) <- TRUE
        }
        tmp <- dummymat[lower.tri(dummymat, TRUE)]
        if(opts$dentype %in% c('Davidian', 'EHW')){
            for(g in 2L:Data$ngroups)
                pars[[g]][[nitems + 1L]]@est[2L] <- TRUE
        } else {
            for(g in 2L:Data$ngroups){
                pars[[g]][[nitems + 1L]]@est <- pars[[g]][[nitems + 1L]]@est |
                    c(pars[[g]][[nitems + 1L]]@est[1L:pars[[g]][[nitems + 1L]]@nfact], tmp)
                names(pars[[g]][[nitems + 1L]]@est) <- names(pars[[g]][[nitems + 1L]]@par)
                pars[[g]][[nitems + 1L]]@parnames <- names(pars[[g]][[nitems + 1L]]@est)
            }
        }
    }
    if(opts$dentype == 'mixture'){
        tmp <- length(pars[[1L]][[nitems + 1L]]@par)
        pars[[1L]][[nitems + 1L]]@est[tmp] <- FALSE
        for(g in 1L:Data$ngroups)
            pars[[g]][[nitems + 1L]]@par[tmp] <- g - 1
        names(PrepList) <- Data$groupNames
    }
    if(RETURNVALUES){
        for(g in seq_len(Data$ngroups))
            PrepList[[g]]$pars <- pars[[g]]
        return(ReturnPars(PrepList, PrepList[[1L]]$itemnames, lr.random=latent.regression$lr.random,
                          random=mixed.design$random, lrPars=lrPars, MG = TRUE))
    }
    constrain <- UpdateConstrain(pars=pars, constrain=constrain, invariance=invariance, nfact=Data$nfact,
                                 nLambdas=nLambdas, J=nitems, ngroups=Data$ngroups, PrepList=PrepList,
                                 method=opts$method, itemnames=PrepList[[1L]]$itemnames, model=model,
                                 groupNames=Data$groupNames)
    startlongpars <- c()
    if(opts$NULL.MODEL){
        constrain <- list()
        for(g in seq_len(length(pars))){
            for(i in seq_len(nitems)){
                pars[[g]][[i]] <- set_null_model(pars[[g]][[i]])
            }
        }
    }
    rr <- 0L
    for(g in seq_len(Data$ngroups)){
        if(g > 1L && opts$dentype == 'mixture') break
        r <- Data$Freq[[g]]
        rr <- rr + r
    }
    df <- Data$ngroups * (prod(PrepList[[1L]]$K) - 1)
    if(df > 1e10) df <- 1e10
    nestpars <- nconstr <- 0L
    for(g in seq_len(Data$ngroups))
        for(i in seq_len(nitems+1L))
            nestpars <- nestpars + sum(pars[[g]][[i]]@est)
    if(!is.null(mixed.design$random)){
        for(i in seq_len(length(mixed.design$random)))
            nestpars <- nestpars + sum(mixed.design$random[[i]]@est)
    }
    if(length(lrPars))
        nestpars <- nestpars + sum(lrPars@est)
    if(!is.null(dots$structure)) nestpars <- nestpars - 1L
    if(!is.null(opts$technical$nconstrain))
        constrain <- c(constrain, opts$technical$nconstrain)
    if(length(constrain) > 0L)
        for(i in seq_len(length(constrain)))
            nconstr <- nconstr + length(constrain[[i]]) - 1L
    if(Data$ngroups > 1L && !length(constrain)){
        if(opts$warn && any(invariance %in% c('free_means', 'free_var')))
            warning('Multiple-group model may not be identified without providing anchor items',
                    call.=FALSE)
        for(j in seq_len(Data$ngroups)){
            if(opts$dentype != 'mixture')
                tmp <- apply(subset(Data$data, Data$group == Data$groupNames[j]), 2L,
                             function(x) length(unique(na.omit(x)))) == Data$K
            for(i in which(!tmp)){
                if(any(PrepList[[j]]$pars[[i]]@est))
                    stop(paste0('Multiple Group model will not be identified without ',
                                'proper constraints (groups contain missing data patterns ',
                                'where item responses have been completely omitted or, alternatively, ',
                                'the number of categories within each group is not equal to the ',
                                'total number of categories)'),
                         call. = FALSE)
            }
        }
    }
    nmissingtabdata <- sum(is.na(rowSums(Data$tabdata)))
    dfsubtr <- nestpars - nconstr
    if(opts$dentype == 'EH') dfsubtr <- dfsubtr + (opts$quadpts - 1L) * Data$ngroups
    else if(opts$dentype == 'EHW') dfsubtr <- dfsubtr + (opts$quadpts - 3L) * Data$ngroups
    if(df <= dfsubtr)
        stop('Too few degrees of freedom. There are only ', df, ' degrees of freedom but ',
             dfsubtr, ' parameters were freely estimated.', call.=FALSE)
    df <- df - dfsubtr
    if(!is.null(customItems)){
        for(g in seq_len(Data$ngroups))
            PrepList[[g]]$exploratory <- FALSE
    }
    G2group <- numeric(Data$ngroups)
    DERIV <- vector('list', Data$ngroups)
    for(g in seq_len(Data$ngroups)){
        DERIV[[g]] <- vector('list', Data$nitems)
        for(i in seq_len(Data$nitems))
            DERIV[[g]][[i]] <- selectMethod(Deriv, c(class(pars[[g]][[i]]), 'matrix'))
    }
    Ls <- makeLmats(pars, constrain, random = mixed.design$random,
                    lr.random=latent.regression$lr.random, lrPars=lrPars,
                    nconstrain=opts$technical$nconstrain)
    CUSTOM.IND <- which(sapply(pars[[1L]], class) %in% Use_R_ProbTrace())
    SLOW.IND <- which(sapply(pars[[1L]], class) %in% Use_R_Deriv())
    if(pars[[1]][[length(pars[[1L]])]]@itemclass %in% c(-1L, -999L))
        SLOW.IND <- c(SLOW.IND, length(pars[[1L]]))
    if(opts$dentype != 'Gaussian' && opts$method %in% c('MHRM', 'MIXED', 'SEM'))
        stop('Non-Gaussian densities not currently supported with MHRM algorithm')
    #warnings
    wmsg <- 'Lower and upper bound parameters (g and u) should use \'norm\' (i.e., logit) prior'
    for(g in seq_len(length(pars))){
        for(i in seq_len(length(pars[[1L]]))){
            if(class(pars[[g]][[i]]) == 'dich'){
                pt <- pars[[g]][[i]]@prior.type
                if(!(pt[length(pt)-1L] %in% c(0L, 1L, 4L))) warning(wmsg, call.=FALSE)
                if(!(pt[length(pt)] %in% c(0L, 1L, 4L))) warning(wmsg, call.=FALSE)
                next
            } else if(class(pars[[g]][[i]]) == 'partcomp'){
                pt <- pars[[g]][[i]]@prior.type
                if(!(pt[length(pt)-1L] %in% c(0L, 1L))) warning(wmsg, call.=FALSE)
                if(!(pt[length(pt)] %in% c(0L, 1L))) warning(wmsg, call.=FALSE)
                next
            } else if(class(pars[[g]][[i]]) == 'nestlogit'){
                pt <- pars[[g]][[i]]@prior.type
                if(!(pt[nfact + 2L] %in% c(0L, 1L))) warning(wmsg, call.=FALSE)
                if(!(pt[nfact + 3L] %in% c(0L, 1L))) warning(wmsg, call.=FALSE)
                next
            }
        }
    }
    SEMconv <- NA
    opts$times$end.time.Data <- proc.time()[3L]

    #EM estimation
    opts$times$start.time.Estimate <- proc.time()[3L]
    logLik <- G2 <- SElogLik <- NaN
    logPrior <- 0
    Pl <- vector('list', Data$ngroups)
    if(opts$method %in% c('EM', 'BL', 'QMCEM', 'MCEM')){
        logLik <- G2 <- SElogLik <- 0
        if(length(lrPars)){
            if(opts$SE && !(opts$SE.type %in% c('complete', 'forward', 'central', 'Richardson')))
                stop('Information matrix method for latent regression estimates not supported',
                     call.=FALSE)
            opts$full <- TRUE
        } else opts$full <- FALSE
        temp <- matrix(0L,nrow=nitems,ncol=nspec)
        sitems <- matrix(0L, nrow=sum(PrepList[[1L]]$K), ncol=nspec)
        specific <- NULL
        if(opts$dentype == "discrete"){
            theta <- 0
            Theta <- opts$technical$customTheta
            opts$quadpts <- nrow(Theta)
            if(pars[[1L]][[1L]]@nfact != ncol(Theta))
                stop("mirt.model definition does not have same number of traits/attributes as customTheta input", call.=FALSE)
        } else {
            if(is.null(opts$quadpts)){
                tmp <- if(opts$dentype == 'bfactor') PrepList[[1L]]$nfact - attr(model, 'nspec') + 1L
                    else nfact
                opts$quadpts <- select_quadpts(tmp)
            }
            if(opts$quadpts < 3 && opts$warn) warning('Should use more than 2 quadpts', call.=FALSE)
            theta <- 1
            if(!(opts$method %in% c('QMCEM', 'MCEM')))
                theta <- as.matrix(seq(opts$theta_lim[1L], opts$theta_lim[2L],
                                       length.out = opts$quadpts))
            if(opts$dentype == 'bfactor'){
                specific <- attr(oldmodel, 'specific')
                specific[is.na(specific)] <- 1L
                for(i in seq_len(nitems)) temp[i, specific[i]] <- 1L
                ind <- 1L
                for(i in seq_len(nitems)){
                    for(j in seq_len(PrepList[[1L]]$K[i])){
                        sitems[ind, ] <- temp[i, ]
                        ind <- ind + 1L
                    }
                }
                nfact2 <- PrepList[[1L]]$nfact - attr(model, 'nspec') + 1L
                Theta <- thetaComb(theta, nfact2)
                Theta <- cbind(Theta[,1L:(nfact2-1L),drop=FALSE],
                               matrix(Theta[,nfact2], nrow=nrow(Theta), ncol=ncol(sitems)))
            } else {
                if(opts$method %in% c('QMCEM', 'MCEM')){
                    Theta <- NULL
                } else {
                    if(opts$quadpts^nfact <= opts$MAXQUAD){
                        if(is.null(opts$technical$customTheta))
                            Theta <- thetaComb(theta, nfact)
                    } else stop('Greater than ', opts$MAXQUAD, ' quadrature points.', call.=FALSE)
                    if(opts$message && nfact > 3L && !(opts$odentype %in% c('custom', 'discrete')))
                        message('EM quadrature for high dimensional models are better handled
                                 \twith the \"QMCEM\" or \"MCEM\" method')
                }
            }
            if(!is.null(opts$technical$customTheta)){
                Theta <- opts$technical$customTheta
                if(!is.matrix(Theta)) stop('customTheta input must be a matrix', call.=FALSE)
                opts$quadpts <- nrow(Theta)
            }
            pars <- loadSplinePars(pars, Theta)
        } #end Theta def
        ESTIMATE <- EM.group(pars=pars, constrain=constrain, Ls=Ls, PrepList=PrepList, Data=Data,
                             list = list(NCYCLES=opts$NCYCLES, TOL=opts$TOL, MSTEPTOL=opts$MSTEPTOL,
                                         nfactNames=PrepList[[1L]]$nfactNames, theta=theta,
                                         itemloc=PrepList[[1L]]$itemloc, dentype=opts$dentype,
                                         sitems=sitems, specific=specific, NULL.MODEL=opts$NULL.MODEL,
                                         nfact=nfact, constrain=constrain, verbose=opts$verbose,
                                         SE = opts$SE, SE.type=opts$SE.type, delta=opts$delta, quadpts=opts$quadpts,
                                         accelerate=opts$accelerate, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                                         customPriorFun=opts$customPriorFun, Moptim=opts$Moptim, warn=opts$warn,
                                         message=opts$message, method=opts$method, full=opts$full,
                                         lrPars=lrPars, SE=opts$SE && opts$SE.type == 'numerical', Etable=opts$Etable,
                                         NULL.MODEL=opts$NULL.MODEL, PLCI=opts$PLCI, Norder=opts$Norder,
                                         keep_vcov_PD=opts$keep_vcov_PD, symmetric=opts$technical$symmetric,
                                         MCEM_draws=opts$MCEM_draws, omp_threads=opts$omp_threads),
                             Theta=Theta, DERIV=DERIV, solnp_args=opts$solnp_args, control=control,
                             nconstrain=opts$technical$nconstrain)
        if(opts$method == 'MCEM')
            opts$quadpts <- opts$MCEM_draws(ESTIMATE$cycles)
        opts$Moptim <- ESTIMATE$Moptim
        lrPars <- ESTIMATE$lrPars
        startlongpars <- ESTIMATE$longpars
        rlist <- ESTIMATE$rlist
        logPrior <- ESTIMATE$logPrior
        for(g in seq_len(Data$ngroups)){
            if(g > 1L && opts$dentype == 'mixture') break
            Pl[[g]] <- rlist[[g]]$expected
            Pltmp <- Pl[[g]]
            if(opts$full){
                rg <- 1
                G2group[g] <- NaN
            } else {
                rg <- Data$Freq[[g]]
                Pltmp <- Pltmp[rg != 0]
                rg <- rg[rg != 0]
                Ng <- sum(rg)
                G2group[g] <- 2 * sum(rg * log(rg/(Ng*Pltmp)))
            }
            G2 <- G2 + G2group[g]
            logLik <- logLik + sum(rg*log(Pltmp))
        }
    } else if(opts$method %in% c('MHRM', 'SEM')){ #MHRM estimation
        Theta <- matrix(0, Data$N, nitems)
        if(opts$method == 'SEM') opts$NCYCLES <- NA
        ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, Ls=Ls, PrepList=PrepList, Data=Data,
                               list = list(NCYCLES=opts$NCYCLES, BURNIN=opts$BURNIN,
                                           SEMCYCLES=opts$SEMCYCLES, gain=opts$gain,
                                           KDRAWS=opts$KDRAWS, MHDRAWS=opts$MHDRAWS,
                                           TOL=opts$TOL, SE=FALSE, SE.type = 'none',
                                           nfactNames=PrepList[[1L]]$nfactNames,
                                           itemloc=PrepList[[1L]]$itemloc,
                                           nfact=nfact, constrain=constrain, verbose=opts$verbose,
                                           CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                                           startlongpars=startlongpars,
                                           cand.t.var=opts$technical$MHcand, warn=opts$warn,
                                           message=opts$message, expl=PrepList[[1L]]$exploratory,
                                           plausible.draws=opts$plausible.draws,
                                           MSTEPTOL=opts$MSTEPTOL, Moptim=opts$Moptim,
                                           keep_vcov_PD=opts$keep_vcov_PD),
                               DERIV=DERIV, solnp_args=opts$solnp_args, control=control)
        if(opts$plausible.draws != 0) return(ESTIMATE)
        if(opts$SE && (ESTIMATE$converge || !opts$info_if_converged)){
            if(opts$verbose)
                cat('\nCalculating information matrix...\n')
            tmp <- MHRM.group(pars=ESTIMATE$pars, constrain=constrain, Ls=Ls, PrepList=PrepList, Data=Data,
                                   list = list(NCYCLES=opts$MHRM_SE_draws, BURNIN=1L,
                                               SEMCYCLES=opts$SEMCYCLES, gain=opts$gain,
                                               KDRAWS=opts$KDRAWS, MHDRAWS=opts$MHDRAWS,
                                               TOL=opts$SEtol, SE=TRUE, SE.type=opts$SE.type,
                                               nfactNames=PrepList[[1L]]$nfactNames,
                                               itemloc=PrepList[[1L]]$itemloc,
                                               nfact=nfact, constrain=constrain, verbose=FALSE,
                                               CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                                               startlongpars=ESTIMATE$longpars, plausible.draws=0L,
                                               cand.t.var=opts$technical$MHcand, warn=opts$warn,
                                               message=opts$message, expl=PrepList[[1L]]$exploratory,
                                               MSTEPTOL=opts$MSTEPTOL, Moptim='NR1',
                                               keep_vcov_PD=opts$keep_vcov_PD),
                                   DERIV=DERIV, solnp_args=opts$solnp_args, control=control)
            ESTIMATE$pars <- tmp$pars
            ESTIMATE$info <- tmp$info
            ESTIMATE$fail_invert_info <- tmp$fail_invert_info
            ESTIMATE$time <- c(ESTIMATE$time, SE=sum(tmp$time))
            rm(tmp)
        }
        rlist <- vector('list', Data$ngroups)
        for(g in seq_len(Data$ngroups))
            rlist[[g]]$expected = numeric(1L)
    } else if(opts$method == 'MIXED'){
        if(is.null(opts$technical$RANDSTART)) opts$technical$RANDSTART <- 100L
        if(is.null(opts$technical$BURNIN) && length(mixed.design$random)) opts$BURNIN <- 200L
        Theta <- matrix(0, Data$N, nitems)
        ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, Ls=Ls,
                                    PrepList=PrepList, random=mixed.design$random, Data=Data,
                                    lrPars=lrPars, lr.random=latent.regression$lr.random,
                               list = list(NCYCLES=opts$NCYCLES, BURNIN=opts$BURNIN,
                                           SEMCYCLES=opts$SEMCYCLES, gain=opts$gain,
                                           KDRAWS=opts$KDRAWS, MHDRAWS=opts$MHDRAWS,
                                           TOL=opts$TOL, SE.type = 'none',
                                           nfactNames=PrepList[[1L]]$nfactNames,
                                           itemloc=PrepList[[1L]]$itemloc,
                                           nfact=nfact, constrain=constrain, verbose=opts$verbose,
                                           CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                                           startlongpars=startlongpars, SE=FALSE,
                                           cand.t.var=opts$technical$MHcand, warn=opts$warn,
                                           message=opts$message, expl=FALSE, plausible.draws=0L,
                                           RANDSTART=opts$technical$RANDSTART,
                                           MSTEPTOL=opts$MSTEPTOL, Moptim=opts$Moptim,
                                           keep_vcov_PD=opts$keep_vcov_PD),
                               DERIV=DERIV, solnp_args=opts$solnp_args, control=control)
        if(opts$SE && (ESTIMATE$converge || !opts$info_if_converged)){
            if(opts$verbose)
                cat('\nCalculating information matrix...\n')
            tmp <- MHRM.group(pars=ESTIMATE$pars, constrain=constrain, Ls=Ls,
                              PrepList=PrepList, random=mixed.design$random, Data=Data,
                              lrPars=ESTIMATE$lrPars, lr.random=latent.regression$lr.random,
                              list = list(NCYCLES=opts$MHRM_SE_draws, BURNIN=1L,
                                          SEMCYCLES=opts$SEMCYCLES, gain=opts$gain,
                                          KDRAWS=opts$KDRAWS, MHDRAWS=opts$MHDRAWS,
                                          TOL=opts$SEtol, SE=TRUE, SE.type=opts$SE.type,
                                          nfactNames=PrepList[[1L]]$nfactNames,
                                          itemloc=PrepList[[1L]]$itemloc,
                                          nfact=nfact, constrain=constrain, verbose=FALSE,
                                          CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                                          startlongpars=ESTIMATE$longpars, plausible.draws=0L,
                                          cand.t.var=opts$technical$MHcand, warn=opts$warn,
                                          message=opts$message, expl=FALSE,
                                          RANDSTART=1L,
                                          MSTEPTOL=opts$MSTEPTOL, Moptim='NR1',
                                          keep_vcov_PD=opts$keep_vcov_PD),
                              DERIV=DERIV, solnp_args=opts$solnp_args, control=control)
            ESTIMATE$pars <- tmp$pars
            ESTIMATE$random <- tmp$random
            ESTIMATE$lrPars <- tmp$lrPars
            ESTIMATE$lr.random <- tmp$lr.random
            ESTIMATE$info <- tmp$info
            ESTIMATE$fail_invert_info <- tmp$fail_invert_info
            ESTIMATE$time <- c(ESTIMATE$time, SE=sum(tmp$time))
            rm(tmp)
        }
        rlist <- vector('list', Data$ngroups)
        for(g in seq_len(Data$ngroups))
            rlist[[g]]$expected = numeric(1L)
    }
    for(g in seq_len(length(pars))){
        for(i in seq_len(length(pars[[1L]]))){
            if(class(pars[[g]][[i]]) == 'dich'){
                tmp <- ESTIMATE$pars[[g]][[i]]@par
                nms <- ESTIMATE$pars[[g]][[i]]@parnames
                if(tmp[nms == 'g'] > tmp[nms == 'u']){
                    if(opts$warn)
                        warning('g parameter greater than u detected. Model did not converge', call.=FALSE)
                    ESTIMATE$converge <- FALSE
                }
            }
        }
    }
    opts$times$end.time.Estimate <- proc.time()[3L]
    if(opts$logLik_if_converged && !ESTIMATE$converge) opts$draws <- 0
    opts$times$start.time.SE <- proc.time()[3L]
    if(!opts$NULL.MODEL && opts$SE){
        tmp <- ESTIMATE
        if(opts$verbose && !(opts$method %in% c('MHRM', 'MIXED', 'SEM')))
            cat('\n\nCalculating information matrix...\n')
        if(opts$SE.type %in% c('complete', 'Oakes') && opts$method %in% c('EM', 'QMCEM')){
            opts$times$start.time.SE <- ESTIMATE$start.time.SE
            ESTIMATE <- loadESTIMATEinfo(info=-ESTIMATE$hess, ESTIMATE=ESTIMATE, constrain=constrain,
                                         warn=opts$warn)
        } else if(opts$SE.type == 'SEM' && opts$method == 'EM'){
            collectLL <- as.numeric(ESTIMATE$collectLL)
            collectLL <- exp(c(NA, collectLL) - c(collectLL, NA))
            from <- suppressWarnings(max(which(collectLL <= opts$SEM_from)))
            if(from < 1L) from <- 1L
            to <- min(which(collectLL >= opts$SEM_to))
            dontrun <- FALSE
            if(from == to){
                if(opts$warn)
                    warning('SEM window is too small to compute information matrix.
                            Consider changing the starting values', call.=FALSE)
                dontrun <- TRUE
            }
            lengthsplit <- do.call(c, lapply(strsplit(names(ESTIMATE$correct), 'COV_'), length))
            lengthsplit <- lengthsplit + do.call(c, lapply(strsplit(names(ESTIMATE$correct), 'MEAN_'), length))
            is.latent <- lengthsplit > 2L
            if(!dontrun){
                if(ESTIMATE$cycles <= 10L)
                    if(opts$message)
                        message('Very few EM cycles performed. Consider decreasing TOL further to
                            increase EM iteration count or starting farther away from ML estimates by
                            passing the \'GenRandomPars = TRUE\' argument')
                estmat <- matrix(FALSE, length(ESTIMATE$correction), length(ESTIMATE$correction))
                DM <- estmat + 0
                diag(estmat) <- TRUE
                if(!opts$technical$parallel){
                    ncores <- .mirtClusterEnv$ncores
                    .mirtClusterEnv$ncores <- 1L
                }
                DM <- myLapply(1L:ncol(estmat), FUN=SE.SEM, estmat=estmat, pars=ESTIMATE$pars, constrain=constrain, Data=Data,
                              list = list(NCYCLES=opts$NCYCLES, TOL=opts$SEtol, MSTEPTOL=opts$MSTEPTOL,
                                          nfactNames=PrepList[[1L]]$nfactNames, theta=theta,
                                          itemloc=PrepList[[1L]]$itemloc, keep_vcov_PD=opts$keep_vcov_PD,
                                          sitems=sitems, specific=specific, NULL.MODEL=opts$NULL.MODEL,
                                          nfact=nfact, constrain=constrain, verbose=opts$verbose,
                                          CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND, Moptim=ESTIMATE$Moptim,
                                          EHPrior=ESTIMATE$Prior, warn=opts$warn, dentype=opts$dentype,
                                          message=opts$message, full=opts$full, lrPars=lrPars, method=opts$method),
                              Theta=Theta, theta=theta, ESTIMATE=ESTIMATE, from=from, to=to,
                              DERIV=DERIV, is.latent=is.latent, Ls=Ls, PrepList=PrepList,
                              solnp_args=opts$solnp_args, control=control)
                SEMconv <- sapply(DM, function(x) all(attr(x, 'converged')))
                if(!all(SEMconv)){
                    warning(sprintf(c('%i parameters did not converge in numerical SEM derivative.\n',
                                    'Try using different starting values or passing GenRandomPars=TRUE'),
                                    sum(!SEMconv)),
                            call.=FALSE)
                    SEMconv <- FALSE
                } else SEMconv <- TRUE
                DM <- do.call(rbind, DM)
                if(!opts$technical$parallel)
                    .mirtClusterEnv$ncores <- ncores
                ESTIMATE$pars <- reloadPars(longpars=ESTIMATE$longpars, pars=ESTIMATE$pars,
                                            ngroups=Data$ngroups, J=Data$nitems)
                DM[, is.latent] <- DM[is.latent, ]
                DM[is.latent, is.latent] <- 0
                info <- try(solve(-solve(ESTIMATE$hess) %*% solve(diag(ncol(DM)) - DM)), silent=TRUE)
                info[,is.latent] <- t(info[is.latent, ,drop=FALSE])
                if(opts$technical$symmetric) info <- (info + t(info)) / 2
                if(is(info, 'try-error')){
                    if(opts$warn)
                        warning('Information matrix in SEM could not be computed due to instability',
                                call.=FALSE)
                } else ESTIMATE <- loadESTIMATEinfo(info=info, ESTIMATE=ESTIMATE, constrain=constrain,
                                                    warn=opts$warn)
            }
        } else if(opts$SE.type == 'numerical' && opts$method == 'BL'){
            ESTIMATE <- loadESTIMATEinfo(info=-ESTIMATE$hess, ESTIMATE=ESTIMATE, constrain=constrain,
                                         warn=opts$warn)
        } else if(opts$SE.type %in% c('Richardson', 'forward', 'central') &&
                  !(opts$method %in% c('MHRM', 'SEM', 'MIXED'))){
            ESTIMATE <- SE.Numerical(pars=ESTIMATE$pars, Theta=ESTIMATE$Theta, theta=theta, PrepList=PrepList, Data=Data,
                              dentype=opts$dentype, itemloc=PrepList[[1L]]$itemloc, ESTIMATE=ESTIMATE,
                              constrain=constrain, Ls=Ls, specific=oldmodel, sitems=sitems,
                              CUSTOM.IND=CUSTOM.IND, EHPrior=ESTIMATE$Prior, warn=opts$warn, type=opts$SE.type,
                              delta=opts$delta, lrPars=ESTIMATE$lrPars, omp_threads=opts$omp_threads)
        } else if(opts$SE.type %in% c('MHRM', 'FMHRM') && opts$method == 'EM'){
            if(opts$dentype %in% c('EH', 'EHW'))
                stop('MHRM standard error methods not available when using empirical histograms', call.=FALSE)
            ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, Ls=Ls, PrepList=PrepList, Data=Data,
                                   list = list(NCYCLES=1000L, BURNIN=1L, SEMCYCLES=opts$SEMCYCLES,
                                               KDRAWS=opts$KDRAWS, MHDRAWS=opts$MHDRAWS,
                                               TOL=opts$SEtol, SE=TRUE, SE.type=opts$SE.type,
                                               gain=opts$gain, nfactNames=PrepList[[1L]]$nfactNames,
                                               itemloc=PrepList[[1L]]$itemloc,
                                               nfact=nfact, constrain=constrain, verbose=FALSE, expl=FALSE,
                                               CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND, message=opts$message,
                                               startlongpars=startlongpars, SE=opts$SE, warn=opts$warn,
                                               plausible.draws=0L, MSTEPTOL=opts$MSTEPTOL, Moptim='NR1',
                                               keep_vcov_PD=opts$keep_vcov_PD),
                                   DERIV=DERIV, solnp_args=opts$solnp_args, control=control)
        } else if(any(opts$SE.type %in% c('crossprod', 'Louis', 'sandwich.Louis', 'sandwich')) &&
                  !(opts$method %in% c('MHRM', 'SEM', 'MIXED'))){
            ESTIMATE <- SE.simple(PrepList=PrepList, ESTIMATE=ESTIMATE, Theta=Theta, Data=Data,
                                  constrain=constrain, Ls=Ls, N=nrow(data), type=opts$SE.type,
                                  CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND, warn=opts$warn,
                                  message=opts$message, complete=ESTIMATE$hess)
        } else if(opts$SE.type == 'Fisher' && !(opts$method %in% c('MHRM', 'SEM', 'MIXED'))){
            if(logPrior != 0 && opts$warn)
                warning('Information matrix with the Fisher method does not
                        account for prior parameter distribution information')
            ESTIMATE <- SE.Fisher(PrepList=PrepList, ESTIMATE=ESTIMATE, Theta=Theta, Data=Data,
                                  constrain=constrain, Ls=Ls, full=opts$full,
                                  CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND, warn=opts$warn,
                                  omp_threads=opts$omp_threads)
        }
        ESTIMATE$cycles <- tmp$cycles
        ESTIMATE$Prior <- tmp$Prior
        rm(tmp)
    }
    opts$times$end.time.SE <- proc.time()[3L]
    opts$times$start.time.post <- proc.time()[3L]
    cmods <- vector('list', Data$ngroups)
    if(is.null(opts$theta_lim))
        opts$theta_lim <- numeric(1)
    lrPars <- ESTIMATE$lrPars
    class(lrPars) <- 'S4'
    for(g in seq_len(Data$ngroups)){
        if(opts$method == 'MIXED' || opts$dentype == "discrete"){
            F <- matrix(NA)
            h2 <- numeric(1)
        } else {
            F <- Lambdas(ESTIMATE$pars[[g]], Names=colnames(data))
            colnames(F) <- PrepList[[1L]]$factorNames
            h2 <- rowSums(F^2)
        }
        cmods[[g]] <- new('SingleGroupClass', ParObjects=list(pars=ESTIMATE$pars[[g]], lrPars=lrPars,
                                                              random=ESTIMATE$random,
                                                              lr.random=ESTIMATE$lr.random),
                          Data = list(K=Data$K, nitems=Data$nitems),
                          Model = list(itemloc=PrepList[[1L]]$itemloc, nfact=nfact, constrain=constrain,
                                       factorNames=PrepList[[1L]]$factorNames,
                                       itemtype=PrepList[[1L]]$itemtype,
                                       prodlist=PrepList[[1L]]$prodlist),
                          Options = list(method = 'MHRM', exploratory=PrepList[[1L]]$exploratory,
                                         theta_lim=opts$theta_lim, dentype=opts$dentype),
                          Fit = list(G2=G2group[g], F=F, h2=h2),
                          Internals = list(Pl = rlist[[g]]$expected, CUSTOM.IND=CUSTOM.IND,
                                           SLOW.IND=SLOW.IND))
        if(opts$dentype %in% c("discrete", 'EH', 'EHW', 'Davidian')){
            cmods[[g]]@Model$Theta <- Theta
            cmods[[g]]@Internals$Prior <- list(ESTIMATE$Prior[[g]])
        }
    }
    #missing stats for MHRM
    if(opts$method %in% c('MHRM', 'MIXED', 'SEM') &&
       (!opts$logLik_if_converged || !(!ESTIMATE$converge && opts$logLik_if_converged))){
        logLik <- G2 <- SElogLik <- 0
        if(opts$draws > 0L){
            if(opts$verbose) cat("\nCalculating log-likelihood...\n")
            flush.console()
            if(!opts$technical$parallel){
                ncores <- .mirtClusterEnv$ncores
                .mirtClusterEnv$ncores <- 1L
            }
            for(g in seq_len(Data$ngroups)){
                cmods[[g]]@Data <- list(data=Data$data[Data$group == Data$groupName[g], ],
                                        fulldata=Data$fulldata[[g]], tabdata=Data$tabdata,
                                        Freq=list(Data$Freq[[g]]), K=Data$K)
                cmods[[g]] <- calcLogLik(cmods[[g]], opts$draws, G2 = 'return',
                                         lrPars=ESTIMATE$lrPars)
                cmods[[g]]@Data <- list(K=Data$K, nitems=Data$nitems)
                logLik <- logLik + cmods[[g]]@Fit$logLik
                logPrior <- logPrior + cmods[[g]]@Fit$logPrior
                SElogLik <- SElogLik + cmods[[g]]@Fit$SElogLik
                G2 <- G2 + cmods[[g]]@Fit$G2
                Pl[[g]] <- cmods[[g]]@Internals$Pl
            }
            if(!opts$technical$parallel)
                .mirtClusterEnv$ncores <- ncores
        }
    }

    ####post estimation stats
    if(opts$Moptim %in% c('solnp', 'nloptr')){
        if(!is.null(opts$solnp_args$eqfun))
            df <- df + length(opts$solnp_args$eqfun(ESTIMATE$shortpars, list()))
        if(!is.null(opts$solnp_args$eval_g_eq))
            df <- df + length(opts$solnp_args$eval_g_eq(ESTIMATE$shortpars, list()))
    }
    r <- rr
    N <- sum(r)
    tmp <- dfsubtr
    AIC <- (-2) * logLik + 2 * tmp
    BIC <- (-2) * logLik + tmp*log(N)
    SABIC <- (-2) * logLik + tmp*log((N+2)/24)
    HQ <- (-2) * logLik + 2*tmp*log(log(N))
    p.G2 <- 1 - pchisq(G2,df)
    RMSEA.G2 <- rmsea(X2=G2, df=df, N=N)
    null.mod <- unclass(new('SingleGroupClass'))
    TLI.G2 <- CFI.G2 <- NaN
    if(opts$calcNull && length(r) * 3L < prod(Data$K) && opts$warn)
        warning(c('Full table of responses is very sparse. ',
                'Goodness-of-fit statistics may be very inaccurate'), call.=FALSE)
    if(!opts$NULL.MODEL && opts$method != 'MIXED' && opts$calcNull && nmissingtabdata == 0L){
        null.mod <- try(unclass(computeNullModel(data=data, itemtype=itemtype, key=key,
                                                 group=if(length(pars) > 1L) group else NULL)))
        if(is(null.mod, 'try-error')){
            if(opts$warn)
                warning('Null model calculation did not converge.')
            null.mod <- unclass(new('SingleGroupClass'))
        } else if(!is.nan(G2)) {
            TLI.G2 <- tli(X2=G2, X2.null=null.mod@Fit$G2, df=df, df.null=null.mod@Fit$df)
            CFI.G2 <- cfi(X2=G2, X2.null=null.mod@Fit$G2, df=df, df.null=null.mod@Fit$df)
        }
    }
    if(nmissingtabdata > 0L)
        p.G2 <- RMSEA.G2 <- G2 <- TLI.G2 <- CFI.G2 <-  NaN
    if(is.null(parprior)) parprior <- list()
    if(is.null(opts$quadpts)) opts$quadpts <- NaN
    opts$times$end.time.post <- proc.time()[3L]
    time <- opts$times
    opts$times <- NULL
    # set macro objects
    Options <- opts
    Options$exploratory <- PrepList[[1L]]$exploratory
    Fit <- list(G2=G2, p=p.G2, TLI=TLI.G2, CFI=CFI.G2, RMSEA=RMSEA.G2, df=df,
                AIC=AIC, BIC=BIC, SABIC=SABIC, HQ=HQ, logLik=logLik,
                logPrior=logPrior, SElogLik=SElogLik, F=F, h2=h2)
    pis <- if(opts$dentype == 'mixture')
        ExtractMixtures(lapply(cmods, function(x) x@ParObjects$pars)) else NULL
    Model <- list(model=oldmodel, factorNames=PrepList[[1L]]$factorNames, itemtype=PrepList[[1L]]$itemtype,
                  itemloc=PrepList[[1L]]$itemloc, nfact=nfact, pis=pis,
                  Theta=Theta, constrain=constrain, parprior=parprior, nest=as.integer(dfsubtr),
                  invariance=invariance, lrPars=lrPars, formulas=attr(mixed.design, 'formula'),
                  prodlist=PrepList[[1L]]$prodlist, nestpars=nestpars)
    if(!is.null(opts$technical$Etable)){
        Model$Etable <- ESTIMATE$rlist
        Model$Etable$Theta <- Theta
    }
    Data$covdata <- if(length(lrPars)) lrPars@df else attr(mixed.design, 'covdata')
    Data$itemdesign <- attr(mixed.design, 'itemdesign')
    ParObjects <- list(pars=cmods, lrPars=lrPars, random=ESTIMATE$random, lr.random=ESTIMATE$lr.random)
    OptimInfo <- list(iter=ESTIMATE$cycles, converged=ESTIMATE$converge, cand.t.var=ESTIMATE$cand.t.var,
                      condnum=NA, secondordertest=NA, SEMconv=SEMconv, aveAR=ESTIMATE$aveAR)
    vcov <- matrix(NA, 1, 1)
    if(Options$SE){
        information <- ESTIMATE$info
        if(!is.null(opts$technical$infoAsVcov)){
            vcov <- information
        } else {
            if(!ESTIMATE$fail_invert_info){
                isna <- is.na(diag(information))
                info <- information[!isna, !isna]
                vcov <- matrix(NA, ncol(information), ncol(information))
                rownames(vcov) <- colnames(vcov) <- colnames(information)
                vcov2 <- try(solve(info), silent=TRUE)
                vcov[!isna, !isna] <- vcov2
                if(!is(vcov2, 'try-error')){
                    OptimInfo$condnum <- kappa(info, exact=TRUE)
                    OptimInfo$secondordertest <- all(eigen(info)$values > 0)
                    } else OptimInfo$secondordertest <- FALSE
            } else {
                OptimInfo$secondordertest <- FALSE
            }
        }
    } else OptimInfo$secondordertest <- NA
    Internals <- list(collectLL=ESTIMATE$collectLL, Prior=ESTIMATE$Prior, Pl=Pl,
                      shortpars=as.numeric(ESTIMATE$shortpars), key=key,
                      bfactor=list(), CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                      survey.weights=survey.weights)
    if(opts$storeEtable)
        Internals$Etable <- ESTIMATE$Etable
    if(opts$method == 'SEM') Options$TOL <- NA
    if(opts$odentype == "discrete"){
        Fit$F <- Fit$h2 <- NULL
        mod <- new('DiscreteClass',
                   Data=Data,
                   Options=Options,
                   Fit=Fit,
                   Model=Model,
                   ParObjects=ParObjects,
                   OptimInfo=OptimInfo,
                   Internals=Internals,
                   vcov=vcov)
    } else {
        if(Data$ngroups == 1L){
            ParObjects$pars <- cmods[[1L]]@ParObjects$pars
            Internals$Pl <- Internals$Pl[[1L]]
            if(opts$method == 'MIXED'){
                Fit$p <- NaN
                mod <- new('MixedClass',
                           Data=Data,
                           Options=Options,
                           Fit=Fit,
                           Model=Model,
                           ParObjects=ParObjects,
                           OptimInfo=OptimInfo,
                           Internals=Internals,
                           vcov=vcov)
            } else {
                if(Options$exploratory){
                    FF <- F %*% t(F)
                    V <- eigen(FF)$vector[ ,1L:nfact]
                    L <- eigen(FF)$values[1L:nfact]
                    if (nfact == 1L) F <- as.matrix(V * sqrt(L))
                    else F <- V %*% sqrt(diag(L))
                    if (sum(F[ ,1L] < 0)) F <- (-1) * F
                    colnames(F) <- paste("F", 1L:ncol(F), sep="")
                    h2 <- rowSums(F^2)
                } else {
                    if(opts$method == 'EM')
                        Internals$bfactor <- list(prior=ESTIMATE$prior,
                                                  Priorbetween=ESTIMATE$Priorbetween,
                                                  sitems=ESTIMATE$sitems, specific=specific)
                }
                mod <- new('SingleGroupClass',
                           Data=Data,
                           Options=Options,
                           Fit=Fit,
                           Model=Model,
                           ParObjects=ParObjects,
                           OptimInfo=OptimInfo,
                           Internals=Internals,
                           vcov=vcov)
            }
        } else {
            if(opts$method == 'EM')
                Internals$bfactor <- list(prior=ESTIMATE$prior,
                                          Priorbetween=ESTIMATE$Priorbetween,
                                          sitems=ESTIMATE$sitems, specific=specific)
            cls <- ifelse(opts$odentype == 'mixture', 'MixtureClass', 'MultipleGroupClass')
            mod <- new(cls,
                       Data=Data,
                       Options=Options,
                       Fit=Fit,
                       Model=Model,
                       ParObjects=ParObjects,
                       OptimInfo=OptimInfo,
                       Internals=Internals,
                       vcov=vcov)
        }
    }
    mod@time <- c("TOTAL:" = as.numeric(proc.time()[3L] - time$start.time),
                  Data = as.numeric(time$end.time.Data - time$start.time.Data),
                  ESTIMATE$time,
                  SE = as.numeric(time$end.time.SE - time$start.time.SE),
                  Post = as.numeric(time$end.time.post - time$start.time.post))
    return(mod)
}
