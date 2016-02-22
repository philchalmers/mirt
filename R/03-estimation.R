ESTIMATION <- function(data, model, group, itemtype = NULL, guess = 0, upper = 1,
                       invariance = '', pars = NULL, constrain = NULL, key = NULL,
                       parprior = NULL, mixed.design = NULL, customItems = NULL,
                       GenRandomPars = FALSE, large = FALSE,
                       survey.weights = NULL, discrete=FALSE, latent.regression = NULL,
                       gpcm_mats=list(), control = list(), ...)
{
    start.time <- proc.time()[3L]
    dots <- list(...)
    if(missing(data)) missingMsg('data')
    if(missing(model)) missingMsg('model')
    if(missing(group)) missingMsg('group')
    if(is.logical(large) && large){
        Data <- opts <- list()
        data <- as.matrix(data)
        if(is.numeric(data))
            data <- matrix(as.integer(data), nrow(data), dimnames=list(rownames(data), colnames(data)))
        rownames(data) <- 1L:nrow(data)
        if(is.null(colnames(data)))
            colnames(data) <- paste0('Item.', 1L:ncol(data))
        Data$data <- data
        Data$group <- factor(group)
        Data$groupNames <- unique(Data$group)
        Data$ngroups <- length(Data$groupNames)
        Data$nitems <- ncol(data)
        K <- apply(Data$data, 2L, function(x) length(unique(na.omit(x))))
        tmp <- list(...)
        if(!is.null(tmp$technical$customK)) K <- tmp$technical$customK
        itemloc <- c(1L, cumsum(K) + 1L)
        PrepListFull <- list(K=K, itemloc=itemloc, Names=NULL, itemnames=colnames(Data$data))
    } else {
        if(missing(data) || is.null(nrow(data)))
            stop('data argument is required', call.=FALSE)
        if(missing(model) || !is(model, 'numeric') && !is(model, 'mirt.model') &&
           !is(model, 'character'))
            stop('model argument (numeric, character, or from mirt.model() function) is required',
                 call.=FALSE)
        if(!(is.factor(group) || is.character(group)) || length(group) != nrow(data))
            stop('group input provided is not valid', call.=FALSE)
        if(!is.null(itemtype))
            stopifnot(is(itemtype, 'character'))
        if(!is.null(constrain))
            stopifnot(is(constrain, 'list'))
        if(!is.null(parprior))
            stopifnot(is(parprior, 'list'))
        if(!is.null(customItems))
            stopifnot(is(customItems, 'list'))
        stopifnot(is(invariance, 'character'))
        stopifnot(is(GenRandomPars, 'logical'))
        stopifnot(is(large, 'logical') || is(large, 'list'))
        opts <- makeopts(GenRandomPars=GenRandomPars, ...)
        if(!is.null(survey.weights)){
            stopifnot(opts$method == 'EM')
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
        if(opts$method == 'MHRM' || opts$method == 'MIXED')
            on.exit(set.seed((as.numeric(Sys.time()) - floor(as.numeric(Sys.time()))) * 1e8))
        #change itemtypes if NULL.MODEL
        if(opts$NULL.MODEL){
            constrain <- NULL
            if(!is.null(itemtype)){
                itemtype[itemtype == 'grsm'] <- 'graded'
                itemtype[itemtype == 'rsm'] <- 'gpcm'
                itemtype[itemtype == '3PL' | itemtype == '3PLu' | itemtype == '4PL'] <- '2PL'
                itemtype[itemtype == '3PLNRM' | itemtype == '3PLuNRM' | itemtype == '4PLNRM'] <- '2PLNRM'
            }
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
        if(nrow(data) > 1L){
            data <- apply(data, 2L, function(x, message){
                s <- sort(unique(x))
                se <- min(s, na.rm = TRUE):max(x, na.rm = TRUE)
                if(length(s) != length(se)){
                    if(message)
                        message('Item re-scored so that all values are within a distance of 1')
                    for(i in 2L:length(s))
                        x <- ifelse(x == s[i], se[i], x)
                }
                x
            }, message = opts$message)
        }
        if(any(rowSums(is.na(data)) == ncol(data))){
            if(!opts$removeEmptyRows)
                stop('data contains completely empty response patterns.',
                     'Please remove manually or pass removeEmptyRows=TRUE to the technical list',
                     call.=FALSE)
            else {
                pick <- rowSums(is.na(data)) != ncol(data)
                data <- subset(data, pick)
                group <- subset(group, pick)
                if(!is.null(latent.regression) || !is.null(mixed.design))
                    stop('removeEmptyRows input not supported for latent regression/mixed effect models.',
                         'Please remove the require rows manually for each object.', call.=FALSE)
            }
        }
        Data$data <- data

        if(is.null(opts$grsm.block)) Data$grsm.block <- rep(1L, ncol(data))
        else Data$grsm.block <- opts$grsm.block
        # if(is.null(opts$rsm.block)) Data$rsm.block <- rep(1L, ncol(data))
        Data$group <- factor(group)
        Data$groupNames <- unique(Data$group)
        if(any(grepl('-', Data$groupNames)))
            stop('Group names cannot contain a dash (-) character', call.=FALSE)
        Data$ngroups <- length(Data$groupNames)
        Data$nitems <- ncol(Data$data)
        Data$N <- nrow(Data$data)
        Data$mins <- apply(data, 2L, min, na.rm=TRUE)
        if(is.character(model)){
            tmp <- any(sapply(colnames(data), grepl, x=model))
            model <- mirt.model(model, itemnames = if(tmp) colnames(data) else NULL)
        }
        oldmodel <- model
        if(length(model) == 1L){
            newmodel <- list()
            for(g in 1L:Data$ngroups)
                newmodel[[g]] <- model
            names(newmodel) <- Data$groupNames
            model <- newmodel
        }
        Data$model <- model
        PrepList <- vector('list', Data$ngroups)
        names(PrepList) <- Data$groupNames
        tmp <- 1L:Data$ngroups
        selectmod <- Data$model[[tmp[names(Data$model) == Data$groupNames[1L]]]]
        if(!is.null(dots$PrepList)) {
            PrepListFull <- PrepList[[1L]] <- dots$PrepList
        } else {
            PrepListFull <- PrepList[[1L]] <-
                PrepData(data=Data$data, model=selectmod, itemtype=itemtype, guess=guess,
                         upper=upper, parprior=parprior, verbose=opts$verbose,
                         technical=opts$technical, parnumber=1L, BFACTOR=opts$BFACTOR,
                         grsm.block=Data$grsm.block, rsm.block=Data$rsm.block,
                         mixed.design=mixed.design, customItems=customItems,
                         fulldata=opts$PrepList[[1L]]$fulldata, key=key,
                         gpcm_mats=gpcm_mats, internal_constraints=opts$internal_constraints)
            if(!is.null(dots$Return_PrepList)) return(PrepListFull)
        }
        parnumber <- max(PrepList[[1L]]$pars[[Data$nitems+1L]]@parnum) + 1L
        attr(PrepListFull$pars, 'nclasspars') <- attr(PrepList[[1L]]$pars, 'nclasspars') <-
            sapply(PrepListFull$pars, function(y) length(y@parnum))
        for(g in 1L:Data$ngroups){
            if(g != 1L){
                PrepList[[g]] <- list(pars=PrepList[[1L]]$pars)
                for(i in 1L:length(PrepList[[g]]$pars)){
                    PrepList[[g]]$pars[[i]]@parnum <- parnumber:(parnumber + length(PrepList[[g]]$pars[[i]]@parnum) - 1L)
                    parnumber <- max(PrepList[[g]]$pars[[i]]@parnum) + 1L
                }
            }
        }
        if(length(mixed.design$random) > 0L){
            for(i in 1L:length(mixed.design$random)){
                mixed.design$random[[i]]@parnum <- parnumber:(parnumber - 1L +
                                                                  length(mixed.design$random[[i]]@par))
                parnumber <- max(mixed.design$random[[i]]@parnum) + 1L
            }
        }
        if(!is.null(latent.regression)){
            if(length(PrepListFull$prodlist))
                stop('Polynominal combinations currently not supported when latent regression effects are used', call.=FALSE)
            lrPars <- make.lrdesign(df=latent.regression$df, formula=latent.regression$formula,
                                    factorNames=PrepListFull$factorNames, EM=latent.regression$EM,
                                    TOL=opts$TOL)
            lrPars@parnum <- parnumber:(parnumber - 1L + length(lrPars@par))
            parnumber <- max(lrPars@parnum) + 1L
        } else lrPars <- list()
        if(length(latent.regression$lr.random) > 0L){
            stop('lr.random input not yet supported', call.=FALSE)
        }
    }
    if(!is.null(opts$PrepList)){
        PrepList <- opts$PrepList
    } else {
        if(!is.list(large)){
            tmpdata <- Data$data
            tmpdata[is.na(tmpdata)] <- 99999L
            stringfulldata <- apply(tmpdata, 1L, paste, sep='', collapse = '/')
            stringtabdata <- unique(stringfulldata)
            tmptabdata <- maketabData(stringfulldata=stringfulldata, stringtabdata=stringtabdata,
                                      group=Data$group, groupNames=Data$groupNames, nitem=Data$nitems,
                                      K=PrepListFull$K, itemloc=PrepListFull$itemloc,
                                      Names=PrepListFull$Names, itemnames=PrepListFull$itemnames,
                                      survey.weights=survey.weights)
            if(is.logical(large) && large){
                return(tmptabdata)
            } else large <- tmptabdata
        }
        Data$tabdatalong <- large$tabdata
        Data$tabdata <- large$tabdata2
        for(g in 1L:Data$ngroups){
            select <- Data$group == Data$groupNames[g]
            Data$fulldata[[g]] <- PrepListFull$fulldata[select, , drop=FALSE]
            Data$Freq[[g]] <- large$Freq[[g]]
        }
    }
    if(opts$returnPrepList) return(PrepList)
    if(opts$BFACTOR){
        #better start values
        if((PrepList[[1L]]$nfact - attr(model[[1L]], 'nspec')) == 1L){
            nfact <- PrepListFull$nfact
            for(g in 1L:Data$ngroups){
                for(i in 1L:Data$nitems){
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
    PrepList <- UpdatePrior(PrepList, model, groupNames=Data$groupNames, warn=opts$warn)
    if(GenRandomPars){
        for(g in 1L:Data$ngroups)
            for(i in 1L:length(PrepList[[g]]$pars))
                PrepList[[g]]$pars[[i]] <- GenRandomPars(PrepList[[g]]$pars[[i]])
    }
    if(discrete){
        PrepList[[1L]]$exploratory <- FALSE
        if(is.null(opts$technical$customTheta))
            opts$technical$customTheta <- diag(PrepList[[1L]]$nfact)
    }
    RETURNVALUES <- FALSE
    if(!is.null(pars)){
        if(is(pars, 'data.frame')){
            PrepList <- UpdatePrepList(PrepList, pars, random=mixed.design$random,
                                       lrPars=lrPars, MG = TRUE)
            mixed.design$random <- attr(PrepList, 'random')
            if(any(pars$class == 'lrPars')) lrPars <- update.lrPars(pars, lrPars)
            attr(PrepList, 'random') <- NULL
        }
        if(!is.null(attr(pars, 'values')) || (is.character(pars) && pars == 'values'))
            RETURNVALUES <- TRUE
    }
    pars <- vector('list', Data$ngroups)
    for(g in 1L:Data$ngroups)
        pars[[g]] <- PrepList[[g]]$pars
    nitems <- Data$nitems
    Data$K <- PrepList[[1L]]$K
    nfact <- PrepList[[1L]]$pars[[nitems+1L]]@nfact
    if(nfact != 1L && any(c('Rasch') %in% itemtype ) && PrepList[[1L]]$exploratory)
       stop('Rasch itemtypes are for confimatory models only.', call.=FALSE)
    nLambdas <- PrepList[[1L]]$pars[[1L]]@nfact
    if(is.null(constrain)) constrain <- list()
    nspec <- ifelse(!is.null(attr(model[[1L]], 'nspec')), attr(model[[1L]], 'nspec'), 1L)
    #default MG uses configural model (independent groups but each identified)
    if('free_means' %in% invariance ){ #Free factor means (means 0 for ref)
        if(opts$BFACTOR){
            for(g in 2L:Data$ngroups)
                pars[[g]][[nitems + 1L]]@est[1L:(nfact-nspec)] <- TRUE
        } else {
            for(g in 2L:Data$ngroups)
                pars[[g]][[nitems + 1L]]@est[1L:nfact] <- TRUE
        }
    }
    dummymat <- matrix(FALSE, pars[[1L]][[nitems + 1L]]@nfact, pars[[1L]][[nitems + 1L]]@nfact)
    if(any('free_var' %in% invariance)){ #Free factor vars (vars 1 for ref)
        if(opts$BFACTOR){
            tmp <- dummymat[1L:(nfact-nspec),1L:(nfact-nspec), drop=FALSE]
            diag(tmp) <- TRUE
            dummymat[1L:(nfact-nspec),1L:(nfact-nspec)] <- tmp
        } else {
            diag(dummymat) <- TRUE
        }
        tmp <- dummymat[lower.tri(dummymat, TRUE)]
        for(g in 2L:Data$ngroups){
            pars[[g]][[nitems + 1L]]@est <- c(pars[[g]][[nitems + 1L]]@est[1L:pars[[g]][[nitems + 1L]]@nfact], tmp)
            names(pars[[g]][[nitems + 1L]]@est) <- names(pars[[g]][[nitems + 1L]]@par)
        }
    }
    if(RETURNVALUES){
        for(g in 1L:Data$ngroups)
            PrepList[[g]]$pars <- pars[[g]]
        return(ReturnPars(PrepList, PrepList[[1L]]$itemnames,
                          random=mixed.design$random, lrPars=lrPars, MG = TRUE))

    }
    constrain <- UpdateConstrain(pars=pars, constrain=constrain, invariance=invariance, nfact=Data$nfact,
                                 nLambdas=nLambdas, J=nitems, ngroups=Data$ngroups, PrepList=PrepList,
                                 method=opts$method, itemnames=PrepList[[1L]]$itemnames, model=model,
                                 groupNames=Data$groupNames)
    startlongpars <- c()
    if(opts$NULL.MODEL){
        constrain <- list()
        for(g in 1L:length(pars)){
            for(i in 1L:nitems){
                pars[[g]][[i]] <- set_null_model(pars[[g]][[i]])
            }
        }
    }
    rr <- 0L
    for(g in 1L:Data$ngroups){
        r <- Data$Freq[[g]]
        rr <- rr + r
    }
    df <- prod(PrepList[[1L]]$K) - 1
    if(df > 1e10) df <- 1e10
    nestpars <- nconstr <- 0L
    for(g in 1L:Data$ngroups)
        for(i in 1L:(nitems+1L))
            nestpars <- nestpars + sum(pars[[g]][[i]]@est)
    if(!is.null(mixed.design$random)){
        for(i in 1L:length(mixed.design$random))
            nestpars <- nestpars + sum(mixed.design$random[[i]]@est)
    }
    if(length(lrPars))
        nestpars <- nestpars + sum(lrPars@est)
    if(length(constrain) > 0L)
        for(i in 1L:length(constrain))
            nconstr <- nconstr + length(constrain[[i]]) - 1L
    if(discrete) nestpars <- nestpars + nrow(opts$technical$customTheta) - 1L
    nmissingtabdata <- sum(is.na(rowSums(Data$tabdata)))
    dfsubtr <- nestpars - nconstr
    if(df <= dfsubtr)
        stop('Too few degrees of freedom. There are only ', df, ' degrees of freedom but ',
             dfsubtr, ' parameters were freely estimated.', call.=FALSE)
    df <- df - dfsubtr
    if(!is.null(customItems)){
        for(g in 1L:Data$ngroups)
            PrepList[[g]]$exploratory <- FALSE
    }
    G2group <- numeric(Data$ngroups)
    DERIV <- vector('list', Data$ngroups)
    for(g in 1L:Data$ngroups){
        DERIV[[g]] <- vector('list', Data$nitems)
        for(i in 1L:Data$nitems)
            DERIV[[g]][[i]] <- selectMethod(Deriv, c(class(pars[[g]][[i]]), 'matrix'))
    }
    Ls <- makeLmats(pars, constrain, random = mixed.design$random, lrPars=lrPars)
    CUSTOM.IND <- which(sapply(pars[[1L]], class) %in% Use_R_ProbTrace())
    SLOW.IND <- which(sapply(pars[[1L]], class) %in% Use_R_Deriv())
    #warnings
    wmsg <- 'Lower and upper bound parameters (g and u) should use \'norm\' (i.e., logit) prior'
    for(g in 1L:length(pars)){
        for(i in 1L:length(pars[[1L]])){
            if(class(pars[[g]][[i]]) == 'dich'){
                pt <- pars[[g]][[i]]@prior.type
                if(!(pt[length(pt)-1L] %in% c(0L, 1L))) warning(wmsg, call.=FALSE)
                if(!(pt[length(pt)] %in% c(0L, 1L))) warning(wmsg, call.=FALSE)
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
    if(opts$method == 'EM' || opts$method == 'BL' || opts$method == 'QMCEM'){
        if(length(lrPars)){
            if(opts$SE && !(opts$SE.type %in% c('complete', 'forward', 'central', 'Richardson'))) ## TODO
                stop('Information matrix method for latent regression estimates not supported',
                     call.=FALSE)
            opts$full <- TRUE
        } else opts$full <- FALSE
        temp <- matrix(0L,nrow=nitems,ncol=nspec)
        sitems <- matrix(0L, nrow=sum(PrepList[[1L]]$K), ncol=nspec)
        specific <- NULL
        if(discrete){
            theta <- 0
            Theta <- opts$technical$customTheta
            opts$quadpts <- nrow(Theta)
            opts$customPriorFun <- lca_prior
        } else {
            if(is.null(opts$quadpts)){
                tmp <- if(opts$BFACTOR) PrepList[[1L]]$nfact - attr(model[[1L]], 'nspec') + 1L
                    else nfact
                opts$quadpts <- select_quadpts(tmp)
            }
            if(opts$quadpts < 3 && opts$warn) warning('Should use more than 2 quadpts', call.=FALSE)
            theta <- 1
            if(opts$method != 'QMCEM')
                theta <- as.matrix(seq(opts$theta_lim[1L], opts$theta_lim[2L],
                                       length.out = opts$quadpts))
            if(opts$BFACTOR){
                for(g in 1L:length(pars)){
                    tmp <- pars[[g]][[nitems+1L]]
                    gp <- ExtractGroupPars(tmp)
                    tmp@par <- as.numeric(tmp@est)
                    gp2 <- ExtractGroupPars(tmp)
                    ind <- (length(gp$gmeans)-ncol(sitems)+1L):length(gp$gmeans)
                    gmeans <- gp$gmeans[ind]
                    gcov <- gp$gcov[ind, ind, drop=FALSE]
                    lgmeans <- gp2$gmeans[ind]
                    lgcov <- gp2$gcov[ind, ind, drop=FALSE]
                    if(any(lgmeans != 0) || any(lgcov != 0))
                        stop('Cannot freely estimate specific factor parameters in
                             quadrature reduction method', call.=FALSE)
                    if(!all(gmeans == 0))
                        stop('Means for specific factors must be 0', call.=FALSE)
                    if(!all(gcov == diag(ncol(sitems))))
                        stop('Covariance matrix for specific factors must be an identity matrix',
                             call.=FALSE)
                }
                specific <- attr(oldmodel, 'specific')
                specific[is.na(specific)] <- 1L
                for(i in 1L:nitems) temp[i, specific[i]] <- 1L
                ind <- 1L
                for(i in 1L:nitems){
                    for(j in 1L:PrepList[[1L]]$K[i]){
                        sitems[ind, ] <- temp[i, ]
                        ind <- ind + 1L
                    }
                }
                nfact2 <- PrepList[[1L]]$nfact - attr(model[[1L]], 'nspec') + 1L
                Theta <- thetaComb(theta, nfact2)
                Theta <- cbind(Theta[,1L:(nfact2-1L),drop=FALSE],
                               matrix(Theta[,nfact2], nrow=nrow(Theta), ncol=ncol(sitems)))
            } else {
                if(opts$method == 'QMCEM'){
                    Theta <- QMC_quad(npts=opts$quadpts, nfact=nfact, lim=opts$theta_lim,
                                      norm=TRUE)
                } else {
                    if(opts$quadpts^nfact <= opts$MAXQUAD){
                        if(is.null(opts$technical$customTheta))
                            Theta <- thetaComb(theta, nfact)
                    } else stop('Greater than ', opts$MAXQUAD, ' quadrature points.', call.=FALSE)
                    if(opts$message && nfact > 3L)
                        message('EM quadrature for high dimensional models are better handled
                                 \twith the \"QMCEM\" method')
                }
            }
            if(!is.null(opts$technical$customTheta)){
                Theta <- opts$technical$customTheta
                if(!is.matrix(Theta)) stop('customTheta input must be a matrix', call.=FALSE)
                opts$quadpts <- nrow(Theta)
            }
        } #end Theta def
        ESTIMATE <- EM.group(pars=pars, constrain=constrain, Ls=Ls, PrepList=PrepList, Data=Data,
                             list = list(NCYCLES=opts$NCYCLES, TOL=opts$TOL, MSTEPTOL=opts$MSTEPTOL,
                                         nfactNames=PrepList[[1L]]$nfactNames, theta=theta, EH=opts$empiricalhist,
                                         itemloc=PrepList[[1L]]$itemloc, BFACTOR=opts$BFACTOR,
                                         sitems=sitems, specific=specific, NULL.MODEL=opts$NULL.MODEL,
                                         nfact=nfact, constrain=constrain, verbose=opts$verbose,
                                         SEM=any(opts$SE.type %in% c('SEM', 'complete')) && opts$SE,
                                         accelerate=opts$accelerate, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                                         customPriorFun=opts$customPriorFun, Moptim=opts$Moptim, warn=opts$warn,
                                         message=opts$message, BL=opts$method == 'BL', full=opts$full,
                                         lrPars=lrPars, SE=opts$SE && opts$SE.type == 'numerical', Etable=opts$Etable),
                             Theta=Theta, DERIV=DERIV, solnp_args=opts$solnp_args, control=control)
        opts$Moptim <- ESTIMATE$Moptim
        lrPars <- ESTIMATE$lrPars
        startlongpars <- ESTIMATE$longpars
        rlist <- ESTIMATE$rlist
        logLik <- G2 <- SElogLik <- 0
        logPrior <- ESTIMATE$logPrior
        for(g in 1L:Data$ngroups){
            Pl <- rlist[[g]]$expected
            if(opts$full){
                rg <- 1
                G2group[g] <- NaN
            } else {
                rg <- Data$Freq[[g]]
                Pl <- Pl[rg != 0]
                rg <- rg[rg != 0]
                Ng <- sum(rg)
                G2group[g] <- 2 * sum(rg * log(rg/(Ng*Pl)))
            }
            G2 <- G2 + G2group[g]
            logLik <- logLik + sum(rg*log(Pl))
        }
        Pl <- list(Pl)
    } else if(opts$method == 'MHRM'){ #MHRM estimation
        Theta <- matrix(0, Data$N, nitems)
        ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, Ls=Ls, PrepList=PrepList, Data=Data,
                               list = list(NCYCLES=opts$NCYCLES, BURNIN=opts$BURNIN,
                                           SEMCYCLES=opts$SEMCYCLES, gain=opts$gain,
                                           KDRAWS=opts$KDRAWS, MHDRAWS=opts$MHDRAWS,
                                           TOL=opts$TOL, SE=FALSE, SE.type = 'none',
                                           nfactNames=PrepList[[1L]]$nfactNames,
                                           itemloc=PrepList[[1L]]$itemloc, BFACTOR=opts$BFACTOR,
                                           nfact=nfact, constrain=constrain, verbose=opts$verbose,
                                           CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                                           startlongpars=startlongpars,
                                           cand.t.var=opts$technical$MHcand, warn=opts$warn,
                                           message=opts$message, expl=PrepList[[1L]]$exploratory),
                               DERIV=DERIV)
        if(opts$SE){
            if(opts$verbose)
                cat('\nCalculating information matrix...\n')
            tmp <- MHRM.group(pars=ESTIMATE$pars, constrain=constrain, Ls=Ls, PrepList=PrepList, Data=Data,
                                   list = list(NCYCLES=opts$MHRM_SE_draws, BURNIN=1L,
                                               SEMCYCLES=opts$SEMCYCLES, gain=opts$gain,
                                               KDRAWS=opts$KDRAWS, MHDRAWS=opts$MHDRAWS,
                                               TOL=opts$SEtol, SE=TRUE, SE.type=opts$SE.type,
                                               nfactNames=PrepList[[1L]]$nfactNames,
                                               itemloc=PrepList[[1L]]$itemloc, BFACTOR=opts$BFACTOR,
                                               nfact=nfact, constrain=constrain, verbose=FALSE,
                                               CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                                               startlongpars=ESTIMATE$longpars,
                                               cand.t.var=opts$technical$MHcand, warn=opts$warn,
                                               message=opts$message, expl=PrepList[[1L]]$exploratory),
                                   DERIV=DERIV)
            ESTIMATE$pars <- tmp$pars
            ESTIMATE$info <- tmp$info
            ESTIMATE$fail_invert_info <- tmp$fail_invert_info
            ESTIMATE$time <- c(ESTIMATE$time, SE=sum(tmp$time))
        }
        rlist <- vector('list', Data$ngroups)
        for(g in 1L:Data$ngroups)
            rlist[[g]]$expected = numeric(1L)
    } else if(opts$method == 'MIXED'){
        if(is.null(opts$technical$RANDSTART)) opts$technical$RANDSTART <- 100L
        if(is.null(opts$technical$BURNIN) && length(mixed.design$random)) opts$BURNIN <- 200L
        Theta <- matrix(0, Data$N, nitems)
        ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, Ls=Ls,
                                    PrepList=PrepList, random=mixed.design$random, Data=Data,
                                    lrPars=lrPars,
                               list = list(NCYCLES=opts$NCYCLES, BURNIN=opts$BURNIN,
                                           SEMCYCLES=opts$SEMCYCLES, gain=opts$gain,
                                           KDRAWS=opts$KDRAWS, MHDRAWS=opts$MHDRAWS,
                                           TOL=opts$TOL, SE.type = 'none',
                                           nfactNames=PrepList[[1L]]$nfactNames,
                                           itemloc=PrepList[[1L]]$itemloc, BFACTOR=opts$BFACTOR,
                                           nfact=nfact, constrain=constrain, verbose=opts$verbose,
                                           CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                                           startlongpars=startlongpars, SE=FALSE,
                                           cand.t.var=opts$technical$MHcand, warn=opts$warn,
                                           message=opts$message, expl=FALSE,
                                           RANDSTART=opts$technical$RANDSTART),
                               DERIV=DERIV)
        if(opts$SE){
            if(opts$verbose)
                cat('\nCalculating information matrix...\n')
            tmp <- MHRM.group(pars=ESTIMATE$pars, constrain=constrain, Ls=Ls,
                              PrepList=PrepList, random=mixed.design$random, Data=Data,
                              lrPars=ESTIMATE$lrPars,
                              list = list(NCYCLES=opts$MHRM_SE_draws, BURNIN=1L,
                                          SEMCYCLES=opts$SEMCYCLES, gain=opts$gain,
                                          KDRAWS=opts$KDRAWS, MHDRAWS=opts$MHDRAWS,
                                          TOL=opts$SEtol, SE=TRUE, SE.type=opts$SE.type,
                                          nfactNames=PrepList[[1L]]$nfactNames,
                                          itemloc=PrepList[[1L]]$itemloc, BFACTOR=opts$BFACTOR,
                                          nfact=nfact, constrain=constrain, verbose=FALSE,
                                          CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                                          startlongpars=ESTIMATE$longpars,
                                          cand.t.var=opts$technical$MHcand, warn=opts$warn,
                                          message=opts$message, expl=FALSE,
                                          RANDSTART=1L),
                              DERIV=DERIV)
            ESTIMATE$pars <- tmp$pars
            ESTIMATE$random <- tmp$random
            ESTIMATE$lrPars <- tmp$lrPars
            ESTIMATE$info <- tmp$info
            ESTIMATE$fail_invert_info <- tmp$fail_invert_info
            ESTIMATE$time <- c(ESTIMATE$time, SE=sum(tmp$time))
        }
        rlist <- vector('list', Data$ngroups)
        for(g in 1L:Data$ngroups)
            rlist[[g]]$expected = numeric(1L)
    }
    opts$times$end.time.Estimate <- proc.time()[3L]
    opts$times$start.time.SE <- proc.time()[3L]
    if(!opts$NULL.MODEL && opts$SE){
        tmp <- ESTIMATE
        if(opts$verbose && !(opts$method == 'MHRM' || opts$method == 'MIXED'))
            cat('\n\nCalculating information matrix...\n')
        if(opts$SE.type == 'complete' && opts$method == 'EM'){
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
                                          itemloc=PrepList[[1L]]$itemloc, BFACTOR=opts$BFACTOR,
                                          sitems=sitems, specific=specific, NULL.MODEL=opts$NULL.MODEL,
                                          nfact=nfact, constrain=constrain, verbose=opts$verbose,
                                          CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND, Moptim=ESTIMATE$Moptim,
                                          EH=opts$empiricalhist, EHPrior=ESTIMATE$Prior, warn=opts$warn,
                                          message=opts$message, full=opts$full, lrPars=lrPars),
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
                if(opts$technical$symmetric_SEM) info <- (info + t(info)) / 2
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
        } else if(opts$SE.type %in% c('Richardson', 'forward', 'central') && opts$method != 'MIXED'){
            ESTIMATE <- SE.Numerical(pars=ESTIMATE$pars, Theta=Theta, theta=theta, PrepList=PrepList, Data=Data,
                              BFACTOR=opts$BFACTOR, itemloc=PrepList[[1L]]$itemloc, ESTIMATE=ESTIMATE,
                              constrain=constrain, Ls=Ls, specific=oldmodel, sitems=sitems, EH=opts$empiricalhist,
                              CUSTOM.IND=CUSTOM.IND, EHPrior=ESTIMATE$Prior, warn=opts$warn, type=opts$SE.type,
                              delta=opts$delta, lrPars=ESTIMATE$lrPars)
        } else if(opts$SE.type == 'MHRM' && opts$method == 'EM'){
            if(opts$empiricalhist)
                stop('MHRM standard error not available when using empirical histograms', call.=FALSE)
            ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, Ls=Ls, PrepList=PrepList, Data=Data,
                                   list = list(NCYCLES=1000L, BURNIN=1L, SEMCYCLES=opts$SEMCYCLES,
                                               KDRAWS=opts$KDRAWS, MHDRAWS=opts$MHDRAWS,
                                               TOL=opts$SEtol, SE=TRUE, SE.type=opts$SE.type,
                                               gain=opts$gain, nfactNames=PrepList[[1L]]$nfactNames,
                                               itemloc=PrepList[[1L]]$itemloc, BFACTOR=opts$BFACTOR,
                                               nfact=nfact, constrain=constrain, verbose=FALSE, expl=FALSE,
                                               CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND, message=opts$message,
                                               startlongpars=startlongpars, SE=opts$SE, warn=opts$warn),
                                   DERIV=DERIV)
        } else if(any(opts$SE.type %in% c('crossprod', 'Louis', 'sandwich')) && opts$method != 'MIXED'){
            if(logPrior != 0 && opts$warn)
                warning('Information matrix with the crossprod, Louis, and sandwich method
                        do not account for prior parameter distribution information')
            ESTIMATE <- SE.simple(PrepList=PrepList, ESTIMATE=ESTIMATE, Theta=Theta, Data=Data,
                                  constrain=constrain, Ls=Ls, N=nrow(data), type=opts$SE.type,
                                  CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND, warn=opts$warn,
                                  message=opts$message)
        } else if(opts$SE.type == 'Fisher' && opts$method != 'MIXED'){
            if(logPrior != 0 && opts$warn)
                warning('Information matrix with the Fisher method does not
                        account for prior parameter distribution information')
            ESTIMATE <- SE.Fisher(PrepList=PrepList, ESTIMATE=ESTIMATE, Theta=Theta, Data=Data,
                                  constrain=constrain, Ls=Ls, N=nrow(data), full=opts$full,
                                  CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND, warn=opts$warn)
        }
        ESTIMATE$cycles <- tmp$cycles
        ESTIMATE$Prior <- tmp$Prior
    }
    opts$times$end.time.SE <- proc.time()[3L]
    opts$times$start.time.post <- proc.time()[3L]
    cmods <- vector('list', Data$ngroups)
    if(is.null(opts$theta_lim))
        opts$theta_lim <- numeric(1)
    lrPars <- ESTIMATE$lrPars
    class(lrPars) <- 'S4'
    for(g in 1L:Data$ngroups){
        if(opts$method == 'MIXED'){
            F <- matrix(NA)
            h2 <- numeric(1)
        } else {
            F <- Lambdas(ESTIMATE$pars[[g]], Names=colnames(data))
            colnames(F) <- PrepList[[1L]]$factorNames
            h2 <- rowSums(F^2)
        }
        cmods[[g]] <- new('SingleGroupClass', ParObjects=list(pars=ESTIMATE$pars[[g]], lrPars=lrPars,
                                                              random=ESTIMATE$random),
                          Data = list(K=Data$K, nitems=Data$nitems),
                          Model = list(itemloc=PrepList[[1L]]$itemloc, nfact=nfact, constrain=constrain,
                                       factorNames=PrepList[[1L]]$factorNames,
                                       itemtype=PrepList[[1L]]$itemtype,
                                       prodlist=PrepList[[1L]]$prodlist),
                          Options = list(method = 'MHRM', exploratory=PrepList[[1L]]$exploratory,
                                         theta_lim=opts$theta_lim),
                          Fit = list(G2=G2group[g], F=F, h2=h2),
                          Internals = list(Pl = rlist[[g]]$expected, CUSTOM.IND=CUSTOM.IND,
                                           SLOW.IND=SLOW.IND))
        if(discrete){
            cmods[[g]]@Model$Theta <- Theta
            cmods[[g]]@Internals$Prior <- list(ESTIMATE$Prior[[g]])
        }
    }
    #missing stats for MHRM
    if(opts$method =='MHRM' || opts$method == 'MIXED'){
        if(opts$verbose) cat("\nCalculating log-likelihood...\n")
        flush.console()
        logLik <- G2 <- SElogLik <- 0
        Pl <- list()
        if(!opts$technical$parallel){
            ncores <- .mirtClusterEnv$ncores
            .mirtClusterEnv$ncores <- 1L
        }
        logPrior <- 0
        for(g in 1L:Data$ngroups){
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

    ####post estimation stats
    if(opts$Moptim %in% c('solnp', 'alabama')){
        if(!is.null(opts$solnp_args$eqfun))
            df <- df + length(opts$solnp_args$eqfun(ESTIMATE$shortpars, list()))
        if(!is.null(opts$solnp_args$heq))
            df <- df + length(opts$solnp_args$heq(ESTIMATE$shortpars, list()))
    }
    r <- rr
    N <- sum(r)
    tmp <- dfsubtr
    if(logPrior != 0){
        AIC <- AICc <- BIC <- SABIC <- NaN
        DIC <- -2 * (logLik + logPrior) + 2 * tmp
    } else {
        AIC <- (-2) * logLik + 2 * tmp
        AICc <- AIC + 2 * tmp * (tmp + 1) / (N - tmp - 1L)
        BIC <- (-2) * logLik + tmp*log(N)
        SABIC <- (-2) * logLik + tmp*log((N+2)/24)
        DIC <- AIC
    }
    p.G2 <- 1 - pchisq(G2,df)
    RMSEA.G2 <- rmsea(X2=G2, df=df, N=N)
    null.mod <- unclass(new('SingleGroupClass'))
    TLI.G2 <- CFI.G2 <- NaN
    if(length(r) * 3L < prod(Data$K)){
        G2 <- NaN; p.G2 <- NaN
        opts$calcNull <- FALSE
    }
    if(!opts$NULL.MODEL && opts$method != 'MIXED' && opts$calcNull && nmissingtabdata == 0L){
        null.mod <- try(unclass(computeNullModel(data=data, itemtype=itemtype,
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
                AIC=AIC, AICc=AICc, BIC=BIC, SABIC=SABIC, DIC=DIC, logLik=logLik,
                logPrior=logPrior, SElogLik=SElogLik, F=F, h2=h2)
    Model <- list(model=oldmodel, factorNames=PrepList[[1L]]$factorNames, itemtype=PrepList[[1L]]$itemtype,
                  itemloc=PrepList[[1L]]$itemloc, nfact=nfact,
                  Theta=Theta, constrain=constrain, parprior=parprior, nest=as.integer(dfsubtr),
                  invariance=invariance, lrPars=lrPars, formulas=attr(mixed.design, 'formula'),
                  prodlist=PrepList[[1L]]$prodlist)
    Data$covdata <- if(length(lrPars)) lrPars@df else attr(mixed.design, 'covdata')
    Data$itemdesign <- attr(mixed.design, 'itemdesign')
    ParObjects <- list(pars=cmods, lrPars=lrPars, random=ESTIMATE$random)
    OptimInfo <- list(iter=ESTIMATE$cycles, converged=ESTIMATE$converge, cand.t.var=ESTIMATE$cand.t.var,
                      condnum=NA, secondordertest=NA, SEMconv=SEMconv)
    vcov <- matrix(NA, 1, 1)
    if(Options$SE){
        information <- ESTIMATE$info
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
        } else OptimInfo$secondordertest <- FALSE
    }
    Internals <- list(collectLL=ESTIMATE$collectLL, Prior=ESTIMATE$Prior, Pl=Pl,
                      shortpars=as.numeric(ESTIMATE$shortpars), key=key,
                      bfactor=list(), CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND)
    if(discrete){
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
            mod <- new('MultipleGroupClass',
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
    mod@time <- c(TOTAL = as.numeric(proc.time()[3L] - time$start.time),
                  DATA = as.numeric(time$end.time.Data - time$start.time.Data),
                  ESTIMATE = as.numeric(time$end.time.Estimate - time$start.time.Estimate),
                  ESTIMATE$time,
                  SE = as.numeric(time$end.time.SE - time$start.time.SE),
                  POST = as.numeric(time$end.time.post - time$start.time.post))
    return(mod)
}
