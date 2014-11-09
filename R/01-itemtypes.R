#Classes

setClass("GroupPars",
         representation(par='numeric',
                        SEpar='numeric',
                        est='logical',
                        parnum='numeric',
                        nfact='integer',
                        gradient='numeric',
                        hessian='matrix',
                        itemtrace='matrix',
                        lbound='numeric',
                        ubound='numeric',
                        any.prior='logical',
                        prior.type='integer',
                        prior_1='numeric',
                        prior_2='numeric')
)

setClass("RandomPars",
         representation(par='numeric',
                        SEpar='numeric',
                        est='logical',
                        between='logical',
                        parnum='numeric',
                        ndim='integer',
                        gframe='data.frame',
                        gdesign='matrix',
                        lbound='numeric',
                        ubound='numeric',
                        cand.t.var='numeric',
                        drawvals='matrix',
                        mtch='numeric',
                        any.prior='logical',
                        prior.type='integer',
                        prior_1='numeric',
                        prior_2='numeric')
)

setClass("lrPars",
         representation(par='numeric',
                        SEpar='numeric',
                        est='logical',
                        beta='matrix',
                        sigma = 'matrix',
                        mus='matrix',
                        df='data.frame',
                        parnum='numeric',
                        nfact='integer',
                        nfixed='integer',
                        X='matrix',
                        tXX='matrix',
                        inv_tXX='matrix',
                        lbound='numeric',
                        ubound='numeric',
                        any.prior='logical',
                        prior.type='integer',
                        prior_1='numeric',
                        prior_2='numeric',
                        formula='formula',
                        EM='logical')
)

setClass("AllItemsClass",
         representation(par='numeric',
                        SEpar='numeric',
                        est='logical',
                        constr='logical',
                        itemclass='integer',
                        parnum='numeric',
                        nfact='integer',
                        nfixedeffects='numeric', #number of fixed effect predictors
                        fixed.design='matrix',
                        dat='matrix',
                        ncat='integer',
                        gradient='numeric',
                        hessian='matrix',
                        itemtrace='matrix',
                        lbound='numeric',
                        ubound='numeric',
                        any.prior='logical',
                        prior.type='integer',
                        prior_1='numeric',
                        prior_2='numeric',
                        'VIRTUAL')
)

setClass("dich", contains = 'AllItemsClass')

setClass("graded", contains = 'AllItemsClass')

setClass("rating", contains = 'AllItemsClass')

setClass("gpcm", contains = 'AllItemsClass')

setClass("rsm", contains = 'AllItemsClass')

setClass("nominal", contains = 'AllItemsClass')

setClass("partcomp", contains = 'AllItemsClass')

setClass("nestlogit", contains = 'AllItemsClass',
         representation = representation(correctcat='integer'))

setClass("ideal", contains = 'AllItemsClass')

setClass("lca", contains = 'AllItemsClass',
         representation = representation(score='numeric'))

setClass('custom', contains = 'AllItemsClass',
         representation = representation(name='character',
                                         P='function',
                                         gr='function',
                                         usegr='logical',
                                         hss='function',
                                         usehss='logical',
                                         userdata='matrix',
                                         useuserdata='logical'))

#--------------------------------------------------------------------------

#Generics

setGeneric('ProbTrace', function(x, Theta, ...) standardGeneric("ProbTrace"))

setGeneric('ExtractLambdas', function(x) standardGeneric("ExtractLambdas"))

setGeneric('ExtractZetas', function(x) standardGeneric("ExtractZetas"))

setGeneric('Deriv', function(x, Theta, ...) standardGeneric("Deriv"))

setGeneric('DerivTheta', function(x, Theta) standardGeneric("DerivTheta"))

setGeneric('dP', function(x, Theta) standardGeneric("dP"))

setGeneric('calcLogLik', function(object, ...) standardGeneric("calcLogLik"))

setGeneric("itemplot.internal",  function(object, ...) standardGeneric("itemplot.internal"))

setGeneric("fscores.internal", function(object, ...) standardGeneric("fscores.internal"))

setGeneric("summary", function(object, ...) standardGeneric("summary"))

setGeneric('DrawValues', function(x, Theta, ...) standardGeneric("DrawValues"))

setGeneric('RandomDeriv', function(x, ...) standardGeneric("RandomDeriv"))

setGeneric('GenRandomPars', function(x) standardGeneric("GenRandomPars"))
