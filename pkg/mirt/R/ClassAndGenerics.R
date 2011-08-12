########################################## 
#Class defs

setClass(
	Class = 'mirtClass',
	representation = representation(EMiter = 'numeric', pars = 'matrix', guess = 'numeric', 
		X2 = 'numeric', df = 'numeric', p = 'numeric', AIC = 'numeric', log.lik = 'numeric',
		F = 'matrix', h2 = 'numeric', tabdata = 'matrix', Theta = 'matrix', Pl = 'numeric',
		fulldata = 'matrix', cormat = 'matrix', facility = 'numeric', converge = 'numeric', 
		quadpts = 'numeric', BIC = 'numeric', Call = 'call'),	
	validity = function(object) return(TRUE)
)	

setClass(
	Class = 'bfactorClass',
	representation = representation(EMiter = 'numeric', pars = 'matrix', guess = 'numeric', 
		AIC = 'numeric', X2 = 'numeric', df = 'numeric', log.lik = 'numeric', p = 'numeric', 
		F = 'matrix', h2 = 'numeric', itemnames = 'character', tabdata = 'matrix', 
		N = 'numeric', Pl = 'numeric', Theta = 'matrix', fulldata = 'matrix', 
		logicalfact = 'matrix', facility = 'numeric', specific = 'numeric', BIC = 'numeric',
		cormat = 'matrix', converge = 'numeric', par.prior = 'matrix', quadpts = 'numeric', 
		Call = 'call'),	
	validity = function(object) return(TRUE)
)	

setClass(
	Class = 'polymirtClass',
	representation = representation(pars = 'matrix', guess = 'numeric', SEpars = 'matrix', 
		cycles = 'numeric', Theta = 'matrix', fulldata = 'matrix', data = 'matrix', 
		K = 'numeric', F = 'matrix', h2 = 'numeric', itemloc = 'numeric', AIC = 'numeric',
		converge = 'numeric', logLik = 'numeric', SElogLik = 'numeric', df = 'integer', 
		G2 = 'numeric', p = 'numeric', tabdata = 'matrix', BIC = 'numeric', Call = 'call'),	
	validity = function(object) return(TRUE)
)	

setClass(
	Class = 'confmirtClass',
	representation = representation(pars = 'matrix', guess = 'numeric', SEpars = 'matrix', 
		SEg = 'numeric', gpars = 'list', SEgpars = 'list', estpars = 'list',cycles = 'numeric', 
		Theta = 'matrix', fulldata = 'matrix', data = 'matrix', K = 'numeric', itemloc = 'numeric',
		h2 = 'numeric',F = 'matrix', converge = 'numeric', logLik = 'numeric',SElogLik = 'numeric',
		df = 'integer', AIC = 'numeric', nconstvalues = 'integer', G2 = 'numeric', p = 'numeric',
		tabdata = 'matrix', BIC = 'numeric', estComp = 'logical', Call = 'call'),	
	validity = function(object) return(TRUE)
)	

########################################## 
#Generic defs

setGeneric("fscores", 
	def = function(object, ...) standardGeneric("fscores")
)

setGeneric("itemplot", 
	def = function(object, item, ...) standardGeneric("itemplot")
)

