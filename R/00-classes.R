setClass("AllModelClass",
         representation(pars='list', 
                        K='numeric', 
                        G2='numeric', 
                        df='numeric', 
                        p='numeric', 
                        AIC='numeric', 
                        F='matrix', 
                        h2='numeric', 
                        tabdata='matrix', 
                        tabdatalong='matrix', 
                        Theta='matrix', 
                        Pl='numeric',
                        data='matrix',                         
                        converge='numeric', 
                        itemloc = 'numeric',                        
                        BIC='numeric', 
                        RMSEA='numeric',                         
                        null.mod = 'S4', 
                        TLI = 'numeric',                          
                        logLik='numeric',                                    
                        SElogLik='numeric',
                        Call='call',
                        esttype='character',
                        iter='numeric',                        
                        quadpts='numeric',
                        nfact='numeric',                        
                        prodlist='list',
                        constrain='list',
                        fulldata='matrix',
                        factorNames='character',
                        'VIRTUAL'),    
         validity = function(object) return(TRUE)
)                       

#------------------------------------------------------------------------------
#' @exportClass ExploratoryClass
# @keywords classes
setClass(
    Class = 'ExploratoryClass', contains = 'AllModelClass',
    representation = representation(Target='numeric',
                                    rotate='character'),	
    validity = function(object) return(TRUE)
)	

#------------------------------------------------------------------------------
#' @exportClass ConfirmatoryClass
# @keywords classes
setClass(
    Class = 'ConfirmatoryClass', contains = 'AllModelClass',    	
    validity = function(object) return(TRUE)
)	

#------------------------------------------------------------------------------
#' @exportClass MultipleGroupClass
setClass(
    Class = 'MultipleGroupClass', contains = 'AllModelClass',
    representation = representation(group='factor',
                                    groupNames='factor',
                                    invariance='character',
                                    cmods='list'),    
    validity = function(object) return(TRUE)
)
