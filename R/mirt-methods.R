#Methods 
setMethod(
    f = "print",
    signature = signature(x = 'mirtClass'),
    definition = function(x){  
        printFunction(x=x)		
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'mirtClass'),
    definition = function(object){  
        printFunction(x=object)
    }
)

setMethod(
    f = "summary",
    signature = 'mirtClass',
    definition = function(object, rotate = '', Target = NULL, suppress = 0, digits = 3, 
                          print = TRUE, ...){        
        summaryFunction(object=object, rotate=rotate, Target=Target, suppress=suppress, 
                        digits=digits, print=print, ...)
    }
)

setMethod(
    f = "coef",
    signature = 'mirtClass',
    definition = function(object, rotate = '', Target = NULL, SE = TRUE, digits = 3, ...){  
        coefFunction(object=object, rotate=rotate, Target=Target, SE=SE, digits=digits, ...)
    }
)

setMethod(
    f = "anova",
    signature = signature(object = 'mirtClass'),
    definition = function(object, object2){        
        anovaFunction(object=object, object2=object2)
    }
)

setMethod(
    f = "residuals",
    signature = signature(object = 'mirtClass'),
    definition = function(object, restype = 'LD', digits = 3, printvalue = NULL)
    {   
        residualsFunction(object=object, restype=restype, digits=digits, printvalue=printvalue)
        
    }
)

setMethod(
    f = "plot",
    signature = signature(x = 'mirtClass', y = 'missing'),
    definition = function(x, y, type = 'info', npts = 50, theta_angle = 45, 
                          rot = list(xaxis = -70, yaxis = 30, zaxis = 10))
    {          
        plotFunction(x=x, y=y, type=type, npts=npts, theta_angle=theta_angle, 
                     rot=rot)
    }		
)	

setMethod(
    f = "fitted",
    signature = signature(object = 'mirtClass'),
    definition = function(object, digits = 3){  
        fittedFunction(object=object, digits=digits)
    }
)
