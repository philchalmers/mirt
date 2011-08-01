\name{confmirt.model}
\alias{confmirt.model}
\title{Specify confirmatory loadings}
\description{
The \code{confmirt.model} function scans user input to specify the confirmatory model.
}

\usage{
confmirt.model(file = '')
}

\arguments{
  \item{file}{a string specifying an external file that declares the input.}
}  

\value{
  Returns a model specification object to be used in \code{\link{confmirt}}.
}

\details{
Factors are named and specify which numerical items they affect, seperated either by commas or by - to indicate a range of items. Four optional keyword to specify relationships between factors or by specifying constraints on the model. Note that all independent terms must be seperated with a comma. These are:
  
\describe{
  \item{COV}{Specify the realatinship between the latent traits. Estimating a correalation is declared by pairing two factors with an asterik (e.g., F1*F2), while fixing this value to a constant is performed by stating 'eq' and some numeric value (e.g., F1*F2 eq .5).}
  \item{MEAN}{Specify the constant values of the latent intercepts. For instance, F1 eq 1, gives a latent mean of 1 for factor F1. }
  \item{SLOPE}{Specify fixed and equal relationsionships between the item slopes. These constraints require the declaration of the factor, followed by @, then the respective item, and finally the constraint behavior. Paramters may be fixed to be equal using 'eq', or equal to a constant value using the same method. For instance, F1@2 eq .5, fixes the first factor slope on item 2 to be .5, and F1@2 eq  F2@4 eq F2@5, fixes the three slope parameters to be equal during estimation.}
  \item{INT}{Specify fixed and equal relationsionships between the item intercepts. These constraints require the declaration of the item, followed by @, then the respective intercept (from lowest to highest), and finally the constraint behavior. Paramters may be fixed to be equal using 'eq', or equal to a constant value using the same method. For instance, 1@2 eq .5, fixes the first item's second intercept to be .5, and 1@1 eq 2@1 eq 3@1, fixes the three intercept parameters to be equal during estimation.}  
  }

}

\examples{
\dontrun{
model <- confmirt.model()
  F1  = 1,2,3,4-10
  F2 = 11-20
  COV = F1*F2             
  SLOPE = F1@1 eq F1@2   
  MEAN = F1 eq 1, F2 eq -1          
  INT = 1@1 eq -0.5
    
    }
}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}


