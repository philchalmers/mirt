#' Compute numerical derivatives
#'
#' Compute numerical derivatives using forward/backword difference,
#' central difference, or Richardson extropolation.
#'
#' @param f the objective function being evaluated
#' @param par a vector of parameters
#' @param ... additional arguments to be passed to \code{f}
#' @param delta the term used to perturb the \code{f} function. Default is 1e-5
#' @param gradient logical; compute the gradient terms? If FALSE then the Hessian is computed instead
#' @param type type of difference to compute. Can be either \code{'forward'} for the forward difference,
#'   \code{'central'} for the central difference (default), or \code{'Richardson'} for the Richardson extropolation.
#'   Backword difference is acheived by supplying a negative \code{delta} value with \code{'forward'}.
#'   When \code{type = 'Richardson'}, the default value of \code{delta} is increased to \code{delta * 100}
#'   to provide a reasonable perterbation starting location (each \code{delta} is halved at each iteration).
#' @export numerical_deriv
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords numerical derivatives
#'
#' @examples
#'
#' \dontrun{
#' f <- function(x) 3*x[1]^3 - 4*x[2]^2
#' par <- c(3,8)
#'
#' # grad = 9 * x^2 , -8 * y
#' (actual <- c(9 * par[1]^2, -8 * par[2]))
#' numerical_deriv(f, par, type = 'forward')
#' numerical_deriv(f, par, type = 'central')
#' numerical_deriv(f, par, type = 'Richardson')
#'
#' # hessian = h11 -> 18 * x, h22 -> -8, h12 -> h21 -> 0
#' (actual <- matrix(c(18 * par[1], 0, 0, -8), 2, 2))
#' numerical_deriv(f, par, type = 'forward', gradient = FALSE)
#' numerical_deriv(f, par, type = 'central', gradient = FALSE)
#' numerical_deriv(f, par, type = 'Richardson', gradient = FALSE)
#'
#' }
numerical_deriv <- function(f, par, ...,  delta = 1e-5, gradient = TRUE, type = 'central'){
    forward_difference <- function(par, f, delta, ...){
        dots <- list(...)
        np <- length(par)
        g <- numeric(np)
        if(is.null(dots$ObJeCtIvE)) fx <- f(par, ...) else fx <- dots$ObJeCtIvE
        for(i in seq_len(np)){
            p <- par
            p[i] <- p[i] + delta
            g[i] <- (f(p, ...) - fx) / delta
        }
        g
    }
    forward_difference2 <- function(par, f, delta, ...){
        dots <- list(...)
        np <- length(par)
        hess <- matrix(0, np, np)
        if(is.null(dots$ObJeCtIvE)) fx <- f(par, ...) else fx <- dots$ObJeCtIvE
        fx1 <- numeric(np)
        for(i in seq_len(np)){
            tmp <- par
            tmp[i] <- tmp[i] + delta
            fx1[i] <- f(tmp, ...)
        }
        for(i in seq_len(np)){
            for(j in i:np){
                fx1x2 <- par
                fx1x2[i] <- fx1x2[i] + delta
                fx1x2[j] <- fx1x2[j] + delta
                hess[i,j] <- hess[j, i] <- (f(fx1x2, ...) - fx1[i] - fx1[j] + fx) / (delta^2)
            }
        }
        (hess + t(hess))/2
    }
    central_difference <- function(par, f, delta, ...){
        np <- length(par)
        g <- numeric(np)
        for(i in seq_len(np)){
            p1 <- p2 <- par
            p1[i] <- p1[i] + delta
            p2[i] <- p2[i] - delta
            g[i] <- (f(p1, ...) - f(p2, ...)) / (2 * delta)
        }
        g
    }
    central_difference2 <- function(par, f, delta, ...){
        np <- length(par)
        hess <- matrix(0, np, np)
        fx <- f(par, ...)
        for(i in seq_len(np)){
            for(j in i:np){
                if(i == j){
                    p <- par
                    p[i] <- p[i] + 2 * delta; s1 <- f(p, ...)
                    p[i] <- p[i] - 4 * delta; s3 <- f(p, ...)
                    hess[i, i] <- (s1 - 2*fx + s3) / (4 * delta^2)
                } else {
                    p <- par
                    p[i] <- p[i] + delta; p[j] <- p[j] + delta; s1 <- f(p, ...)
                    p[j] <- p[j] - 2*delta; s2 <- f(p, ...)
                    p[i] <- p[i] - 2*delta; s4 <- f(p, ...)
                    p[j] <- p[j] + 2*delta; s3 <- f(p, ...)
                    hess[i,j] <- hess[j,i] <- (s1 - s2 - s3 + s4) / (4 * delta^2)
                }
            }
        }
        (hess + t(hess))/2
    }
    richardson <- function(par, f, delta, r = 4L, eps = 1e-12, ...){
        R0 <- R1 <- matrix(0, length(par), r)
        R0[, 1L] <- central_difference(par=par, f=f, delta=delta, ...)
        for(i in 1L:(r-1L)){
            delta <- delta/2
            R1[ ,1L] <- central_difference(par=par, f=f, delta=delta, ...)
            for (j in 1L:i)
                R1[ ,j + 1] <- (4^j * R1[ , j] - R0[, j]) / (4^j - 1)
            if(all(abs(R1[ ,j + 1] - R0[ ,j]) < eps)) break
            R0 <- R1
        }
        R1[ , i+1]
    }
    richardson2 <- function(par, f, delta, r = 4L, eps = 1e-12, ...){
        R0 <- R1 <- matrix(0, length(par)^2, r)
        R0[, 1L] <- as.vector(central_difference2(par=par, f=f, delta=delta, ...))
        for(i in 1L:(r-1L)){
            delta <- delta/2
            R1[ ,1L] <- as.vector(central_difference2(par=par, f=f, delta=delta, ...))
            for (j in 1L:i)
                R1[ ,j + 1] <- (4^j * R1[ , j] - R0[, j]) / (4^j - 1)
            if(all(abs(R1[ ,j + 1] - R0[ ,j]) < eps)) break
            R0 <- R1
        }
        hess <- matrix(R1[ , i+1], length(par), length(par))
        (hess + t(hess))/2
    }

    if(!length(par)){
        if(gradient) return(numeric())
        else return(matrix(numeric()))
    }
    if(type == 'central'){
        ret <- if(gradient) central_difference(par=par, f=f, delta=delta, ...)
        else central_difference2(par=par, f=f, delta=delta, ...)
    } else if(type == 'forward'){
        ret <- if(gradient) forward_difference(par=par, f=f, delta=delta, ...)
        else forward_difference2(par=par, f=f, delta=delta, ...)
    } else if(type == 'Richardson'){
        ret <- if(gradient) richardson(par=par, f=f, delta=delta*100, ...)
        else richardson2(par=par, f=f, delta=delta*100, ...)
    }
    ret
}