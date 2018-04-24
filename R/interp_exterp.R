interpolateWoods <- function (point, std_point, nq1, nq2, delta) {
    term1 <- (point - std_point)/delta
    term2 <- nq2 - nq1

    term1 * term2 + nq1
}

extrapolateWoods <- function (point, std_point, nq1, nq2, delta, tail) {
    if (tail == "left") {
        ratio <- nq1/nq2
        if(is.nan(ratio) || ratio < .0000001) ratio <- .001
        (ratio ^ ((std_point - point) / delta)) * nq1
    } else if (tail == "right") {
        ratio <- nq2/nq1
        if(is.nan(ratio) || ratio < .0000001) ratio <- .001
        (ratio ^ ((point - std_point) / delta)) * nq2
    }
}

standardizeQuadrature <- function (qp, nq, estmean=FALSE, estsd=FALSE) {

    # Standardize qp:
    qp_m <- sum(nq*qp) / sum(nq)
    qp_sd <- sqrt(sum(nq * (qp - qp_m)^2) / sum(nq))
    attr_ret <- c(mean = qp_m, var=qp_sd^2)
    if(estmean) qp_m <- 0
    if(estsd) qp_sd <- 1
    std_qp <- (qp - qp_m) / qp_sd

    if(estmean && estsd){
        attr(nq, 'mean_var') <- attr_ret
        return(nq)
    }

    min_stdqp <- min(std_qp)
    max_stdqp <- max(std_qp)

    res <- numeric(length(qp))
    delta <- qp[2] - qp[1] # is the distance between any two q or, equivalently, between any two q*

    for (i in 1:length(qp)) {
        if (qp[i] <= min_stdqp) {
            res[i] <- extrapolateWoods(qp[i], min_stdqp, nq[1], nq[2], delta, "left")
        } else if (max_stdqp <= qp[i]) {
            res[i] <- extrapolateWoods(qp[i], max_stdqp, nq[length(nq)-1], nq[length(nq)], delta, "right")
        } else {
            std_ind <- max(which(qp[i] > std_qp))
            res[i] <- interpolateWoods(qp[i], std_qp[std_ind], nq[std_ind], nq[std_ind+1], delta)
        }
    }

    res <- res / sum(res) * sum(nq)
    attr(res, 'mean_var') <- attr_ret
    res
}
