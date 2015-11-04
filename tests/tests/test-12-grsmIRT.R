# Response data : num of response categories = 5

R <-
  structure(c(4, 3, 3, 2, 2, 3, 2, 3, 4, 4, 3, 3, 2, 3, 2, 5, 4,
              2, 2, 4, 2, 1, 4, 3, 4, 1, 4, 4, 4, 4, 4, 2, 1, 4, 2, 2, 1, 1,
              3, 5, 4, 2, 2, 2, 3, 3, 5, 2, 2, 2, 2, 2, 5, 4, 4, 4, 2, 4, 3,
              4, 4, 2, 3, 4, 4, 3, 4, 4, 4, 2, 4, 5, 4, 4, 1, 4, 3, 5, 4, 3,
              3, 3, 4, 2, 2, 4, 5, 5, 2, 3, 2, 4, 1, 5, 4, 5, 1, 2, 3, 2, 5,
              3, 3, 3, 4, 4, 2, 3, 4, 4, 4, 2, 3, 3, 4, 2, 4, 3, 3, 2, 4, 5,
              1, 3, 1, 4, 4, 4, 4, 2, 1, 3, 4, 2, 1, 3, 5, 5, 2, 4, 1, 4, 2,
              4, 4, 2, 1, 4, 3, 2, 4, 3, 4, 4, 4, 3, 2, 2, 4, 3, 3, 2, 4, 3,
              4, 3, 4, 4, 3, 2, 4, 5, 2, 3, 1, 5, 4, 2, 4, 1, 2, 3, 3, 3, 1,
              4, 5, 5, 2, 4, 2, 5, 2, 4, 4, 2, 1, 4, 3, 1, 5, 4, 1, 3, 5, 2,
              2, 2, 2, 4, 3, 2, 3, 3, 4, 1, 4, 4, 4, 3, 4, 4, 2, 3, 1, 5, 2,
              5, 4, 2, 2, 2, 2, 3, 1, 4, 4, 4, 2, NA, 1, 4, 2, 2, 3, 1, 1,
              4, 3, 2), .Dim = c(50L, 5L),
            .Dimnames = list(c("3", "4", "5",
                               "6", "7", "10", "11", "12", "17", "18", "19", "25", "26", "28",
                               "29", "33", "34", "35", "39", "40", "42", "44", "45", "46", "47",
                               "48", "50", "53", "56", "61", "63", "66", "67", "68", "71", "72",
                               "73", "77", "80", "81", "82", "83", "84", "85", "87", "88", "89",
                               "90", "91", "92"), c("TA1-1", "TA1-2", "TA1-3", "TA1-4", "TA1-5"
                               )))

coef2mat = function(cf) {
  stopifnot(is.list(cf))
  do.call(rbind, cf[1:(length(cf)-1)])
}

Qncat = function(R) {
  stopifnot(is.matrix(R))
  max(apply(R, 2, function(x) length(unique(x[!is.na(x)]))))
}

n.item <- ncol(R); ncat <- Qncat(R)

model <- paste("F1 = 1-",n.item,
               "\nCONSTRAIN = (1-",n.item,",d1)",
               paste(",(1-",n.item,",d",2:(ncat-1),")",sep="", collapse=""),
               "\nSTART = (1, c, 0.0)",
               "\nFIXED = (1, c)",
               sep=""); cat(model)

fit.grsmIRT <- mirt(R, mirt.model(model), itemtype="grsmIRT")
coef2mat(coef(fit.grsmIRT))

R[R==5] = 4 # when the num. of response categories is 4
n.item <- ncol(R); ncat <- Qncat(R)
model <- paste("F1 = 1-",n.item,
               "\nCONSTRAIN = (1-",n.item,",d1)",
               paste(",(1-",n.item,",d",2:(ncat-1),")",sep="", collapse=""),
               "\nSTART = (1, c, 0.0)",
               "\nFIXED = (1, c)",
               sep=""); cat(model)
fit.grsmIRT <- mirt(R, mirt.model(model), itemtype="grsmIRT")
coef2mat(coef(fit.grsmIRT))


R[R==1] = 2 # when the num. of response categories is 3
n.item <- ncol(R); ncat <- Qncat(R)
model <- paste("F1 = 1-",n.item,
               "\nCONSTRAIN = (1-",n.item,",d1)",
               paste(",(1-",n.item,",d",2:(ncat-1),")",sep="", collapse=""),
               "\nSTART = (1, c, 0.0)",
               "\nFIXED = (1, c)",
               sep=""); cat(model)
fit.grsmIRT <- mirt(R, mirt.model(model), itemtype="grsmIRT")
coef2mat(coef(fit.grsmIRT))

