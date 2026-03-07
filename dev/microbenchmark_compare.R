#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(microbenchmark))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (length(hit) == 0L) return(default)
  sub(key, "", hit[[1L]], fixed = TRUE)
}

arg_int <- function(name, default) {
  as.integer(arg_value(name, as.character(default)))
}

mode <- arg_value("mode", "baseline")
repo <- normalizePath(arg_value("repo", "."), mustWork = TRUE)
out <- arg_value("out", "")
lib_loc <- arg_value("lib", "")
dataset <- arg_value("dataset", "sim")

N_prepare <- arg_int("n_prepare", 50000L)
N_fit <- arg_int("n_fit", 12000L)
J <- arg_int("j", 40L)
theta_n <- arg_int("theta_n", 20000L)
omp_threads <- arg_int("omp_threads", 1L)
times_prepare <- arg_int("times_prepare", 5L)
times_fit <- arg_int("times_fit", 3L)
times_eapsum <- arg_int("times_eapsum", 5L)
times_itemtrace <- arg_int("times_itemtrace", 5L)

if (!mode %in% c("baseline", "head")) {
  stop("--mode must be one of: baseline, head")
}

if (!dataset %in% c("sim", "sat12")) {
  stop("--dataset must be one of: sim, sat12")
}

if (mode == "baseline") {
  if (nzchar(lib_loc)) {
    suppressPackageStartupMessages(library("mirt", lib.loc = lib_loc, character.only = TRUE))
  } else {
    suppressPackageStartupMessages(library(mirt))
  }
} else {
  if (nzchar(lib_loc)) {
    suppressPackageStartupMessages(library("mirt", lib.loc = lib_loc, character.only = TRUE))
  } else {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("pkgload is required for --mode=head without --lib")
    }
    pkgload::load_all(repo, compile = FALSE, quiet = TRUE)
  }
}

set.seed(20260227)
if (dataset == "sim") {
  a <- matrix(abs(rnorm(J, mean = 1.1, sd = 0.2)), ncol = 1L)
  d <- rnorm(J, mean = -0.1, sd = 0.9)
  dat_prepare <- simdata(a, d, N = N_prepare, itemtype = "2PL")
  dat_fit <- simdata(a, d, N = N_fit, itemtype = "2PL")
  itemtype <- "2PL"
} else {
  data("SAT12", package = "mirt", envir = environment())
  key <- c(1, 4, 5, 2, 3, 1, 2, 1, 3, 1, 2, 4, 2, 1, 5, 3,
           4, 4, 1, 4, 3, 3, 4, 1, 3, 5, 1, 3, 1, 5, 4, 5)
  SAT12[SAT12 == 8] <- NA
  dat_prepare <- dat_fit <- key2binary(SAT12, key)
  itemtype <- NULL
  N_prepare <- N_fit <- nrow(dat_fit)
  J <- ncol(dat_fit)
}

bench_prepare <- microbenchmark(
  prepare_large = mirt(dat_prepare, 1, itemtype = itemtype, large = "return", verbose = FALSE),
  times = times_prepare,
  unit = "ms"
)

bench_fit <- microbenchmark(
  fit_large = mirt(dat_fit, 1, itemtype = itemtype, verbose = FALSE,
                   technical = list(NCYCLES = 40L)),
  times = times_fit,
  unit = "ms"
)

mod_for_scores <- mirt(dat_fit, 1, itemtype = itemtype, verbose = FALSE,
                       technical = list(NCYCLES = 40L))
Theta_itemtrace <- matrix(seq(-6, 6, length.out = theta_n), ncol = 1L)

bench_itemtrace <- microbenchmark(
  itemtrace_large = mirt:::computeItemtrace(
    mod_for_scores@ParObjects$pars,
    Theta_itemtrace,
    mod_for_scores@Model$itemloc,
    CUSTOM.IND = mod_for_scores@Internals$CUSTOM.IND,
    omp_threads = omp_threads
  ),
  times = times_itemtrace,
  unit = "ms"
)

bench_eapsum <- microbenchmark(
  eapsum_large = fscores(mod_for_scores, method = "EAPsum", full.scores = TRUE),
  times = times_eapsum,
  unit = "ms"
)

to_row <- function(name, bench, n_obs, n_items) {
  data.frame(
    mode = mode,
    dataset = dataset,
    benchmark = name,
    N = n_obs,
    J = n_items,
    theta_n = theta_n,
    omp_threads = omp_threads,
    times = nrow(bench),
    median_ms = unname(median(bench$time) / 1e6),
    mean_ms = unname(mean(bench$time) / 1e6),
    stringsAsFactors = FALSE
  )
}

res <- rbind(
  to_row("prepare_large", bench_prepare, N_prepare, J),
  to_row("fit_large", bench_fit, N_fit, J),
  to_row("itemtrace_large", bench_itemtrace, theta_n, J),
  to_row("eapsum_large", bench_eapsum, N_fit, J)
)

if (nzchar(out)) {
  write.csv(res, file = out, row.names = FALSE)
} else {
  print(res)
}
