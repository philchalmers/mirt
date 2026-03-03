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

N_prepare <- arg_int("n_prepare", 50000L)
N_fit <- arg_int("n_fit", 12000L)
J <- arg_int("j", 40L)
times_prepare <- arg_int("times_prepare", 5L)
times_fit <- arg_int("times_fit", 3L)
times_eapsum <- arg_int("times_eapsum", 5L)

if (!mode %in% c("baseline", "head")) {
  stop("--mode must be one of: baseline, head")
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
a <- matrix(abs(rnorm(J, mean = 1.1, sd = 0.2)), ncol = 1L)
d <- rnorm(J, mean = -0.1, sd = 0.9)

dat_prepare <- simdata(a, d, N = N_prepare, itemtype = "2PL")
dat_fit <- simdata(a, d, N = N_fit, itemtype = "2PL")

bench_prepare <- microbenchmark(
  prepare_large = mirt(dat_prepare, 1, itemtype = "2PL", large = "return", verbose = FALSE),
  times = times_prepare,
  unit = "ms"
)

bench_fit <- microbenchmark(
  fit_large = mirt(dat_fit, 1, itemtype = "2PL", verbose = FALSE,
                   technical = list(NCYCLES = 40L)),
  times = times_fit,
  unit = "ms"
)

mod_for_scores <- mirt(dat_fit, 1, itemtype = "2PL", verbose = FALSE,
                       technical = list(NCYCLES = 40L))

bench_eapsum <- microbenchmark(
  eapsum_large = fscores(mod_for_scores, method = "EAPsum", full.scores = TRUE),
  times = times_eapsum,
  unit = "ms"
)

to_row <- function(name, bench, n_obs, n_items) {
  data.frame(
    mode = mode,
    benchmark = name,
    N = n_obs,
    J = n_items,
    times = nrow(bench),
    median_ms = unname(median(bench$time) / 1e6),
    mean_ms = unname(mean(bench$time) / 1e6),
    stringsAsFactors = FALSE
  )
}

res <- rbind(
  to_row("prepare_large", bench_prepare, N_prepare, J),
  to_row("fit_large", bench_fit, N_fit, J),
  to_row("eapsum_large", bench_eapsum, N_fit, J)
)

if (nzchar(out)) {
  write.csv(res, file = out, row.names = FALSE)
} else {
  print(res)
}
