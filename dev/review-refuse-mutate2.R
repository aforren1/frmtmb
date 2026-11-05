# Item 9, second set of mutants: one per process.
id <- commandArgs(trailingOnly = TRUE)[1L]
source("dev/review-refuse-lib.R")
library(frmtmb.sample)
ns <- asNamespace("frmtmb")
mutate <- function(name, fn, nsp = "frmtmb") assignInNamespace(name, fn, nsp)
tf <- switch(id,
  sds_reach_re = "tests/testthat/test-prior-sds.R",
  key_off = "tests/testthat/test-mm-id-cov.R",
  cor_key = "tests/testthat/test-refuse-pairs.R",
  mis_inv = "tests/testthat/test-mi-families.R",
  mis_inv2 = "tests/testthat/test-refuse-pairs.R",
  place0 = "tests/testthat/test-mi-families.R",
  mm_cov_arg = "tests/testthat/test-mm-id-cov.R",
  stop("unknown"))
switch(id,
  sds_reach_re = {
    orig <- ns$sd_spec_reach
    mutate("sd_spec_reach", function(spec, frame, bk, s) {
      if (!identical(bk[["covstruct"]], "smooth") && identical(s$class, "sds")) {
        s$class <- "sd"
      }
      orig(spec, frame, bk, s)
    })
  },
  key_off = mutate("bar_read_id", function(bar, env) {
    rhs <- frmtmb:::strip_parens(bar[[3L]])
    if (is.call(rhs) && identical(rhs[[1L]], as.name("mm"))) {
      rhs[["id"]] <- NULL
      bar[[3L]] <- rhs
    }
    bar
  }),
  cor_key = {
    orig <- ns$bar_read_cor
    mutate("bar_read_cor", function(bar, env) {
      lhs <- frmtmb:::strip_parens(bar[[2L]])
      if (is.call(lhs) && identical(lhs[[1L]], as.name("|")) &&
            grepl("cor = FALSE", deparse1(bar[[3L]]), fixed = TRUE)) {
        stop("an |ID| key is refused")
      }
      orig(bar, env)
    })
  },
  mis_inv = , mis_inv2 = {
    # the reader alone keeps the chart: mi_values() returns u
    mutate("mi_values", function(fit, vn) {
      xv <- fit$frame[["y"]][[vn]]
      mm_ <- fit$frame[["mi_map"]][[vn]]
      if (!is.null(mm_)) xv[mm_$rows] <- fit$estimates[["miss"]][mm_$idx]
      xv
    })
  },
  place0 = mutate("mi_placeholder", function(obs, chart) 0),
  mm_cov_arg = {
    orig <- ns$gr_cov_check_within
    mutate("gr_cov_check_within", function(re, cnms) invisible(NULL))
  })
pkg <- if (grepl("frmtmb.sample", tf)) "frmtmb.sample" else "frmtmb"
env <- testthat::test_env(pkg)
library(pkg, character.only = TRUE)
res <- as.data.frame(testthat::test_file(tf, package = pkg, env = env, reporter = "silent", stop_on_failure = FALSE))
pr("MUT %-20s %-40s pass=%d fail=%d error=%d skip=%d; failing: %s", id, basename(tf), sum(res$passed), sum(res$failed), sum(res$error), sum(res$skipped),
   paste(substr(res$test[res$failed > 0 | res$error], 1, 60), collapse = " | "))
