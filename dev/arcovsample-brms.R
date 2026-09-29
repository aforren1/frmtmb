# Lane wt-arcovsample, validation (b): frmtmb.sample's log_lik() row by
# row against brms 2.23.0's, at the SAME parameter values.
#
# Construction. brms is fitted briefly (it only has to exist, and its
# posterior is irrelevant), then each of its draws is READ OUT and
# written into a frmtmb draws object by parameter name. Both packages
# then report a row log-density at one parameter vector, and the two
# matrices are compared cell by cell. Nothing here relies on the two
# samplers agreeing.
#
# Two orderings have to line up and are checked rather than assumed:
# brms sorts its data by order(gr, time) and undoes that with
# reorder_obs(old_order) before returning log_lik, so both matrices are
# in the USER's row order.
#
#   Rscript dev/arcovsample-brms.R > dev/arcovsample-log/brms.txt

LIB <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
  library(brms)
})
cat("frmtmb        ", format(packageVersion("frmtmb")), "\n")
cat("frmtmb.sample ", format(packageVersion("frmtmb.sample")), "\n")
cat("brms          ", format(packageVersion("brms")), "\n")
cat("rstan         ", format(packageVersion("rstan")), "\n")

SEED <- 4022L
set.seed(SEED)
ng <- 8L
nt <- 9L
dd <- expand.grid(t = seq_len(nt), g = factor(seq_len(ng)))
dd <- dd[-c(4L, 20L, 33L), ]           # interior gaps: row-counted lags
dd$x <- stats::rnorm(nrow(dd))
u <- stats::rnorm(ng, 0, 0.6)
dd$y <- 1 + 0.5 * dd$x + u[as.integer(dd$g)] + stats::rnorm(nrow(dd))
dd <- dd[sample(nrow(dd)), ]           # shuffled: the sort is exercised
rownames(dd) <- NULL
cat("N =", nrow(dd), " groups =", ng, " seed =", SEED, "\n")

# frmtmb's stored draws matrix carries brms's spelling and brms's scale
# where brms has the same parameter (`b_*`, `sigma`, `nu`, `r_g[...]`),
# and frmtmb's internal one where it does not: `theta_1` is a LOG
# standard deviation and `thetaac_k` the raw ARMA coefficients, AR
# before MA. `src` names brms's column for each of the second kind and
# `fun` the map onto it. Every frmtmb column has to be accounted for or
# the run stops: a silently unfilled column would compare one model
# against a different one, at zero.
transplant <- function(ds, bfit) {
  bm <- as.matrix(posterior::as_draws_matrix(bfit))
  want <- setdiff(colnames(ds$draws), "lp__")
  ac <- ds$fit$frame[["autocor"]][[1L]]
  actgt <- c(if (ac[["p"]]) paste0("ar[", seq_len(ac[["p"]]), "]"),
             if (ac[["q"]]) paste0("ma[", seq_len(ac[["q"]]), "]"))
  src <- stats::setNames(
    actgt, paste0("thetaac_", seq_along(actgt)))
  if ("theta_1" %in% want) src <- c(src, theta_1 = "sd_g__Intercept")
  direct <- setdiff(want, names(src))
  miss <- c(setdiff(direct, colnames(bm)),
            setdiff(unname(src), colnames(bm)))
  cat("  frmtmb columns : ", paste(want, collapse = ", "), "\n", sep = "")
  cat("  mapped         : ",
      paste(names(src), "<-", src, collapse = ", "), "\n", sep = "")
  cat("  unmatched      : ",
      if (length(miss)) paste(miss, collapse = ", ") else "(none)", "\n",
      sep = "")
  stopifnot(length(miss) == 0L)
  nd <- nrow(bm)
  out <- matrix(NA_real_, nd, ncol(ds$draws),
                dimnames = list(NULL, colnames(ds$draws)))
  out[, direct] <- bm[, direct, drop = FALSE]
  for (nm in names(src)) {
    v <- bm[, src[[nm]]]
    out[, nm] <- if (identical(nm, "theta_1")) log(v) else v
  }
  out[, "lp__"] <- 0
  stopifnot(!anyNA(out))
  ds$draws <- out
  ds
}

one <- function(label, form, bform, family, bfamily) {
  cat("\n================ ", label, " ================\n", sep = "")
  # a short brms run: only its draws and its log_lik() are wanted
  bfit <- suppressWarnings(suppressMessages(
    brms::brm(bform, data = dd, family = bfamily, chains = 1,
              iter = 300, warmup = 200, refresh = 0, seed = 5,
              silent = 2)))
  # frm_sample() supplies the column spelling and the draws object shape
  ds <- suppressWarnings(suppressMessages(
    frm_sample(form, family = family, data = dd, chains = 1,
               iter = 300, refresh = 0, seed = 5)))
  dsb <- transplant(ds, bfit)

  # The transplant is checked BEFORE the claim rests on it: the one-step
  # mean is a function of every transplanted parameter, so agreeing
  # epreds mean the vectors are the same vector. If this line is large,
  # the log_lik comparison below says nothing.
  epb <- brms::posterior_epred(bfit)
  epf <- posterior_epred(dsb)
  cat("  max |posterior_epred diff| : ",
      format(max(abs(epf - epb)), digits = 12), " (epred sd ",
      format(stats::sd(epb), digits = 5), ")\n", sep = "")

  llb <- brms::log_lik(bfit)
  llf <- log_lik(dsb)
  cat("  brms log_lik dim : ", paste(dim(llb), collapse = " x "), "\n",
      sep = "")
  cat("  frmtmb log_lik dim : ", paste(dim(llf), collapse = " x "), "\n",
      sep = "")
  stopifnot(identical(dim(llb), dim(llf)))
  ad <- abs(llf - llb)
  rl <- ad / pmax(abs(llb), .Machine$double.eps)
  cat("  max |diff| per cell        : ", format(max(ad), digits = 12),
      "\n", sep = "")
  cat("  max relative diff per cell : ", format(max(rl), digits = 12),
      "\n", sep = "")
  cat("  cell with the worst diff   : draw ",
      which(ad == max(ad), arr.ind = TRUE)[1L, 1L], " column ",
      which(ad == max(ad), arr.ind = TRUE)[1L, 2L], "\n", sep = "")
  cat("  spread of the cells compared (sd of brms log_lik) : ",
      format(stats::sd(llb), digits = 6), "\n", sep = "")
  # the FIRST rows of each group are the ones with no full lag; report
  # their agreement separately, since that is where brms's convention
  # (e_s = 0 before a group starts) could differ silently
  first <- unlist(lapply(split(seq_len(nrow(dd)), dd$g),
                         function(r) r[order(dd$t[r])][1L]))
  cat("  max |diff| on each group's FIRST row : ",
      format(max(ad[, first]), digits = 12), "\n", sep = "")

  # (c) loose: each package's loo() on its OWN draws
  lb <- suppressWarnings(brms::loo(bfit))
  lf <- suppressWarnings(loo(ds))
  cat("  brms elpd_loo (its own draws)   : ",
      format(lb$estimates["elpd_loo", "Estimate"], digits = 8), " se ",
      format(lb$estimates["elpd_loo", "SE"], digits = 5), "\n", sep = "")
  cat("  frmtmb elpd_loo (its own draws) : ",
      format(lf$estimates["elpd_loo", "Estimate"], digits = 8), " se ",
      format(lf$estimates["elpd_loo", "SE"], digits = 5), "\n", sep = "")
  # and on the SAME draws, where the difference is the log_lik matrix
  # alone: this one IS the cell comparison above, aggregated
  lt <- suppressWarnings(loo(dsb))
  cat("  brms elpd_loo, brms draws       : ",
      format(lb$estimates["elpd_loo", "Estimate"], digits = 8), "\n",
      sep = "")
  cat("  frmtmb elpd_loo, brms draws     : ",
      format(lt$estimates["elpd_loo", "Estimate"], digits = 8), "\n",
      sep = "")
  invisible(NULL)
}

one("gaussian ar(1)", frmtmb::bf(y ~ x + ar(t, g)),
    brms::bf(y ~ x + ar(t, g)), gaussian(), brms::brmsfamily("gaussian"))

one("student arma(1,1)", frmtmb::bf(y ~ x + arma(t, g, p = 1, q = 1)),
    brms::bf(y ~ x + arma(t, g, p = 1, q = 1)), frmtmb::student(),
    brms::brmsfamily("student"))

one("gaussian arma(1,1) + (1 | g)",
    frmtmb::bf(y ~ x + arma(t, g, p = 1, q = 1) + (1 | g)),
    brms::bf(y ~ x + arma(t, g, p = 1, q = 1) + (1 | g)),
    gaussian(), brms::brmsfamily("gaussian"))

cat("\nDONE\n")
