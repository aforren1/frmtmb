## Item 2.5f measurement harness. Runs against whatever frmtmb.sample
## is first on .libPaths(), so the SAME script gives the before and the
## after numbers.
##
## Usage:
##   Rscript dev/brmsmatch-measure.R before
##   Rscript dev/brmsmatch-measure.R after
##
## The draws it measures are cached at dev/stan-cache/brmsmatch-draws.rds:
## 4 chains x 1000 iterations, sampler seed 20260915, data seed 9, the
## same construction dev/sgrev-rhat.R used (that file's cache was gone).
arg <- commandArgs(trailingOnly = TRUE)
tag <- if (length(arg)) arg[[1L]] else "run"
lib <- if (identical(tag, "before")) {
  ## the base build, read-only
  c("C:/Users/adf44/source/r/rellib-r3",
    "C:/Users/adf44/source/r/pinlib",
    "C:/Users/adf44/AppData/Local/R/win-library/4.6")
} else {
  c("C:/Users/adf44/source/r/brmsmatch-lib",
    "C:/Users/adf44/source/r/rellib-r3",
    "C:/Users/adf44/source/r/pinlib",
    "C:/Users/adf44/AppData/Local/R/win-library/4.6")
}
.libPaths(lib)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
q(library(posterior))

cat("== build ==\n")
cat("tag                ", tag, "\n")
cat("frmtmb             ", format(packageVersion("frmtmb")), "\n")
cat("frmtmb.sample      ", format(packageVersion("frmtmb.sample")),
    " from ", dirname(system.file(package = "frmtmb.sample")), "\n")
cat("StanHeaders        ", format(packageVersion("StanHeaders")), "\n")
cat("brms               ", format(packageVersion("brms")), "\n")
cat("posterior          ", format(packageVersion("posterior")), "\n")

cache <- "dev/stan-cache/brmsmatch-draws.rds"
if (file.exists(cache)) {
  ds <- readRDS(cache)
} else {
  set.seed(9)
  dd <- data.frame(x = stats::rnorm(120), g = factor(rep(1:6, 20)))
  dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x +
                         stats::rnorm(6, 0, 0.5)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  ds <- q(frm_sample(fit, chains = 4, iter = 1000, refresh = 0,
                     seed = 20260915))
  saveRDS(ds, cache)
}
cat("draws              ", nrow(ds$draws), " x ", ncol(ds$draws), "\n\n")

arr <- posterior::as_draws_array(ds)
nd <- posterior::ndraws(arr)

## ---- 1. the two diagnostics against brms's definition ---------------
## brms::rhat.brmsfit is summarise_draws(rhat = posterior::rhat) and
## brms::neff_ratio.brmsfit is min(ess_bulk, ess_tail) / ndraws, both
## read off the installed brms in dev/brmsmatch-formals.R.
sdr <- posterior::summarise_draws(arr, rhat = posterior::rhat)
brms_rhat <- stats::setNames(sdr$rhat, sdr$variable)
sde <- posterior::summarise_draws(arr, ess_bulk = posterior::ess_bulk,
                                  ess_tail = posterior::ess_tail)
brms_neff <- stats::setNames(pmin(sde$ess_bulk, sde$ess_tail) / nd,
                             sde$variable)

report <- function(label, got, want) {
  cat("-- ", label, " --\n", sep = "")
  cm <- intersect(names(got), names(want))
  cat("names in got but not in variables(ds): ",
      paste(setdiff(names(got), posterior::variables(ds)),
            collapse = " "), "\n")
  cat("names in variables(ds) but not in got: ",
      paste(setdiff(posterior::variables(ds), names(got)),
            collapse = " "), "\n")
  cat("common names: ", length(cm), " of ", length(want), "\n")
  if (!length(cm)) return(invisible(NULL))
  rel <- abs(got[cm] - want[cm]) / abs(want[cm])
  cat("max relative difference vs brms: ",
      format(max(rel), digits = 6), "\n")
  cat("max absolute difference vs brms: ",
      format(max(abs(got[cm] - want[cm])), digits = 6), "\n")
  cat("identical to brms's:             ",
      identical(unname(got[cm]), unname(want[cm])), "\n")
  invisible(NULL)
}

cat("== 1. rhat() and neff_ratio() against brms's definition ==\n")
r <- rhat(ds)
report("rhat(ds)", r, brms_rhat)
cat("max |brms rhat - 1|, the whole signal: ",
    format(max(abs(brms_rhat - 1)), digits = 6), "\n")
cat("difference as a fraction of (rhat - 1): ",
    format(max(abs(r[intersect(names(r), names(brms_rhat))] -
                     brms_rhat[intersect(names(r), names(brms_rhat))]) /
                 abs(brms_rhat[intersect(names(r),
                                         names(brms_rhat))] - 1)),
           digits = 6), "\n\n")
n <- neff_ratio(ds)
report("neff_ratio(ds)", n, brms_neff)

cat("\n== 2. the name check ==\n")
cat("variables(ds)        : ", paste(posterior::variables(ds),
                                     collapse = " "), "\n")
cat("names(rhat(ds))      : ", paste(names(r), collapse = " "), "\n")
cat("names(neff_ratio(ds)): ", paste(names(n), collapse = " "), "\n")
cat("rhat(ds)[\"x\"]        : ", format(unname(r["x"])), "\n")
cat("neff_ratio(ds)[\"x\"]  : ", format(unname(n["x"])), "\n")
cat("rhat(ds)[\"x\"] is NA  : ", is.na(r["x"]), "\n")

## the pars= argument both brms methods carry
cat("\nrhat(ds, \"^s\") : ")
print(try(rhat(ds, "^s"), silent = TRUE))
cat("neff_ratio(ds, \"^s\") : ")
print(try(neff_ratio(ds, "^s"), silent = TRUE))

## ---- 3. the ten positional signatures -------------------------------
cat("\n== 3. the ten positional signatures ==\n")
nd2 <- data.frame(x = c(-1, 0, 1),
                  g = factor(1, levels = levels(ds$fit$frame[["data"]]$g)))
if (is.null(levels(ds$fit$frame[["data"]]$g))) {
  nd2 <- data.frame(x = c(-1, 0, 1), g = factor(1, levels = as.character(1:6)))
}

shape <- function(v) {
  if (inherits(v, "try-error")) {
    return(paste0("ERROR: ",
                  sub("\n.*$", "", sub("^Error[^:]*: ", "",
                                       as.character(v)))))
  }
  d <- dim(v)
  paste0(class(v)[1L], " ",
         if (is.null(d)) paste0("len=", length(v)) else
           paste(d, collapse = "x"))
}

pos <- list(
  `as.mcmc(ds, TRUE)`               = quote(as.mcmc(ds, TRUE)),
  `as.mcmc(ds, "^s")`               = quote(as.mcmc(ds, "^s")),
  `log_lik(ds, nd2)`                = quote(log_lik(ds, nd2)),
  `mcmc_plot(ds, "^s")`             = quote(mcmc_plot(ds, "^s")),
  `posterior_epred(ds, nd2, NA)`    = quote(posterior_epred(ds, nd2, NA)),
  `posterior_interval(ds, 0.9)`     = quote(posterior_interval(ds, 0.9)),
  `posterior_interval(ds, "^s")`    = quote(posterior_interval(ds, "^s")),
  `posterior_linpred(ds, F, nd2, NA)` =
    quote(posterior_linpred(ds, FALSE, nd2, NA)),
  `posterior_predict(ds, nd2, NA)`  = quote(posterior_predict(ds, nd2, NA)),
  `pp_mixture(ds, nd2)`             = quote(pp_mixture(ds, nd2)),
  `predictive_error(ds, nd2)`       = quote(predictive_error(ds, nd2)),
  `psis(ds, nd2)`                   = quote(psis(ds, nd2))
)
for (nm in names(pos)) {
  v <- try(suppressWarnings(suppressMessages(eval(pos[[nm]]))),
           silent = TRUE)
  cat(sprintf("%-36s %s\n", nm, shape(v)))
}

## the positional slots that now carry a value rather than a refusal
cat("\n== 4. what the newly positional slots DO ==\n")
ndy <- nd2
ndy$y <- c(0.5, 1.0, 1.5)
chk <- function(nm, e) {
  v <- try(suppressWarnings(suppressMessages(e)), silent = TRUE)
  cat(sprintf("%-46s %s\n", nm, shape(v)))
  invisible(v)
}
ep_na <- chk("posterior_epred(ds, nd2, NA)",
             posterior_epred(ds, nd2, NA))
ep_nu <- chk("posterior_epred(ds, nd2, NULL)",
             posterior_epred(ds, nd2, NULL))
if (!inherits(ep_na, "try-error") && !inherits(ep_nu, "try-error")) {
  cat("  the two differ (re_formula = NA really drops the RE): ",
      !isTRUE(all.equal(ep_na, ep_nu)), "\n")
  cat("  max |NA - NULL| column mean: ",
      format(max(abs(colMeans(ep_na) - colMeans(ep_nu))), digits = 6),
      "\n")
}
pe <- chk("predictive_error(ds, ndy)", predictive_error(ds, ndy))
if (!inherits(pe, "try-error")) {
  set.seed(1)
  cat("  equals y - posterior_predict(newdata) columnwise (means): ",
      format(max(abs(colMeans(pe) -
                       (ndy$y - colMeans(posterior_predict(ds,
                                                           newdata = ndy))))),
             digits = 6), "\n")
}
pee <- chk("predictive_error(ds, ndy, method = 'posterior_epred')",
           predictive_error(ds, ndy, method = "posterior_epred"))
if (!inherits(pee, "try-error")) {
  cat("  equals y - posterior_epred(newdata), max abs: ",
      format(max(abs(pee - sweep(-posterior_epred(ds, newdata = ndy),
                                 2L, ndy$y, "+"))), digits = 6), "\n")
}
chk("posterior_predict(ds, nd2, NULL, NULL, exp)",
    posterior_predict(ds, nd2, NULL, NULL, exp))
chk("psis(ds, NULL, NULL, 'lab')", psis(ds, NULL, NULL, "lab"))
chk("log_lik(ds, NULL, NULL, NULL, 10)",
    log_lik(ds, NULL, NULL, NULL, 10))
chk("as.mcmc(ds, NA, FALSE, TRUE)", as.mcmc(ds, NA, FALSE, TRUE))
chk("as.mcmc(ds, NA, FALSE, FALSE, TRUE)",
    as.mcmc(ds, NA, FALSE, FALSE, TRUE))
chk("mcmc_plot(ds, NA, 'trace')", mcmc_plot(ds, NA, "trace"))
chk("posterior_interval(ds, NA, 'x', 0.9)",
    posterior_interval(ds, NA, "x", 0.9))
chk("posterior_linpred(ds, TRUE, nd2, NA)",
    posterior_linpred(ds, TRUE, nd2, NA))

cat("\n== formals of the ten methods, here ==\n")
tb <- get(".__S3MethodsTable__.", envir = asNamespace("frmtmb.sample"),
          inherits = FALSE)
for (k in c("as.mcmc.frmtmb_draws", "log_lik.frmtmb_draws",
            "mcmc_plot.frmtmb_draws", "posterior_epred.frmtmb_draws",
            "posterior_interval.frmtmb_draws",
            "posterior_linpred.frmtmb_draws",
            "posterior_predict.frmtmb_draws", "pp_mixture.frmtmb_draws",
            "predictive_error.frmtmb_draws", "psis.frmtmb_draws",
            "rhat.frmtmb_draws", "neff_ratio.frmtmb_draws")) {
  f <- tryCatch(get(k, envir = tb, inherits = FALSE),
                error = function(e) NULL)
  cat(sprintf("%-34s %s\n", k,
              if (is.null(f)) "MISSING" else
                paste(names(formals(f)), collapse = " ")))
}
cat("DONE\n")
