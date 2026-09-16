## Reviewer: verify the two diagnostics against brms's OWN function
## body, on draws this script builds itself, and re-derive the cached
## draws object bitwise.
##
##   Rscript dev/bmrev-diag.R
##
## Two constructions:
##  A. INDEPENDENT draws: data seed 17, frm_sample(seed = 31337).
##  B. the lane's construction repeated: data seed 9, seed 20260915,
##     compared to dev/stan-cache/brmsmatch-draws.rds with identical().
LIB <- "C:/Users/adf44/source/r/bmrev-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
q(library(posterior)); q(library(brms))
cat("frmtmb        ", format(packageVersion("frmtmb")), "\n")
cat("frmtmb.sample ", format(packageVersion("frmtmb.sample")),
    " from ", dirname(system.file(package = "frmtmb.sample")), "\n")
cat("StanHeaders   ", format(packageVersion("StanHeaders")), "\n")
cat("brms          ", format(packageVersion("brms")), "\n")
cat("posterior     ", format(packageVersion("posterior")), "\n\n")

mkdata <- function(seed, n = 120) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:6, n / 6)))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x +
                         stats::rnorm(6, 0, 0.5)[dd$g], 1)
  dd
}
fitds <- function(dseed, sseed) {
  fit <- frm(frmtmb::bf(y ~ x + (1 | g)), family = gaussian(),
             data = mkdata(dseed))
  q(frm_sample(fit, chains = 4, iter = 1000, refresh = 0, seed = sseed))
}

## ---- B first: is the cached object reproducible bitwise? ------------
cat("== B. the lane's cache, rebuilt ==\n")
cache <- "dev/stan-cache/brmsmatch-draws.rds"
t0 <- proc.time()[["elapsed"]]
dsB <- fitds(9, 20260915)
cat("rebuild took ", round(proc.time()[["elapsed"]] - t0, 1), " s\n")
if (file.exists(cache)) {
  old <- readRDS(cache)
  cat("identical(draws matrix):        ",
      identical(old$draws, dsB$draws), "\n")
  cat("max |old - new| over the matrix:",
      format(max(abs(as.matrix(old$draws) - as.matrix(dsB$draws))),
             digits = 17), "\n")
  cat("dim old / new:                  ",
      paste(dim(old$draws), collapse = "x"), " / ",
      paste(dim(dsB$draws), collapse = "x"), "\n")
} else cat("cache ABSENT\n")

## ---- A: independent draws -------------------------------------------
cat("\n== A. independent draws, data seed 17, sampler seed 31337 ==\n")
ds <- fitds(17, 31337)
cat("draws ", nrow(ds$draws), " x ", ncol(ds$draws), "\n")
cat("variables(ds): ", paste(variables(ds), collapse = " "), "\n")

## ---- run brms's OWN body on OUR array, through a shim ---------------
## brms's rhat.brmsfit / neff_ratio.brmsfit touch `x` only through
## contains_draws(x) and as_draws_array(x, variable = pars, ...).
arr <- posterior::as_draws_array(ds)
shim <- structure(list(), class = "bmrevshim")
as_draws_array.bmrevshim <- function(x, variable = NULL, ...) {
  if (is.null(variable)) arr else
    posterior::subset_draws(arr, variable = variable, ...)
}
registerS3method("as_draws_array", "bmrevshim", as_draws_array.bmrevshim,
                 envir = asNamespace("posterior"))
utils::assignInNamespace("contains_draws",
                         function(x, ...) invisible(TRUE), ns = "brms")
brms_rhat <- brms:::rhat.brmsfit(shim)
brms_neff <- brms:::neff_ratio.brmsfit(shim)

ours_rhat <- rhat(ds)
ours_neff <- neff_ratio(ds)

cmp <- function(lbl, a, b) {
  cat("-- ", lbl, " --\n", sep = "")
  cat("  names equal:      ", identical(names(a), names(b)), "\n")
  cat("  identical():      ", identical(a, b), "\n")
  cat("  identical(unname):", identical(unname(a), unname(b)), "\n")
  k <- intersect(names(a), names(b))
  cat("  max |diff|:       ",
      format(max(abs(a[k] - b[k])), digits = 17), "\n")
  cat("  brms's names:     ", paste(names(b), collapse = " "), "\n")
  cat("  ours' names:      ", paste(names(a), collapse = " "), "\n")
}
cat("\n== brms's OWN rhat.brmsfit body, run on our draws array ==\n")
cmp("rhat", ours_rhat, brms_rhat)
cat("\n== brms's OWN neff_ratio.brmsfit body, run on our array ==\n")
cmp("neff_ratio", ours_neff, brms_neff)

## a control: the two must NOT agree with rstan's, or the test is empty
sm <- rstan::summary(ds$stanfit)$summary
cat("\n-- control: rstan's Rhat on the same fit --\n")
cat("  rstan Rhat names: ", paste(head(rownames(sm), 12),
                                  collapse = " "), "\n")
cat("  max |ours rhat - 1|:  ",
    format(max(abs(ours_rhat - 1)), digits = 8), "\n")
cat("  max |rstan Rhat - 1|: ",
    format(max(abs(sm[, "Rhat"] - 1)), digits = 8), "\n")
cat("  sorted values differ: ",
    !isTRUE(all.equal(sort(unname(ours_rhat)),
                      sort(unname(sm[, "Rhat"])))), "\n")
cat("  max rel diff, sorted: ",
    format(max(abs(sort(unname(ours_rhat)) - sort(unname(sm[, "Rhat"]))) /
                 sort(unname(sm[, "Rhat"]))), digits = 8), "\n")
nrs <- sm[, "n_eff"] / nrow(ds$draws)
cat("  max rel diff neff, sorted: ",
    format(max(abs(sort(unname(ours_neff)) - sort(unname(nrs))) /
                 sort(unname(nrs))), digits = 8), "\n")

## ---- 2. summary(ds) vs rhat(ds): does the package disagree? ---------
cat("\n== 2. summary(ds) against rhat(ds) ==\n")
s <- summary(ds)
cat("class(summary(ds)): ", paste(class(s), collapse = " "), "\n")
str(s, max.level = 1)
cat("DONE-A\n")
saveRDS(ds, "dev/stan-cache/bmrev-draws.rds")
cat("saved dev/stan-cache/bmrev-draws.rds\n")
