## REVIEW claim 1: no in-package refit recounts thresholds.
## Construct bootstrap replicates whose response loses the top category
## (and, grouped, loses it in one level) and compare the replicate's
## parameter vector, names and log-likelihood against a thres(K)-pinned
## fit on the replicate's own data. Four families, both re_formula arms,
## and the hand-written simulate() + frm() path beside the internal one.
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, " ver =",
    as.character(utils::packageVersion("frmtmb")), "\n\n")
options(width = 130)
rel <- function(a, b) abs(a - b) / max(abs(b), .Machine$double.eps)

## ---- ungrouped, no random effect -----------------------------------
set.seed(202)
x <- stats::rnorm(40)
cp <- cbind(stats::plogis(-0.6 - 0.5 * x), stats::plogis(0.5 - 0.5 * x),
            stats::plogis(2.6 - 0.5 * x))
dd <- data.frame(x = x, y = 1L + rowSums(stats::runif(40) > cp))
cat("data table:", paste(table(dd$y), collapse = "/"), "\n\n")

fams <- list(cumulative = cumulative(), sratio = sratio(),
             cratio = cratio(), acat = acat())
for (fn in names(fams)) {
  fam <- fams[[fn]]
  fit <- suppressWarnings(frm(bf(y ~ x), family = fam, data = dd))
  set.seed(11)
  sims <- simulate(fit, nsim = 60, re_formula = NA)
  lost <- which(vapply(sims, function(v) max(as.integer(v)), 1L) < 4L)
  cat("==", fn, ": fitted tau_raw =", length(fit$estimates[["tau_raw"]]),
      "; replicates losing category 4:", length(lost),
      if (length(lost)) paste0("(first = ", lost[1L], ")") else "", "\n")
  if (!length(lost)) { cat("   NO LOST REPLICATE at this seed\n"); next }
  yb <- as.integer(sims[[lost[1L]]])
  cat("   replicate table:", paste(table(yb), collapse = "/"), "\n")

  ## the in-package refit path
  rf <- suppressWarnings(refit(fit, yb))
  ## the pinned reference on the replicate's own data
  pin <- suppressWarnings(frm(bf(y | thres(3) ~ x), family = fam,
                             data = data.frame(x = dd$x, y = yb)))
  ## the path a USER writes: simulate() then frm() by hand
  hand <- suppressWarnings(frm(bf(y ~ x), family = fam,
                              data = data.frame(x = dd$x, y = yb)))
  cr <- frmtmb:::get_coef.frmtmb_fit(rf)
  cp2 <- frmtmb:::get_coef.frmtmb_fit(pin)
  ch <- frmtmb:::get_coef.frmtmb_fit(hand)
  cat("   refit  n_tau =", length(rf$estimates[["tau_raw"]]),
      " names =", paste(names(cr), collapse = ","), "\n")
  cat("   pinned n_tau =", length(pin$estimates[["tau_raw"]]),
      " names =", paste(names(cp2), collapse = ","), "\n")
  cat("   hand   n_tau =", length(hand$estimates[["tau_raw"]]),
      " names =", paste(names(ch), collapse = ","), "\n")
  cat("   names identical refit vs pinned:",
      identical(names(cr), names(cp2)), "\n")
  l1 <- as.numeric(logLik(rf)); l2 <- as.numeric(logLik(pin))
  l3 <- as.numeric(logLik(hand))
  cat(sprintf("   logLik refit  = %.13f\n", l1))
  cat(sprintf("   logLik pinned = %.13f  rel = %.4g\n", l2, rel(l1, l2)))
  cat(sprintf("   logLik hand   = %.13f  rel to pinned = %.4g\n",
              l3, rel(l3, l2)))
  ident <- setdiff(names(cr), "tau_raw_3")
  cat("   max rel diff over identified coefs:",
      signif(max(rel(cr[ident], cp2[ident])), 4), "\n")
  cat(sprintf("   tau_raw_3 refit = %.8g  pinned = %.8g\n",
              cr[["tau_raw_3"]], cp2[["tau_raw_3"]]))

  ## the bootstrap matrix itself, default re_formula and NULL
  FUNt <- function(f) c(fixef(f, flatten = TRUE),
                        tau = f$estimates[["tau_raw"]])
  for (rfm in list(NA, NULL)) {
    bs <- suppressWarnings(frm_bootstrap(fit, FUN = FUNt, nsim = 60,
                                         seed = 11, re_formula = rfm))
    cat("   frm_bootstrap(re_formula =",
        if (is.null(rfm)) "NULL" else "NA", "): dim =",
        paste(dim(bs$t), collapse = "x"), " cols =",
        paste(colnames(bs$t), collapse = ","),
        " anyNA =", anyNA(bs$t), " NA cells =", sum(is.na(bs$t)), "\n")
  }
}

## ---- grouped thresholds --------------------------------------------
cat("\n#### grouped thres(gr = g)\n")
set.seed(503)
n <- 96
g <- factor(rep(c("a", "b", "c"), length.out = n))
xg <- stats::rnorm(n)
tau <- list(a = c(-0.6, 0.6, 2.4), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
yg <- vapply(seq_len(n), function(i) {
  1L + sum(stats::runif(1) >
             stats::plogis(tau[[as.character(g[i])]] - 0.5 * xg[i]))
}, 1L)
dg <- data.frame(x = xg, g = g, y = yg)
for (fn in names(fams)) {
  fam <- fams[[fn]]
  fit <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x), family = fam,
                             data = dg))
  th <- fit$spec$responses$y$family[["thres"]]
  nraw <- length(fit$estimates[["tau_raw"]])
  set.seed(13)
  sims <- simulate(fit, nsim = 40, re_formula = NA)
  lost <- which(vapply(sims, function(v) {
    any(tapply(as.integer(v), dg$g, max) - 1L < th[["nthres"]])
  }, TRUE))
  cat("==", fn, ": nthres =", paste(th[["nthres"]], collapse = "/"),
      " nraw =", nraw, " replicates losing a level's top:", length(lost),
      "\n")
  if (!length(lost)) next
  yb <- as.integer(sims[[lost[1L]]])
  cat("   per-level max in replicate:",
      paste(tapply(yb, dg$g, max), collapse = "/"), "\n")
  rf <- suppressWarnings(refit(fit, yb))
  db <- data.frame(x = dg$x, g = dg$g, y = yb,
                   k = as.integer(th[["nthres"]][match(as.character(dg$g),
                                                       th[["levels"]])]))
  pin <- suppressWarnings(frm(bf(y | thres(k, gr = g) ~ x), family = fam,
                             data = db))
  hand <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x), family = fam,
                              data = db))
  cr <- frmtmb:::get_coef.frmtmb_fit(rf)
  cp2 <- frmtmb:::get_coef.frmtmb_fit(pin)
  ch <- frmtmb:::get_coef.frmtmb_fit(hand)
  cat("   refit n_tau =", length(rf$estimates[["tau_raw"]]),
      " pinned =", length(pin$estimates[["tau_raw"]]),
      " hand =", length(hand$estimates[["tau_raw"]]), "\n")
  cat("   names identical refit vs pinned:",
      identical(names(cr), names(cp2)),
      " refit vs hand:", identical(names(cr), names(ch)), "\n")
  l1 <- as.numeric(logLik(rf)); l2 <- as.numeric(logLik(pin))
  cat(sprintf("   logLik refit = %.13f  pinned = %.13f  rel = %.4g\n",
              l1, l2, rel(l1, l2)))
  cat(sprintf("   logLik hand  = %.13f  rel to pinned = %.4g\n",
              as.numeric(logLik(hand)), rel(as.numeric(logLik(hand)), l2)))
  FUNt <- function(f) c(fixef(f, flatten = TRUE),
                        tau = f$estimates[["tau_raw"]])
  bs <- suppressWarnings(frm_bootstrap(fit, FUN = FUNt, nsim = 40, seed = 13))
  cat("   frm_bootstrap dim =", paste(dim(bs$t), collapse = "x"),
      " NA cells =", sum(is.na(bs$t)), "\n")
}
cat("\nDONE rev-01\n")
