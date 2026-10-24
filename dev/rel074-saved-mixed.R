# A saved fit whose FAMILY is from one extension (cox(), frmtmb.survival)
# and whose STRUCTURE is from another (car(), frmtmb.spatial), read back
# in sessions that lack one, the other or both: the cross-lane check of
# the two absent-package refusals (require_model_packages() and
# require_structure_packages(), R/extension-home.R). Every method must
# refuse with class frmtmb_extension_missing and name a missing package,
# never "could not find function"; with both present every method works.
#   Rscript dev/rel074-saved-mixed.R save <rds>       (rellib-r12)
#   Rscript dev/rel074-saved-mixed.R read <rds> <lib;lib>
a <- commandArgs(trailingOnly = TRUE)
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
if (a[1] == "save") {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r12", user))
  suppressMessages({
    library(frmtmb); library(frmtmb.survival); library(frmtmb.spatial)
  })
  g <- expand.grid(r = 1:4, c = 1:4)
  n <- nrow(g)
  W <- matrix(0, n, n)
  for (i in 1:n) for (j in 1:n) {
    if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) W[i, j] <- 1
  }
  dimnames(W) <- list(paste0("L", 1:n), paste0("L", 1:n))
  set.seed(20261023)
  dd <- data.frame(loc = factor(rep(paste0("L", 1:n), each = 12),
                                levels = paste0("L", 1:n)),
                   x = rnorm(12 * n))
  u <- rnorm(n, 0, 0.5)
  dd$time <- rexp(nrow(dd), exp(-0.5 + 0.6 * dd$x + u[as.integer(dd$loc)]))
  dd$out <- rbinom(nrow(dd), 1, 0.25)
  fit <- frm(bf(time | cens(out) ~ x + car(W, gr = loc, type = "escar")),
             family = cox(df = 5), data = dd, data2 = list(W = W))
  cat("fit: logLik", format(as.numeric(logLik(fit)), digits = 12),
      "conv", fit$opt$convergence, "\n")
  saveRDS(list(fit = fit, dd = dd, W = W), a[2])
  quit(save = "no")
}
.libPaths(c(strsplit(a[3], ";", fixed = TRUE)[[1]], user))
suppressMessages(library(frmtmb))
inst <- function(p) nzchar(system.file(package = p))
cat("lib:", dirname(find.package("frmtmb")), "| survival installed:",
    inst("frmtmb.survival"), "| spatial installed:", inst("frmtmb.spatial"),
    "\n")
s <- readRDS(a[2])
fit <- s$fit
nd <- s$dd[1:3, ]
calls <- list(
  print = quote(utils::capture.output(print(fit))),
  summary = quote(summary(fit)),
  logLik = quote(logLik(fit)),
  vcov = quote(vcov(fit)),
  VarCorr = quote(VarCorr(fit)),
  ranef = quote(ranef(fit)),
  fitted = quote(fitted(fit)),
  predict_new = quote(frm_linpred(fit, newdata = nd, se.fit = TRUE)),
  simulate = quote(simulate(fit, nsim = 1, seed = 1)),
  confint_varcorr = quote(confint_varcorr(fit)),
  update = quote(update(fit, . ~ . - x)),
  frm_curve = quote(frm_curve(fit, newdata = data.frame(x = c(-1, 0, 1)))),
  frm_curve_re = quote(frm_curve(fit, newdata = nd, re_formula = NULL)))
for (nm in names(calls)) {
  r <- tryCatch({ eval(calls[[nm]]); "ok" }, error = function(e) {
    k <- if (inherits(e, "frmtmb_extension_missing")) "REFUSED[ext]" else
      "ERR"
    paste(k, substr(gsub("\\s+", " ", conditionMessage(e)), 1, 120))
  })
  cat(sprintf("  %-16s %s\n", nm, r))
}
