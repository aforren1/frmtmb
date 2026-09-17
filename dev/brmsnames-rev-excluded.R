## Reviewer, claim 4: blocks without r_ names (rr, smooth) and the
## duplicate fallback (animal model), on fits and draws.
##   Rscript dev/brmsnames-rev-excluded.R      data seed 21, sampler seed 4
.libPaths(c("C:/Users/adf44/source/r/brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
set.seed(21)
G <- 12; n <- 240
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), g = factor(rep(1:G, 20)),
                id = factor(rep(1:G, 20)))
A <- diag(G); A[cbind(1:(G - 1), 2:G)] <- 0.25; A[cbind(2:G, 1:(G - 1))] <- 0.25
rownames(A) <- colnames(A) <- levels(d$id)
d$y <- 1 + d$x1 + rnorm(G)[d$g] * d$x1 + rnorm(G)[d$g] * 0.5 * d$x2 + rnorm(n)
t1 <- function(e) tryCatch(q(e), error = function(err) paste("ERROR:", conditionMessage(err)))
ms <- list(
  rr = list(bf(y ~ x1 + rr(x1 + x2 | g, d = 1)), NULL),
  smooth_re = list(bf(y ~ s(x1) + (1 + x2 | g)), NULL),
  animal = list(bf(y ~ 1 + (1 | gr(id, cov = A)) + (1 | id)), list(A = A))
)
for (nm in names(ms)) {
  cat("\n==", nm, "==\n")
  fit <- t1(frm(ms[[nm]][[1]], family = gaussian(), data = d, data2 = ms[[nm]][[2]]))
  if (is.character(fit)) { cat(fit, "\n"); next }
  cat("variables(fit):", variables(fit), "\n")
  vc <- t1(VarCorr(fit)); cat("VarCorr(fit) names:", if (is.character(vc)) vc else names(vc), "\n")
  if (!is.character(vc)) for (k in names(vc)) cat("  ", k, "sd rows:", rownames(vc[[k]]$sd),
                                                 "Estimate:", signif(vc[[k]]$sd[, 1], 4), "\n")
  cat("varcorr_matrices sqrt diag:", paste(names(varcorr_matrices(fit)),
      sapply(varcorr_matrices(fit), function(m) paste(signif(sqrt(diag(m)), 4), collapse = ","))), "\n")
  ds <- t1(frm_sample(fit, chains = 1, iter = 100, refresh = 0, seed = 4))
  if (is.character(ds)) { cat(ds, "\n"); next }
  lab <- colnames(ds$draws)
  cat("draws labels: ", length(grep("^r_", lab)), "r_,", length(grep("^b\\[", lab)), "b[i];",
      "dups:", unique(lab[duplicated(lab)]), "\n")
  re <- t1(ranef(ds))
  if (is.character(re)) cat("ranef(ds):", re, "\n") else
    for (k in names(re)) cat("  ranef(ds)$", k, " dim ", paste(dim(re[[k]]), collapse = "x"),
                             " NA estimates: ", sum(is.na(re[[k]][, "Estimate", ])), "\n", sep = "")
  mlre <- unclass(ranef(fit))
  cat("  ranef(fit) blocks:", names(mlre), "\n")
  vcd <- t1(VarCorr(ds)); cat("  VarCorr(ds) names:", if (is.character(vcd)) vcd else names(vcd), "\n")
  cf <- t1(coef(ds)); cat("  coef(ds):", if (is.character(cf)) cf else names(cf), "\n")
  if (!is.character(cf)) for (k in names(cf)) cat("   coef(ds)$", k, " NA: ", sum(is.na(cf[[k]])), "\n", sep = "")
}
cat("DONE\n")
