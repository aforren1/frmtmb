# Reviewer: recovery under a 3 x 3 correlated (1 | ID | g) block on a
# K = 4 response, reference not first, for both families.
#   Rscript famsimplex-rev-recov.R <dir|ln> <seed0> <nseeds>
# Output: dev/famsimplex-rev-recov-<fam>.log, -<fam>.rds
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
a <- commandArgs(TRUE)
which <- a[1]
s0 <- as.integer(a[2])
ns <- as.integer(a[3])
cat("lib:", find.package("frmtmb"), " fam:", which, " seeds", s0, "to",
    s0 + ns - 1, "\n")
ng <- 30
n <- 600
sdre <- c(0.5, 0.4, 0.3)
Cre <- matrix(c(1, 0.5, -0.3, 0.5, 1, 0.2, -0.3, 0.2, 1), 3)
b0 <- c(a = 0.3, b = -0.2, c = 0.5)
b1 <- c(a = 0.5, b = -0.4, c = 0.2)
truth <- c(mua_Intercept = 0.3, mub_Intercept = -0.2, muc_Intercept = 0.5,
           mua_x = 0.5, mub_x = -0.4, muc_x = 0.2)
out <- list()
for (s in s0 + seq_len(ns) - 1L) {
  set.seed(s)
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng), length.out = n)))
  U <- matrix(rnorm(ng * 3), ng) %*% chol(diag(sdre) %*% Cre %*% diag(sdre))
  M <- cbind(b0[1] + b1[1] * d$x, b0[2] + b1[2] * d$x, b0[3] + b1[3] * d$x) +
    U[as.integer(d$g), ]
  # reference category d is the LAST column
  if (which == "dir") {
    Y <- fs_rdirichlet(fs_softmax(cbind(M, 0)) * 12)
  } else {
    Sg <- diag(c(0.5, 0.6, 0.7)) %*%
      matrix(c(1, 0.3, 0.1, 0.3, 1, -0.2, 0.1, -0.2, 1), 3) %*%
      diag(c(0.5, 0.6, 0.7))
    Y <- fs_rlogistic_normal(M, Sg, ref = 4L)
  }
  colnames(Y) <- c("a", "b", "c", "d")
  d$Y <- Y
  fam <- if (which == "dir") dirichlet(refcat = "d") else
    logistic_normal(refcat = "d")
  w <- character(0)
  t0 <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
    tryCatch(frm(bf(Y ~ x + (1 | ID | g)), family = fam, data = d),
             error = function(e) conditionMessage(e)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  el <- proc.time()[["elapsed"]] - t0
  if (is.character(fit)) {
    out[[as.character(s)]] <- list(err = fit)
    next
  }
  fe <- fixef(fit)
  vc <- VarCorr(fit)[["g"]]
  out[[as.character(s)]] <- list(
    est = fe[names(truth), "Estimate"], se = fe[names(truth), "Est.Error"],
    sd = vc$sd[, "Estimate"],
    cor = c(vc$cor["mub_Intercept", "Estimate", "mua_Intercept"],
            vc$cor["muc_Intercept", "Estimate", "mua_Intercept"],
            vc$cor["muc_Intercept", "Estimate", "mub_Intercept"]),
    conv = fit$opt$convergence, warn = w, sec = el)
}
saveRDS(out, sub("[.]R$", paste0("-", which, ".rds"),
                 "C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-rev-recov.R"))
ok <- Filter(function(r) is.null(r$err), out)
cat("fitted", length(ok), "of", length(out), "; errors:",
    sum(vapply(out, function(r) !is.null(r$err), NA)), "\n")
cat("nonzero convergence:", sum(vapply(ok, `[[`, 0, "conv") != 0), "\n")
cat("fits with a warning:", sum(vapply(ok, function(r) length(r$warn) > 0,
                                       NA)), "\n")
uw <- unique(unlist(lapply(ok, `[[`, "warn")))
if (length(uw)) cat("warnings seen:", substr(uw, 1, 200), sep = "\n  ")
E <- do.call(rbind, lapply(ok, `[[`, "est"))
S <- do.call(rbind, lapply(ok, `[[`, "se"))
nk <- nrow(E)
tab <- data.frame(truth = truth, mean = colMeans(E),
                  bias_over_mcse = (colMeans(E) - truth) /
                    (apply(E, 2, sd) / sqrt(nk)),
                  sd_est = apply(E, 2, sd), mean_se = colMeans(S),
                  cover95 = colMeans(abs(E - rep(truth, each = nk)) <
                                       1.959964 * S))
print(round(tab, 3))
SD <- do.call(rbind, lapply(ok, `[[`, "sd"))
CR <- do.call(rbind, lapply(ok, `[[`, "cor"))
cat("sd truth", sdre, " mean est", round(colMeans(SD), 3), "\n")
cat("cor truth (b,a) 0.5 (c,a) -0.3 (c,b) 0.2; mean est",
    round(colMeans(CR), 3), "; MC SE", round(apply(CR, 2, sd) / sqrt(nk), 3),
    "\n")
cat("median seconds:", median(vapply(ok, `[[`, 0, "sec")), "\n")
