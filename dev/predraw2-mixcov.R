# Lane predraw2: predict() on the weakly identified mixtures of backlog
# 4882 against the KNOWN predictive law of the simulated data (the
# generators of dev/ppdraw-sweep-gen.R, 300 rows), over seeds. Per fit and
# training row the exact probability that the true predictive law puts in
# predict()'s central 90% and 98% intervals, averaged over rows; for the
# ordinal mixture the mean absolute error of the predicted category
# shares against the true probabilities. Also predict()'s mean against
# fitted()'s (z, Monte Carlo standard errors).
# Usage: Rscript dev/predraw2-mixcov.R <lib or "base"> <design> <seeds>
#        <ndraws> <out.tsv>
#   design: mixture_gaussian:boundary, mixture_poisson:boundary,
#           mixture_cumulative:interior, ... ; seeds like 1:40
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
n <- 300L
source("dev/ppdraw-sweep-gen.R")
design <- strsplit(args[2], ":")[[1]]
seeds <- eval(parse(text = args[3]))
B <- as.integer(args[4])
out <- args[5]
fam <- design[1]
case <- design[2]
truth_row <- function(d) {
  # the generators' own laws, row by row
  x <- d$x
  inner <- case == "interior"
  switch(fam,
    mixture_gaussian = if (inner) {
      list(type = "mixn", w = c(0.6, 0.4), mu = cbind(-1 + 0.5 * x,
                                                       3 + 0.5 * x), sd = 1)
    } else {
      list(type = "mixn", w = 1, mu = cbind(0.5 * x), sd = 1)
    },
    mixture_poisson = if (inner) {
      list(type = "mixp", w = c(0.6, 0.4),
           lam = cbind(exp(0.7 + 0.3 * x), exp(2.3 + 0.3 * x)))
    } else {
      list(type = "mixp", w = 1, lam = cbind(exp(1 + 0.3 * x)))
    },
    mixture_cumulative = {
      cuts <- c(-Inf, -1, 0.5, 2, Inf)
      sh <- if (inner) c(0, 3) else 0
      w <- if (inner) c(0.6, 0.4) else 1
      P <- Reduce(`+`, lapply(seq_along(sh), function(k) {
        w[k] * t(vapply(x, function(xi) diff(stats::plogis(cuts - xi - sh[k])),
                        numeric(4)))
      }))
      list(type = "cat", P = P)
    })
}
mass_in <- function(tr, i, lo, hi) {
  if (tr$type == "mixn") {
    sum(tr$w * (stats::pnorm(hi, tr$mu[i, ], tr$sd) -
                  stats::pnorm(lo, tr$mu[i, ], tr$sd)))
  } else {
    sum(tr$w * (stats::ppois(hi, tr$lam[i, ]) -
                  stats::ppois(ceiling(lo) - 1, tr$lam[i, ])))
  }
}
rows <- list()
for (seed in seeds) {
  d <- gen(fam, case, seed)
  s <- spec(fam)
  fit <- tryCatch(suppressMessages(suppressWarnings(
    frm(s[[1]], family = s[[2]], data = d))), error = function(e) NULL)
  if (is.null(fit)) next
  said <- character(0)
  set.seed(1)
  p <- tryCatch(withCallingHandlers(
    predict(fit, ndraws = B, summary = FALSE),
    warning = function(w) {
      said <<- c(said, class(w)[1], conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      said <<- c(said, class(m)[1], conditionMessage(m))
      invokeRestart("muffleMessage")
    }), error = function(e) NULL)
  if (is.null(p)) next
  chk <- if (is.environment(fit$cache)) fit$cache[["draw_check"]]
  kk <- regmatches(said, regexpr("Pareto k (is )?[0-9.]+|too few", said))
  kk <- if (length(kk)) sub("Pareto k (is )?", "", kk[1]) else ""
  tr <- truth_row(d)
  fe <- suppressWarnings(fitted(fit))
  r <- list(seed = seed, fire = isTRUE(chk$fire),
            excess = chk$excess %||% NA_real_, rel = chk$rel %||% NA_real_,
            k = kk,
            held = any(said == "frmtmb_draw_held"),
            fallback = any(grepl("could not be recovered", said)))
  if (tr$type == "cat") {
    sh <- t(apply(p, 2L, function(col) tabulate(col, 4L) / B))
    r$mae_pred <- mean(abs(sh - tr$P))
    r$mae_fitted <- mean(abs(fe[, "Estimate", ] - tr$P))
    P <- fe[, "Estimate", ]
    r$z <- max(abs((sh - P) / sqrt(P * (1 - P) / B)))
    r$cov90 <- r$cov98 <- r$w90 <- r$w98 <- NA_real_
  } else {
    q <- apply(p, 2L, stats::quantile, c(0.01, 0.05, 0.95, 0.99),
               names = FALSE)
    r$cov90 <- mean(vapply(seq_len(nrow(d)), function(i) {
      mass_in(tr, i, q[2, i], q[3, i])
    }, 0))
    r$cov98 <- mean(vapply(seq_len(nrow(d)), function(i) {
      mass_in(tr, i, q[1, i], q[4, i])
    }, 0))
    r$w90 <- mean(q[3, ] - q[2, ])
    r$w98 <- mean(q[4, ] - q[1, ])
    m <- rowMeans(p)
    r$z <- (mean(m) - mean(fe[, "Estimate"])) / (stats::sd(m) / sqrt(B))
    tm <- if (tr$type == "mixn") tr$mu %*% tr$w else tr$lam %*% tr$w
    r$mae_pred <- mean(abs(colMeans(p) - tm))
    r$mae_fitted <- mean(abs(fe[, "Estimate"] - tm))
  }
  # the plug-in predictive (no estimation error) on the same rows
  set.seed(1)
  p0 <- suppressWarnings(predict(fit, ndraws = B, summary = FALSE,
                                 propagate_error = FALSE))
  if (tr$type == "cat") {
    sh0 <- t(apply(p0, 2L, function(col) tabulate(col, 4L) / B))
    r$mae_plug <- mean(abs(sh0 - tr$P))
    r$cov90_plug <- r$cov98_plug <- NA_real_
  } else {
    q0 <- apply(p0, 2L, stats::quantile, c(0.01, 0.05, 0.95, 0.99),
                names = FALSE)
    r$cov90_plug <- mean(vapply(seq_len(nrow(d)), function(i) {
      mass_in(tr, i, q0[2, i], q0[3, i])
    }, 0))
    r$cov98_plug <- mean(vapply(seq_len(nrow(d)), function(i) {
      mass_in(tr, i, q0[1, i], q0[4, i])
    }, 0))
    r$mae_plug <- mean(abs(colMeans(p0) - tm))
  }
  rows[[length(rows) + 1L]] <- as.data.frame(r)
  cat(sprintf("seed %3d fire %-5s excess %7.3f k %-5s held %-5s fb %-5s cov90 %.3f cov98 %.3f w98 %.3f z %7.2f mae %.4f (fitted %.4f) plug cov90 %.3f cov98 %.3f mae %.4f\n",
              seed, r$fire, r$excess, r$k, r$held, r$fallback, r$cov90,
              r$cov98, r$w98, r$z, r$mae_pred, r$mae_fitted, r$cov90_plug,
              r$cov98_plug, r$mae_plug))
}
res <- do.call(rbind, rows)
utils::write.table(res, out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("rows", nrow(res), "\n")
