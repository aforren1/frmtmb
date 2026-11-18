# Short check of lane squant fix round 1, item 2: the review's
# left-censored fits (squant-review-tape4.R): tape gradient against a
# Richardson finite-difference gradient, and the fit at a maximum (the
# finite-difference gradient relative to the objective, and a BFGS / nlminb
# polish on fn moves neither logLik nor the estimates).
# Run (worktree root):
#   Rscript dev/reviews/squant-short-fit.R > dev/reviews/squant-short-fit.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")),
    "| elf_pos0:", exists("elf_pos0", asNamespace("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
# central difference with one Richardson step; an absolute step, because
# the kernel width lam sigma (about 0.1 here) is far below |intercept|
fd_grad <- function(fn, par, h0 = 1e-4) {
  vapply(seq_along(par), function(k) {
    D <- function(h) {
      u <- par
      d <- par
      u[k] <- u[k] + h
      d[k] <- d[k] - h
      (fn(u) - fn(d)) / (2 * h)
    }
    h <- h0
    (4 * D(h / 2) - D(h)) / 3
  }, 0)
}
rows <- list()
for (L in c(95, 90, 80, 60)) {
  for (p in c(0.1, 0.5)) {
    for (lam in c(0.25, 1)) {
      set.seed(20261118)
      n <- 200
      d <- data.frame(x = stats::runif(n))
      d$y <- 100 + 2 * d$x + stats::rnorm(n)
      d$c <- 0
      d$c[1:6] <- -1
      d$y[1:6] <- L
      ws <- character()
      fit <- withCallingHandlers(
        frm(bf(y | cens(c) ~ x, quantile = p),
            family = quantile_elf(smoothing = lam), data = d),
        warning = function(w) {
          ws <<- c(ws, substr(conditionMessage(w), 1, 80))
          invokeRestart("muffleWarning")
        })
      par <- fit$opt$par
      dpv <- lapply(ns$eval_dpars(fit)[[1L]], as.numeric)
      tt <- (L - dpv$mu[1:6]) / (lam * dpv$sigma[1:6])
      f0 <- fit$obj$fn(par)
      gt <- as.numeric(fit$obj$gr(par))
      gf <- fd_grad(fit$obj$fn, par)
      agree <- max(abs(gt - gf)) / max(abs(gf), 1e-8)
      # tape against fd at a point away from the optimum, where the
      # gradient is not ~0 and the agreement is a real check
      p2 <- par + rep_len(c(0.05, -0.05, 0.02), length(par))
      gt2 <- as.numeric(fit$obj$gr(p2))
      gf2 <- fd_grad(fit$obj$fn, p2)
      agree2 <- max(abs(gt2 - gf2) / pmax(abs(gf2), 1e-8))
      o <- stats::optim(par, fit$obj$fn, function(q) fd_grad(fit$obj$fn, q),
                        method = "BFGS",
                        control = list(reltol = 1e-15, maxit = 1000))
      nb <- stats::nlminb(par, fit$obj$fn, fit$obj$gr,
                          control = list(rel.tol = 1e-15, x.tol = 0,
                                         eval.max = 2000, iter.max = 1000))
      H <- stats::optimHess(par, fit$obj$fn, fit$obj$gr)
      ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
      rows[[length(rows) + 1L]] <- c(max(abs(gf)) / abs(f0), agree2)
      cat(sprintf(paste0(
        "L %g p %.1f lam %.2f: conv %d warn %d | t %.0f to %.0f | ",
        "max|tape g| %.3g max|fd g| %.3g, |fd g|/|f| %.2g, tape-fd %.2g | ",
        "off-opt tape/fd rel %.2g | dlogLik BFGS %.2g nlminb %.2g | ",
        "max|dpar| BFGS %.2g nlminb %.2g | Hess eig min %.3g%s\n"),
        L, p, lam, fit$opt$convergence, length(ws), min(tt), max(tt),
        max(abs(gt)), max(abs(gf)), max(abs(gf)) / abs(f0), agree, agree2,
        f0 - o$value, f0 - nb$objective, max(abs(o$par - par)),
        max(abs(nb$par - par)), min(ev),
        if (length(ws)) paste0(" | ", paste(unique(ws), collapse = " | "))
        else ""))
    }
  }
}
r <- do.call(rbind, rows)
cat(sprintf("SUMMARY: max |fd grad| / |objective| %.3g; max off-optimum tape/fd relative difference %.3g over %d fits\n",
            max(r[, 1]), max(r[, 2]), nrow(r)))
cat("DONE\n")
