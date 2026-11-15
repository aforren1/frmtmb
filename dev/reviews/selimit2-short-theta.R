# Short review of lane selimit2's fix round, F2: can a finite-MLE mixing
# weight pass se_limit_joint_cert()'s far-out check? Two gaussian
# components 16 sd apart, theta1 ~ z: rows with z > 0 from one component,
# z < 0 from the other, plus a reversed pair (each row from the other
# side's component) at z = -delta/2 and +delta/2. No boundary puts both
# on their own side, so the MLE is finite; but the pair moves by about
# delta of the largest move, below the 1e-3 the far check saturates.
#   Rscript selimit2-short-theta.R <base|lane> [delta ...]
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
deltas <- if (length(a) > 1) as.numeric(a[-1]) else c(1e-4, 1e-6)
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "lib", dirname(find.package("frmtmb")), "\n")
fitw <- function(...) {
  ws <- character(0)
  ms <- character(0)
  fit <- withCallingHandlers(frm(...), warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    ms <<- c(ms, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  list(fit = fit, ws = ws, ms = ms)
}
for (delta in deltas) {
  set.seed(3)
  z <- c(-stats::runif(60, 0.01, 1), stats::runif(60, 0.01, 1))
  y <- c(stats::rnorm(60, -4, 0.5), stats::rnorm(60, 4, 0.5))
  # the reversed pair: a +4 row just below 0, a -4 row just above
  z <- c(z, -delta / 2, delta / 2)
  y <- c(y, 4, -4)
  d <- data.frame(y = y, z = z)
  for (st in "default") {
    r <- if (st == "default") {
      fitw(bf(y ~ 1, theta1 ~ z), data = d,
           family = mixture(gaussian(), gaussian()))
    } else {
      tpl <- frmtmb::par_template(bf(y ~ 1, theta1 ~ z), data = d,
                                  family = mixture(gaussian(), gaussian()))
      fitw(bf(y ~ 1, theta1 ~ z), data = d,
           family = mixture(gaussian(), gaussian()), start = tpl)
    }
    fit <- r$fit
    p <- fit$opt$par
    cat(sprintf("\n== delta %g start %s | code %d | obj %.10f\n", delta, st,
                fit$opt$convergence, fit$opt$objective))
    cat("   par:", paste(frmtmb:::outer_par_names(fit), format(p, digits = 6),
                         collapse = "; "),
        "\n")
    for (w in r$ws) cat("   W:", substr(w, 1, 240), "\n")
    for (m in r$ms) cat("   M:", substr(m, 1, 200), "\n")
    li <- fit$cache$se_limit
    if (!is.null(li)) {
      cat("   limit:", paste(names(li$kind), li$kind, collapse = "; "),
          "| joint", paste(li$joint, collapse = ","), "\n")
    }
    # the objective along the theta coefficients scaled up, the rest held
    it <- grep("^theta1_", frmtmb:::outer_par_names(fit))
    if (length(it) == 2) {
      f0 <- fit$obj$fn(p)
      for (k in c(2, 10, 100, 1e3, 1e4, 1e6)) {
        q <- p
        q[it] <- p[it] * k
        cat(sprintf("   theta x %g: objective change %.4g\n", k,
                    fit$obj$fn(q) - f0))
      }
    }
  }
}
