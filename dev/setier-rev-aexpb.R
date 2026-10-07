# Reviewer of lane setier: curved ridges that se_curvature_real() can
# "rescue". Lane fixes' table rows, with the data seeded (the table drew
# them unseeded, which is lane nanse review's m10):
#   aexpb  y ~ a + exp(b), a ~ 1 + x, b ~ 1: only a_(Intercept) + exp(b) is
#          identified (a curved ridge in (a_(Intercept), b))
#   ab     y ~ a * b, a ~ 1, b ~ 1: only a * b is identified
#   aexpbre the aexpb ridge with (1 | g) on a (finite-difference path)
# Per seed: warnings, finite SEs, lost, and the objective along the
# exact ridge (it must not move).
#   Rscript dev/setier-rev-aexpb.R lane|base [seeds]
args <- commandArgs(TRUE)
arm <- args[1]
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:20
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", arm, find.package("frmtmb"), "\n")
catch <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) invokeRestart("muffleMessage"))
  list(value = v, w = w)
}
kind <- function(m) {
  if (grepl("^Standard errors are not available", m)) return("SE")
  if (grepl("are not identified: at the optimum", m)) return("NLFLAT")
  if (grepl("^Optimizer did not report|gradient", m)) return("CONV")
  "OTHER"
}
tab <- list()
for (s in seeds) {
  set.seed(s)
  n <- 60
  dd <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:6, 10)))
  dd$y <- 3 + 0.5 * dd$x + 0.3 * dd$z + rnorm(6, 0, 0.5)[dd$g] +
    rnorm(n, 0, 0.4)
  for (des in c("aexpb", "ab", "aexpbre")) {
    r <- tryCatch(catch(switch(des,
      aexpb = frm(bf(y ~ a + exp(b), a ~ 1 + x, b ~ 1, nl = TRUE), data = dd,
                  start = list(beta = c(1, 0, -1))),
      ab = frm(bf(y ~ a * b, a ~ 1, b ~ 1, nl = TRUE), data = dd,
               start = list(beta = c(1, 1))),
      aexpbre = frm(bf(y ~ a + exp(b), a ~ 1 + x + (1 | g), b ~ 1,
                       nl = TRUE), data = dd,
                    start = list(beta = c(1, 0, -1))))),
      error = function(e) e)
    if (inherits(r, "error")) {
      cat(des, s, "ERROR", conditionMessage(r), "\n"); next
    }
    f <- r$value
    nm <- ns$outer_par_names(f)
    p <- f$opt$par
    se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
    coef <- grepl("^(a|b)_", nm)
    f0 <- f$obj$fn(p)
    # the exact ridge, moved by 0.5 and by 2 along it
    mv <- vapply(c(0.5, 2), function(h) {
      q <- p
      if (des == "ab") {
        ia <- match("a_(Intercept)", nm); ib <- match("b_(Intercept)", nm)
        q[ia] <- p[ia] * h; q[ib] <- p[ib] / h
      } else {
        ia <- match("a_(Intercept)", nm); ib <- match("b_(Intercept)", nm)
        # keep a + exp(b): b moves by log(h), a takes up the difference
        q[ib] <- p[ib] + log(h)
        q[ia] <- p[ia] + exp(p[ib]) - exp(q[ib])
      }
      f$obj$fn(q) - f0
    }, 0)
    lost <- ns$sdr_of(f)$se_lost
    k <- if (length(r$w)) paste(vapply(r$w, kind, ""), collapse = "+") else
      "none"
    cat(sprintf(paste0("%-8s seed %2d code %d | warn %-12s | finite coef SE ",
                       "%d/%d (%s) | lost %s | ridge nll change %.2g %.2g\n"),
                des, s, f$opt$convergence, k, sum(is.finite(se[coef])),
                sum(coef), paste(signif(se[coef], 3), collapse = " "),
                if (length(lost)) paste(names(lost), lost, sep = ":",
                                        collapse = ",") else "none",
                mv[1], mv[2]))
    tab[[length(tab) + 1L]] <- data.frame(
      des = des, seed = s, code = f$opt$convergence, warn = k,
      finite = sum(is.finite(se[coef])), ncoef = sum(coef),
      ridge = max(abs(mv)))
  }
}
X <- do.call(rbind, tab)
cat("\nsummary: fits whose ridge is flat (|change| < 1e-6), by design\n")
X$flatridge <- X$ridge < 1e-6
for (des in unique(X$des)) {
  x <- X[X$des == des & X$flatridge & X$code == 0, ]
  cat(sprintf(paste0("%-8s flat-ridge code-0 fits %d | all coef SE finite %d",
                     " | of those, silent %d\n"),
              des, nrow(x), sum(x$finite == x$ncoef),
              sum(x$finite == x$ncoef & x$warn == "none")))
}
