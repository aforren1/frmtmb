# Review of lane wcheck: a sigma ~ 1 + fac model with an ordinary
# interior maximum that the cell rule refuses. frmtmb's class Intercept
# prior is evaluated at the centered intercept Ic = b0 + mean(facb) *
# b_facb (priors.R, "Where an intercept prior lands"), so b_facb -> -Inf
# with b0 fixed drives Ic -> -Inf, where a normal() prior stops it. The
# ray that keeps Ic fixed (b0 -> +Inf) gains E_b / n_b - W_a / n_a per
# unit of eta_b (per row of b): bounded when the mean excess of cell b is
# below the mean weight of cell a.
#   Rscript dev/reviews/wcheck-review-center2.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(20261111)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                fac = factor(rep(c("a", "b"), 40)))
d$y <- 0.5 + 0.4 * d$x + rnorm(n, 0, ifelse(d$fac == "a", 0.6, 1.2))
q <- function(e) suppressWarnings(suppressMessages(e))
P <- function(...) set_prior(..., dpar = "sigma")
fo <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac)

run <- function(label, wa, wb, prior) {
  dd <- d
  dd$w <- ifelse(dd$fac == "a", wa, wb)
  Ea <- sum(dd$w[dd$fac == "a"] - 1); Eb <- sum(dd$w[dd$fac == "b"] - 1)
  cat(sprintf("==== %s: w_a %g w_b %g | E_a %.2f E_b %.2f total %.2f | E_b/n_b %.3f vs mean w_a %.3f\n",
              label, wa, wb, Ea, Eb, Ea + Eb, Eb / 40, wa))
  ws <- character()
  r <- tryCatch(withCallingHandlers(
    frm(fo, family = gaussian(), data = dd, prior = prior),
    warning = function(w) {ws <<- c(ws, conditionMessage(w));
      invokeRestart("muffleWarning")},
    message = function(m) invokeRestart("muffleMessage")),
    error = function(e) e)
  if (inherits(r, "error")) {
    cat("  frm():", substr(gsub("\\s+", " ", conditionMessage(r)), 1, 400),
        "\n")
  } else {
    cat(sprintf("  frm(): logLik %.4f conv %d betad %s theta %s warnings %d\n",
                as.numeric(logLik(r)), r$opt$convergence,
                paste(format(r$estimates$betad, digits = 5), collapse = " "),
                paste(format(r$estimates$theta, digits = 5), collapse = " "),
                length(ws)))
    if (length(ws)) cat("   ", substr(ws, 1, 100), sep = "\n    ")
  }
  # the objective (with the prior), and a profile in b_facb
  ob <- q(frm(fo, family = gaussian(), data = dd, prior = prior,
              dry_run = "objective"))
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p <- ob$par
  isd <- which(names(p) == "betad")
  best <- nlminb(p, ob$fn, ob$gr)
  cat(sprintf("  nlminb on the objective: -fn %.4f conv %d max|grad| %.2e betad %s\n",
              -best$objective, best$convergence, max(abs(ob$gr(best$par))),
              paste(format(best$par[isd], digits = 5), collapse = " ")))
  j <- isd[2L]
  free <- setdiff(seq_along(p), j)
  ts <- c(2, 1, 0, -1, -2, -3, -4, -6, -8, -10, -12, -14)
  p <- best$par
  # from the maximum outward on each side
  res <- numeric(length(ts))
  ord <- order(abs(ts - best$par[j]))
  for (k in ord) {
    pk <- best$par
    pk[j] <- ts[k]
    fr <- function(r) {
      q2 <- pk; q2[free] <- r
      v <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
      if (!is.finite(v)) 1e10 else v
    }
    op <- nlminb(pk[free], fr)
    res[k] <- -op$objective
  }
  cat("  profile in b_facb, the rest re-optimized from the maximum:\n")
  cat(sprintf("    b_facb %5.1f  -fn %12.4f\n", ts, res), sep = "")
}

pn <- P("normal(0, 1)", class = "Intercept")
run("R1 normal(0, 1) on Intercept, bounded", 0.8, 1.05, pn)
run("R1b normal(0, 1) on Intercept, bounded, E_b 8", 0.8, 1.2, pn)
run("R2 no prior (truly unbounded along b_facb)", 0.8, 1.05, NULL)
run("R3 normal(0, 1) on Intercept, E_b/n_b > w_a (unbounded)", 0.3, 1.5, pn)
cat("DONE\n")
