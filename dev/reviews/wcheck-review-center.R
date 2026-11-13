# Review of lane wcheck: is sigma ~ 1 + fac with a normal() prior on the
# (centered) sigma Intercept really unbounded, as the cell rule says?
# frmtmb evaluates a class Intercept density at b0 + means_X'b (priors.R
# "Where an intercept prior lands"), so the cell direction b_facb -> -Inf
# with b0 fixed moves the prior's argument too.
#   Rscript dev/reviews/wcheck-review-center.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(20261021)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                fac = factor(rep(c("a", "b"), 40)),
                w = runif(n, 0.5, 2), g = factor(rep(1:20, each = 4)))
eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$y <- eta + rnorm(n, 0, 0.8)
cat(sprintf("excess: all %.3f | a %.3f | b %.3f | mean w a %.3f b %.3f\n",
            sum(d$w) - n, sum(d$w[d$fac == "a"] - 1),
            sum(d$w[d$fac == "b"] - 1), mean(d$w[d$fac == "a"]),
            mean(d$w[d$fac == "b"])))
q <- function(e) suppressWarnings(suppressMessages(e))
P <- function(...) set_prior(..., dpar = "sigma")

prof <- function(label, fo, prior, data = d) {
  cat("====", label, "\n")
  ob <- q(frm(fo, family = gaussian(), data = data, prior = prior,
              dry_run = "objective"))
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p0 <- ob$par
  nm <- names(p0)
  isd <- which(nm == "betad")
  cat("par:", nm, "\n")
  # the full fit from the default start and a few others
  best <- NULL
  for (s in list(p0, replace(p0, isd, c(-1, -1)), replace(p0, isd, c(0, -3)),
                 replace(p0, isd, c(1, -6)))) {
    op <- tryCatch(nlminb(s, ob$fn, ob$gr, control = list(iter.max = 2000,
                                                         eval.max = 4000)),
                   error = function(e) NULL)
    if (is.null(op)) next
    cat(sprintf("  start betad %s -> -fn %.4f conv %d betad %s\n",
                paste(format(s[isd], digits = 3), collapse = " "),
                -op$objective, op$convergence,
                paste(format(op$par[isd], digits = 5), collapse = " ")))
    if (is.null(best) || op$objective < best$objective) best <- op
  }
  gr <- ob$gr(best$par)
  cat(sprintf("  best -fn %.4f max|grad| %.2e\n", -best$objective,
              max(abs(gr))))
  # profile in b_facb with every other parameter re-optimized
  j <- isd[2L]
  free <- setdiff(seq_along(p0), j)
  p <- best$par
  ts <- c(0, -2, -5, -10, -15, -20, -25, -30)
  out <- vapply(ts, function(t) {
    p[j] <- t
    fr <- function(r) {
      q2 <- p; q2[free] <- r
      v <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
      if (!is.finite(v)) 1e10 else v
    }
    op <- nlminb(p[free], fr)
    if (op$objective < 1e9) p[free] <<- op$par
    c(-op$objective, p[isd[1L]])
  }, c(0, 0))
  cat("  profile in b_facb (others re-optimized):\n")
  cat(sprintf("    b_facb %6.1f  -fn %12.4f  b0 %8.3f  Ic %8.3f\n", ts,
              out[1, ], out[2, ], out[2, ] + mean(d$fac == "b") * ts),
      sep = "")
  # and along b_facb with b0 held at the best fit (the lane's direction)
  p <- best$par
  out2 <- vapply(ts, function(t) {
    p[j] <- t
    fr <- function(r) {
      q2 <- p; q2[setdiff(free, isd[1L])] <- r
      v <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
      if (!is.finite(v)) 1e10 else v
    }
    op <- nlminb(p[setdiff(free, isd[1L])], fr)
    -op$objective
  }, 0)
  cat("  along b_facb, b0 held at the best fit:\n")
  cat(sprintf("    b_facb %6.1f  -fn %12.4f\n", ts, out2), sep = "")
  # the check's verdict
  v <- tryCatch({q(frm(fo, family = gaussian(), data = data, prior = prior,
                       dry_run = NULL)); "fit"},
                error = function(e) paste("ERROR:", substr(gsub("\\s+", " ",
                  conditionMessage(e)), 1, 140)))
  cat("  frm():", v, "\n")
}

fo1 <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac)
prof("sigma ~ 1 + fac, normal(0, 1) on Intercept", fo1,
     P("normal(0, 1)", class = "Intercept"))
prof("sigma ~ 1 + fac, normal(0, 5) on Intercept", fo1,
     P("normal(0, 5)", class = "Intercept"))
fo2 <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac + x)
prof("sigma ~ 1 + fac + x, normal(0, 1) on Intercept", fo2,
     P("normal(0, 1)", class = "Intercept"))
# weights only on cell b: cell a mean weight 1, cell b mean weight 1.8.
# The ray Ic fixed, b0 -> +Inf, b_facb -> -Inf gains E_b - n_b * wbar_a
# per unit of b_facb: 40 * 0.8 - 40 * 1 < 0, still bounded
d2 <- d
d2$w <- ifelse(d2$fac == "b", 1.8, 1)
prof("cell b weights 1.8, cell a 1: normal(0, 1) on Intercept", fo1,
     P("normal(0, 1)", class = "Intercept"), data = d2)
# and where that ray does gain: cell b weights 2.5 (E_b 60 > 40)
d3 <- d
d3$w <- ifelse(d3$fac == "b", 2.5, 1)
prof("cell b weights 2.5, cell a 1: normal(0, 1) on Intercept", fo1,
     P("normal(0, 1)", class = "Intercept"), data = d3)
cat("DONE\n")
