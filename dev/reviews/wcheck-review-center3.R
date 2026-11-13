# Review of lane wcheck: case R1 of wcheck-review-center2.R. Is the
# interior maximum (b_facb 0.626) the global one? Maximize the objective
# with b_facb held below a cap, from several starts, and profile b_facb
# far out with the rest re-optimized.
#   Rscript dev/reviews/wcheck-review-center3.R <lib|base>
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
d$w <- ifelse(d$fac == "a", 0.8, 1.05)
q <- function(e) suppressWarnings(suppressMessages(e))
fo <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac)
for (pr in list(set_prior("normal(0, 1)", class = "Intercept",
                          dpar = "sigma"), NULL)) {
  cat("==== prior:", if (is.null(pr)) "none" else "normal(0, 1) on Intercept",
      "\n")
  ob <- q(frm(fo, family = gaussian(), data = d, prior = pr,
              dry_run = "objective"))
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p0 <- ob$par
  nm <- names(p0)
  isd <- which(nm == "betad")
  j <- isd[2L]
  cat("par:", nm, "\n")
  for (cap in c(-5, -10, -20)) {
    up <- rep(Inf, length(p0)); up[j] <- cap
    best <- -Inf
    for (s in list(c(0, cap), c(2, cap - 5), c(5, cap - 10), c(-2, cap),
                   c(10, cap - 20))) {
      st <- p0; st[isd] <- s
      op <- tryCatch(nlminb(st, ob$fn, ob$gr, upper = up,
                            control = list(iter.max = 3000,
                                           eval.max = 6000)),
                     error = function(e) NULL)
      if (is.null(op) || !is.finite(op$objective)) next
      if (-op$objective > best) {
        best <- -op$objective
        bp <- op$par
      }
    }
    cat(sprintf("  b_facb <= %g: best -fn %.4f at betad %s theta %s\n", cap,
                best, paste(format(bp[isd], digits = 5), collapse = " "),
                paste(format(bp[nm == "theta"], digits = 5),
                      collapse = " ")))
  }
  # the profile far out, warm-started outward
  op <- nlminb(p0, ob$fn, ob$gr)
  p <- op$par
  cat(sprintf("  unconstrained: -fn %.4f betad %s\n", -op$objective,
              paste(format(p[isd], digits = 5), collapse = " ")))
  free <- setdiff(seq_along(p0), j)
  for (t in c(-10, -15, -20, -25, -30, -35, -40)) {
    best <- -Inf
    for (s in list(p, replace(p, isd[1L], 5), replace(p, isd[1L], 15))) {
      s[j] <- t
      fr <- function(r) {
        q2 <- s; q2[free] <- r
        v <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
        if (!is.finite(v)) 1e10 else v
      }
      o2 <- nlminb(s[free], fr)
      if (-o2$objective > best) {best <- -o2$objective; bb <- o2$par}
    }
    p[free] <- bb
    p[j] <- t
    cat(sprintf("    b_facb %5.1f  -fn %12.4f  b0 %8.3f  Ic %8.3f\n", t, best,
                p[isd[1L]], p[isd[1L]] + 0.5 * t))
  }
}
cat("DONE\n")
