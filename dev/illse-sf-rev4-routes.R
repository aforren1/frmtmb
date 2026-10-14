# Lane illse: copy of dev/rel070-seflags-rev4-routes.R by dev/illse-sf-make.R (lane arm: the lane library over
# rellib-r8; base arm: rellib-r8).
# Reviewer of lane seflags, check of punch 2b: every route a variance
# takes, on fits with an "illcond" direction, against the same route on
# the well-conditioned reparameterization (poly(x, 8), the same
# likelihood and the same fitted values). A route that misses the
# direction gives a ratio well below 1 with no warning.
# Fits: y ~ x + I(x^2) + ... + I(x^8) [+ z + w] [+ (1 | g)], seeds 8 and
# 4 (the punch round's), and the B5 design seed 3. Rows x = 1, 1.5, 2,
# 2.5 (z = 0, w = "a").
#   Rscript dev/seflags-rev4-routes.R <lane|base|op4>
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-sf-rev4-common.R")
suppressMessages({library(emmeans); library(marginaleffects)})
spl <- requireNamespace("frmtmb.spline", quietly = TRUE)
fmt <- function(a, b) {
  if (is.null(a) || is.null(b) || inherits(a, "rev_err") ||
      inherits(b, "rev_err")) return("ERR")
  paste(format(as.numeric(a) / as.numeric(b), digits = 4), collapse = " ")
}
show <- function(lab, a, b) {
  w <- if (is.list(a) && !is.null(a$w)) length(a$w) else NA
  av <- if (is.list(a) && "v" %in% names(a)) a$v else a
  bv <- if (is.list(b) && "v" %in% names(b)) b$v else b
  cat(sprintf("  %-34s %s | warn %s%s\n", lab, fmt(av, bv), w,
              if (inherits(av, "rev_err")) paste(" |", substr(av, 1, 90)) else ""))
}
se_of <- function(x, col) {
  if (inherits(x$v, "rev_err")) return(x$v)
  x$v[[col]]
}
pw <- paste0(" + I(x^", 2:8, ")", collapse = "")
designs <- list(
  list(s = 8, re = TRUE, zw = FALSE),
  list(s = 4, re = FALSE, zw = FALSE),
  list(s = 3, re = TRUE, zw = TRUE))
for (ds in designs) {
  set.seed(ds$s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)),
                  z = stats::rnorm(200),
                  w = factor(sample(c("a", "b"), 200, TRUE)))
  d$y <- sin(3 * d$x) + if (ds$zw) 0.3 * d$z + 0.5 * (d$w == "b") else 0
  d$y <- d$y + stats::rnorm(20, 0, 0.3)[d$g] + stats::rnorm(200, 0, 0.2)
  extra <- paste0(if (ds$zw) " + z + w", if (ds$re) " + (1 | g)")
  fo <- stats::as.formula(paste0("y ~ x", pw, extra))
  fr <- stats::as.formula(paste0("y ~ poly(x, 8)", extra))
  r <- cap(frm(fo, data = d))
  f <- r$v
  g <- suppressWarnings(frm(fr, data = d))
  cat(sprintf("\n== seed %d re %s zw %s | code %d (ref %d) | lost %s | logLik %.4f (ref %.4f)\n",
              ds$s, ds$re, ds$zw, f$opt$convergence, g$opt$convergence,
              lost_of(f), as.numeric(logLik(f)), as.numeric(logLik(g))))
  for (w in r$w) cat("  FIT W:", substr(w, 1, 200), "\n")
  nd <- data.frame(x = c(1, 1.5, 2, 2.5), z = 0, w = factor("a", c("a", "b")),
                   g = d$g[1])
  rf <- if (ds$re) NA else NULL
  both <- function(expr) list(a = cap(eval(expr, list(m = f))),
                              b = cap(eval(expr, list(m = g))))
  # 1 fitted
  x <- both(quote(fitted(m, newdata = nd, re_formula = rf)[, "Est.Error"]))
  show("fitted()", x$a, x$b$v)
  # 2 frm_linpred
  x <- both(quote(frm_linpred(m, newdata = nd, re_formula = rf,
                              se.fit = TRUE)$se.fit))
  show("frm_linpred(se.fit)", x$a, x$b$v)
  # 3 frm_lp_basis, both forms, and frm_joint_cov read whole
  lb <- function(m) {
    b <- frm_lp_basis(m, newdata = nd, re_formula = rf, extra_cov = TRUE)
    A <- as.matrix(b$A)
    J <- frm_joint_cov(m)$V[b$coef_pos, b$coef_pos]
    list(v1 = sqrt(rowSums((A %*% b$V) * A) + b$extra_var),
         v2 = sqrt(diag(A %*% b$V %*% t(A) + as.matrix(b$extra_cov))),
         v0 = sqrt(pmax(rowSums((A %*% b$V) * A), 0)),
         vj = sqrt(pmax(rowSums((A %*% J) * A), 0)))
  }
  la <- cap(lb(f)); lr <- cap(lb(g))
  if (!inherits(la$v, "rev_err") && !inherits(lr$v, "rev_err")) {
    show("lp_basis A V A' + extra_var", la$v$v1, lr$v$v1)
    show("lp_basis A V A' + extra_cov", la$v$v2, lr$v$v2)
    show("lp_basis A V A' alone (BREAKING)", la$v$v0, lr$v$v0)
    show("A frm_joint_cov() A' (whole)", la$v$vj, lr$v$vj)
  } else show("lp_basis", la$v, lr$v)
  # 4 emmeans at the rows
  x <- both(quote(summary(emmeans::emmeans(m, ~ x, at = list(
    x = c(1, 1.5, 2, 2.5), z = 0, w = "a")))$SE))
  show("emmeans()", x$a, x$b$v)
  x <- both(quote(summary(emmeans::contrast(emmeans::emmeans(m, ~ x,
    at = list(x = c(1, 2.5), z = 0, w = "a")), "revpairwise"))$SE))
  show("emmeans contrast 2.5 - 1", x$a, x$b$v)
  # 5 marginaleffects
  mnd <- nd
  x <- both(quote(marginaleffects::predictions(m, newdata = mnd,
                                               re_formula = rf)$std.error[seq_len(nrow(mnd))]))
  show("marginaleffects predictions()", x$a, x$b$v)
  x <- both(quote(marginaleffects::slopes(m, variables = "x", newdata = mnd,
                                          re_formula = rf)$std.error))
  show("marginaleffects slopes(x)", x$a, x$b$v)
  # 6 conditional_effects over x
  ce <- function(m) {
    o <- conditional_effects(m, effects = "x")[[1]]
    o$se__[c(1, 34, 67, 100)]
  }
  x <- both(quote(ce(m)))
  show("conditional_effects(x) rows", x$a, x$b$v)
  # 7 frmtmb.spline
  if (spl) {
    x <- both(quote(frmtmb.spline::frm_curve(m, newdata = nd,
                                             simultaneous = FALSE)$.se))
    show("frm_curve()", x$a, x$b$v)
    x <- both(quote(frmtmb.spline::frm_curve_deriv(m, var = "x", newdata = nd,
                                                   simultaneous = FALSE)$.se))
    show("frm_curve_deriv(x)", x$a, x$b$v)
    x <- both(quote(frmtmb.spline::frm_curve(m, newdata = nd[4, ],
      contrast = nd[1, ], simultaneous = FALSE)$.se))
    show("frm_curve contrast 2.5 - 1", x$a, x$b$v)
  }
  # 8 parameters beside the direction: sigma, sd, and the B5 coefficients
  vc <- function(m) {
    v <- VarCorr(m)
    c(sigma = v$residual__$sd[1, "Est.Error"],
      sd = if (!is.null(v$g)) v$g$sd[1, "Est.Error"] else NA)
  }
  x <- both(quote(vc(m)))
  show("VarCorr() sigma, sd", x$a, x$b$v)
  if (ds$zw) {
    x <- both(quote(hypothesis(m, "b_wb = 0", class = NULL)$hypothesis$Est.Error))
    show("hypothesis(b_wb = 0)", x$a, x$b$v)
    x <- both(quote(sqrt(diag(vcov(m)))[c("z", "wb")]))
    show("vcov() z, wb", x$a, x$b$v)
    x <- both(quote(summary(emmeans::emmeans(m, ~ w, at = list(x = 1.5,
                                                      z = 0)))$SE))
    show("emmeans(~ w) at x = 1.5", x$a, x$b$v)
    x <- both(quote(marginaleffects::avg_comparisons(m, variables = "w",
                                                     re_formula = rf)$std.error))
    show("avg_comparisons(w)", x$a, x$b$v)
  }
  if (ds$re) {
    x <- both(quote({
      rr <- ranef(m, condVar = TRUE)
      attr(rr$g, "condSD")[1:3, 1]
    }))
    show("ranef() Est.Error (3 levels)", x$a, x$b$v)
    x <- both(quote(fitted(m, newdata = nd)[, "Est.Error"]))
    show("fitted() with the group effect", x$a, x$b$v)
  }
}
