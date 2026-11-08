# Item 1 diagnosis of the rows the sweep leaves short: the linear spelling's
# default, the body's default, and the working start the body reads. Data as
# dev/fitedge-gpbody-sweep.R. Args: arm.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
ns <- asNamespace("frmtmb")
curves <- c("sin", "wig", "exp", "step")
mk <- function(seed, n, ci, body) {
  set.seed(18000 + 100 * seed + n + ci)
  x <- sort(runif(n, 0, 10)); z <- rnorm(n); w <- runif(n, 1, 2)
  f <- switch(curves[ci], sin = sin(x), wig = 2 * sin(3 * x) + x,
              exp = exp(0.3 * x), step = sign(sin(1.5 * x)))
  y <- switch(body,
    A = f + 0.3 * z + 0.1 * rnorm(n),
    B = rpois(n, exp(1 + 0.5 * f / max(abs(f)) + 0.3 * z)),
    C = exp(0.5 + 0.5 * f / max(abs(f))) + 0.3 * z + 0.1 * rnorm(n),
    D = 1 + f * w + 0.1 * rnorm(n))
  data.frame(x, z, w, y)
}
spell <- list(
  A = list(nl = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE), lin = y ~ z + gp(x), fam = gaussian()),
  B = list(nl = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE), lin = y ~ z + gp(x), fam = poisson()),
  C = list(nl = bf(y ~ exp(a) + b * z, a ~ gp(x), b ~ 1, nl = TRUE), lin = NULL, fam = gaussian()),
  D = list(nl = bf(y ~ b + a * w, a ~ gp(x), b ~ 1, nl = TRUE), lin = NULL, fam = gaussian()))
rows <- list(list("A", "sin", 60, 2), list("B", "exp", 120, 3), list("B", "step", 60, 1),
             list("B", "wig", 120, 3), list("C", "wig", 60, 2), list("C", "sin", 60, 2),
             list("C", "wig", 120, 2))
for (rw in rows) {
  b <- rw[[1]]; ci <- match(rw[[2]], curves)
  d <- mk(rw[[4]], rw[[3]], ci, b)
  sp <- spell[[b]]
  nl <- fe_fit(sp$nl, data = d, family = sp$fam)
  lin <- if (!is.null(sp$lin)) fe_fit(sp$lin, data = d, family = sp$fam)
  fr <- nl$fit$frame
  tpl <- fr[["par_template"]]
  w <- tryCatch(ns$gp_start_working(fr, tpl), error = function(e) e)
  cat(sprintf("DG %s %s %s n %d seed %d | nl %.4f theta %s | lin %s theta %s | working %s\n",
              arm, b, rw[[2]], rw[[3]], rw[[4]], logLik(nl$fit),
              paste(signif(nl$fit$estimates$theta, 4), collapse = " "),
              if (is.null(lin)) "-" else sprintf("%.4f", logLik(lin$fit)),
              if (is.null(lin)) "-" else paste(signif(lin$fit$estimates$theta, 4), collapse = " "),
              if (inherits(w, "error")) conditionMessage(w) else
                if (is.null(w)) "NULL" else paste(signif(w$theta, 4), collapse = " ")))
  if (!is.null(lin)) {
    l2 <- fe_fit(sp$nl, data = d, family = sp$fam,
                 start = list(theta = lin$fit$estimates$theta))
    cat(sprintf("DG    body from lin theta %.4f\n", logLik(l2$fit)))
  }
}
cat("DONE\n")
