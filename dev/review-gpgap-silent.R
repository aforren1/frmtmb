# Review of lane gpgap, item 3: where the guard says nothing, is the fit the
# same as brms's? Reference: the same tape non-centered at brms's own line
# (absolute 1e-12), which the closed form confirms on gaussian rows. Designs:
# high-information counts, known se(), a nonlinear body, random effects on
# sigma, a mixture. Also ordinary noise on plain trends, by seed. Seed 5601.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-silent.R silent.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
orig <- get("gp_corr", envir = ns)
put <- function(val) {
  unlockBinding("gp_corr", ns); assign("gp_corr", val, envir = ns)
  lockBinding("gp_corr", ns)
}
abs_corr <- function(theta, blk) {
  Q <- blk[["aux_D2"]][[1]] / (2 * exp(2 * theta[2]))
  exp(-Q) + diag(1, nrow(Q)) * (1e-12 * exp(-2 * theta[1]))
}
gp_sd <- function(fit, th = fit$estimates$theta) {
  bk <- Filter(function(k) k$covstruct %in% c("gp", "hsgp"), fit$frame$re_blocks)[[1]]
  exp(th[bk$theta_idx][[1]])
}
cmp <- function(tag, expr) {
  warned <- character()
  fr <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      cat(sprintf("SIL %-30s ERROR %s\n", tag, substr(conditionMessage(e), 1, 100)))
      NULL
    }),
    warning = function(w) {
      warned <<- c(warned, class(w)[1])
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  if (is.null(fr)) return(invisible())
  frn <- fr$frame
  idx <- which(vapply(frn$re_blocks, function(b) identical(b$covstruct, "gp"), NA))
  frn$ncp_blocks <- idx
  put(abs_corr)
  r <- tryCatch({
    fn <- frmtmb::build_objective(frn)
    tpl <- frn$par_template
    tpl$b[] <- 0
    o <- RTMB::MakeADFun(fn, tpl, random = intersect("b", names(tpl)), map = frn$map,
                         silent = TRUE)
    op <- suppressWarnings(nlminb(o$par, o$fn, o$gr,
                                  control = list(iter.max = 800, eval.max = 1200)))
    list(ll = -op$objective, sd = exp(op$par[names(op$par) == "theta"][[1]]),
         conv = op$convergence, gr = max(abs(o$gr(op$par))))
  }, error = function(e) NULL)
  put(orig)
  cat(sprintf(paste0("SIL %-30s frm: ll %.3f sdgp %.4g conv %d nugwarn %-5s | ",
                     "ncp brms line: ll %s sdgp %s conv %s gr %s | ll gap %s\n"),
              tag, -fr$opt$objective, gp_sd(fr), fr$opt$convergence,
              "frmtmb_gp_nugget" %in% warned,
              if (is.null(r)) "NA" else sprintf("%.3f", r$ll),
              if (is.null(r)) "NA" else sprintf("%.4g", r$sd),
              if (is.null(r)) "NA" else r$conv,
              if (is.null(r)) "NA" else sprintf("%.1e", r$gr),
              if (is.null(r)) "NA" else sprintf("%.3f", r$ll + fr$opt$objective)))
}
set.seed(5601)
n <- 60
x <- round(runif(n, 0, 6), 3)
e <- rnorm(n)
d <- data.frame(x = x, z = round(runif(n, 0, 6), 3),
                g = factor(rep(1:6, length.out = n)))
for (lm in c(4, 8, 12, 16)) {
  d$cnt <- rpois(n, exp(lm + sin(2 * d$x)))
  cmp(sprintf("poisson log mean %d", lm),
      frm(bf(cnt ~ gp(x)), family = poisson(), data = d))
}
for (tr in c(1e3, 1e5, 1e7)) {
  d$tr <- tr
  d$suc <- rbinom(n, d$tr, plogis(sin(2 * d$x)))
  cmp(sprintf("binomial trials %g", tr),
      frm(bf(suc | trials(tr) ~ gp(x)), family = binomial(), data = d))
}
for (se in c(1e-2, 3e-3, 1e-3, 3e-4)) {
  d$se <- se
  d$ys <- sin(2 * d$x) + se * e
  cmp(sprintf("se(%g) known", se), frm(bf(ys | se(se) ~ gp(x)), data = d))
}
d$y <- sin(2 * d$x) + 1e-3 * e
d$yy <- d$y + 0.1 * d$z
cmp("nl b1 ~ gp(x), noise 1e-3", frm(bf(yy ~ b1 + b2 * z, b1 ~ gp(x), b2 ~ 1, nl = TRUE),
                                     data = d))
cmp("sigma ~ (1|g), noise 1e-3", frm(bf(y ~ gp(x), sigma ~ (1 | g)), data = d))
cmp("sigma ~ s(z), noise 1e-3", frm(bf(y ~ gp(x), sigma ~ s(z, k = 5)), data = d))
cmp("mixture, noise 1e-3", frm(bf(y ~ 1, mu1 ~ gp(x), mu2 ~ 1),
                               family = mixture(gaussian(), gaussian()), data = d))
cmp("student, noise 1e-3", frm(bf(y ~ gp(x)), family = student(), data = d))
cmp("gaussian, noise 1e-3 (control)", frm(bf(y ~ gp(x)), data = d))

# ordinary noise, a plain trend, seeds. Closed form at each nugget.
fitm <- function(x, y, nugget) {
  u <- sort(unique(x)); D2 <- outer(u, u, "-")^2; J <- match(x, u)
  nll <- function(p) {
    s2 <- exp(2 * p[2])
    K <- s2 * exp(-D2 / (2 * exp(2 * p[3]))); diag(K) <- diag(K) + nugget(s2)
    V <- K[J, J]; diag(V) <- diag(V) + exp(2 * p[4])
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) return(1e10)
    r <- backsolve(R, y - p[1], transpose = TRUE)
    0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * length(y) * log(2 * pi)
  }
  best <- NULL
  for (st in list(c(mean(y), log(sd(y)), 0, -1),
                  c(mean(y), log(sd(y)) + 2, 1.5, -1),
                  c(mean(y), log(sd(y)) + 4, 2.5, -1),
                  c(mean(y), log(sd(y)) + 1, 0.7, -1))) {
    o <- nlminb(st, nll, control = list(iter.max = 800, eval.max = 1200))
    if (is.null(best) || o$objective < best$objective) best <- o
  }
  best
}
for (seed in c(5701L, 5702L)) {
  set.seed(seed)
  x <- round(runif(100, 0, 10), 3)
  e <- rnorm(100)
  for (lab in c("0.2x", "x", "5x", "0.1x^2", "x^2", "exp(0.5x)", "exp(0.8x)", "logis")) {
    f <- switch(lab, "0.2x" = 0.2 * x, x = x, "5x" = 5 * x, "0.1x^2" = 0.1 * x^2,
                "x^2" = x^2, "exp(0.5x)" = exp(0.5 * x), "exp(0.8x)" = exp(0.8 * x),
                logis = 100 * plogis(x - 5))
    y <- f + 0.3 * e
    warned <- FALSE
    fr <- withCallingHandlers(
      frm(bf(y ~ gp(x)), data = data.frame(x, y)),
      warning = function(w) {
        if (inherits(w, "frmtmb_gp_nugget")) warned <<- TRUE
        invokeRestart("muffleWarning")
      })
    a <- fitm(x, y, function(s2) 1e-6 * s2)
    b <- fitm(x, y, function(s2) 1e-12)
    cat(sprintf(paste0("TRD seed %d %-10s frm sigma %.4f sdgp %.4g warn %-5s | closed frm nugget ",
                       "sigma %.4f sdgp %.4g ll %.3f | closed brms sigma %.4f sdgp %.4g rho %.3g ",
                       "ll %.3f | sigma ratio %.3f gap %.3f\n"),
                seed, lab, as.numeric(sigma(fr)), gp_sd(fr), warned,
                exp(a$par[4]), exp(a$par[2]), -a$objective,
                exp(b$par[4]), exp(b$par[2]), exp(b$par[3]), -b$objective,
                exp(a$par[4] - b$par[4]), a$objective - b$objective))
  }
}
