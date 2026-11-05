# The tape's outer gradient against a central difference of its value, per
# spelling, at a point off the optimum so the gradient is not near zero, under
# the compressed Hessian and without it. A wrong reverse sweep shows here and
# nowhere else. Seed 20261106.
#   Rscript dev/ncgp-grcheck.R
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261106)
n <- 40L
d <- data.frame(x = round(runif(n, 0, 6), 3), z = round(runif(n, 0, 4), 2),
                f = factor(rep(c("a", "b", "c"), length.out = n)),
                w = runif(n, 0.5, 2), g = factor(rep(1:5, length.out = n)))
d$y <- sin(2 * d$x) + 0.1 * rnorm(n)
d$y2 <- cos(d$z) + 0.1 * rnorm(n)
d$cnt <- rpois(n, exp(1 + sin(2 * d$x)))
d$bin <- rbinom(n, 1, plogis(sin(2 * d$x)))
d$ysg <- rnorm(n, 0, exp(-1 + sin(d$x)))
d$ynl <- (2 + sin(d$z)) * exp(0.3 * d$w) + 0.2 * rnorm(n)
d$zm <- d$z; d$zm[c(3, 9, 20)] <- NA
d$g2 <- factor(rep(1:2, length.out = n))
sp <- list(
  list("gp(x)", bf(y ~ gp(x)), gaussian()),
  list("gp(x) REML", bf(y ~ gp(x)), gaussian(), REML = TRUE),
  list("gp(x)+z REML", bf(y ~ z + gp(x)), gaussian(), REML = TRUE),
  list("gp(x, z)", bf(y ~ gp(x, z)), gaussian()),
  list("aniso", bf(y ~ gp(x, z, iso = FALSE)), gaussian()),
  list("by = f", bf(y ~ gp(x, by = f)), gaussian()),
  list("by = w", bf(y ~ gp(x, by = w)), gaussian()),
  list("two gp", bf(y ~ gp(x) + gp(z)), gaussian()),
  list("gp + (1|g)", bf(y ~ gp(x) + (1 | g)), gaussian()),
  list("gp + (1|g) REML", bf(y ~ gp(x) + (1 | g)), gaussian(), REML = TRUE),
  list("gp + s(z)", bf(y ~ gp(x) + s(z, k = 5)), gaussian()),
  list("gp + hsgp", bf(y ~ gp(x) + gp(z, k = 6)), gaussian()),
  list("sigma~gp", bf(ysg ~ 1, sigma ~ gp(x)), gaussian()),
  list("nl", bf(ynl ~ b1 * exp(b2 * w), b1 ~ gp(z), b2 ~ 1, nl = TRUE), gaussian()),
  list("mvbf", mvbf(bf(y ~ gp(x)), bf(y2 ~ gp(z)), rescor = FALSE), NULL),
  list("poisson", bf(cnt ~ gp(x)), poisson()),
  list("bernoulli", bf(bin ~ gp(x)), bernoulli()),
  list("poisson + (1|g)", bf(cnt ~ gp(x) + (1 | g)), poisson()),
  list("student", bf(y ~ gp(x)), student()),
  list("profile", bf(y ~ z + gp(x)), gaussian(), control = frmtmb_control(profile = TRUE)),
  list("mi + gp", bf(y ~ mi(zm) + gp(x)) + bf(zm | mi() ~ 1), gaussian()),
  list("REML 1-level", bf(y ~ z + gp(x) + (1 | g2)), gaussian(), REML = TRUE),
  list("(1|g2) + gp", bf(y ~ (1 | g2) + gp(x)), gaussian())
)
ng <- function(o, p, h = 1e-5) vapply(seq_along(p), function(i) {
  e <- replace(numeric(length(p)), i, h)
  (o$fn(p + e) - o$fn(p - e)) / (2 * h)
}, 0)
mk <- function(fit, comp) {
  fr <- fit$frame
  old <- TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress
  tc <- RTMB::TapeConfig()
  on.exit({TMB::config(tmbad.sparse_hessian_compress = old, DLL = "RTMB"); RTMB::TapeConfig(tc)})
  TMB::config(tmbad.sparse_hessian_compress = comp, DLL = "RTMB")
  RTMB::TapeConfig(matmul = "atomic")
  rnd <- unique(names(fit$obj$env$par)[fit$obj$env$random])
  o <- RTMB::MakeADFun(frmtmb:::build_objective(fr), frmtmb::gp_tape_order(fr, fit$estimates), random = rnd,
                       map = fr$map, silent = TRUE)
  force(o$env$spHess)
  o
}
for (s in sp) {
  args <- c(list(s[[2]], data = d), if (!is.null(s[[3]])) list(family = s[[3]]), s[-(1:3)])
  fit <- suppressWarnings(do.call(frm, args))
  p <- fit$opt$par
  set.seed(1)
  p1 <- p + 0.1 * rnorm(length(p))
  for (comp in 1:0) {
    o <- mk(fit, comp)
    g <- o$gr(p1)
    fd <- ng(o, p1)
    rel <- max(abs(g - fd) / (abs(fd) + 1e-2 * max(abs(fd))))
    cat(sprintf("NCGPG %-18s compress %d fn %.8f max rel gr err %.2e (worst at %s; gr %.5g fd %.5g)\n",
                s[[1]], comp, o$fn(p1), rel,
                names(p1)[which.max(abs(g - fd))], g[which.max(abs(g - fd))], fd[which.max(abs(g - fd))]))
  }
}
