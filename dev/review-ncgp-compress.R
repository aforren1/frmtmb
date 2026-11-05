# Item 3: on real frames, the outer gradient on the compressed Hessian tape (what
# frm() uses) against the plain tape and a central difference, over many models,
# seeds and points, b first as shipped. A wrong gradient here is silent.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
specs <- list(
  gp = list(f = y ~ gp(x)), gp_reml = list(f = y ~ z + gp(x), REML = TRUE),
  re_gp = list(f = y ~ (1 | g) + gp(x)), gp_re = list(f = y ~ gp(x) + (1 | g)),
  re_gp_reml = list(f = y ~ (1 | g) + gp(x), REML = TRUE),
  slope = list(f = y ~ (z | g) + gp(x)),
  two = list(f = y ~ gp(x) + gp(z)), two_re = list(f = y ~ gp(z) + (1 | g) + gp(x)),
  s_gp = list(f = y ~ s(z) + gp(x)), by = list(f = y ~ gp(x, by = f)),
  aniso = list(f = y ~ gp(x, z, iso = FALSE)), sigma = list(f = bf(y ~ gp(x), sigma ~ gp(z))),
  mv = list(f = bf(y ~ gp(x)) + bf(y2 ~ gp(z)) + set_rescor(FALSE)),
  pois = list(f = cnt ~ (1 | g) + gp(x), family = poisson()),
  bern = list(f = bin ~ gp(x), family = bernoulli()),
  prof = list(f = y ~ z + gp(x), control = frmtmb_control(profile = TRUE)),
  mi = list(f = bf(y ~ mi(zm) + gp(x)) + bf(zm | mi() ~ 1) + set_rescor(FALSE)),
  grid = list(f = y ~ gp(xg)), rep = list(f = y ~ gp(xr))
)
tol <- list(); bad <- 0L; tot <- 0L
for (seed in 1:4) {
  set.seed(600 + seed)
  n <- 70
  d <- data.frame(x = runif(n, 0, 6), z = runif(n, 0, 6), g = factor(sample(letters[1:6], n, TRUE)),
                  f = factor(sample(c("u", "v"), n, TRUE)), xg = seq(0, 6, length.out = n),
                  xr = rep(sort(runif(14, 0, 6)), each = 5))
  sg <- c(0.3, 0.01, 1e-3, 0.1)[seed]
  d$y <- sin(2 * d$x) + 0.3 * d$z + rnorm(6)[as.integer(d$g)] * 0.5 + sin(d$xg) + sin(d$xr) + rnorm(n, sd = sg)
  d$y2 <- cos(d$z) + rnorm(n, sd = sg)
  d$cnt <- rpois(n, exp(1 + sin(d$x))); d$bin <- rbinom(n, 1, plogis(sin(d$x)))
  d$zm <- d$z; d$zm[sample(n, 6)] <- NA
  for (nm in names(specs)) {
    s <- specs[[nm]]; a <- s; a$f <- NULL
    fit <- tryCatch(suppressWarnings(suppressMessages(do.call(frm, c(list(s$f, data = d), a)))), error = function(e) e)
    if (inherits(fit, "error")) { cat("RC", nm, seed, "FIT ERROR", substr(conditionMessage(fit), 1, 100), "\n"); next }
    fr <- fit$frame
    rnd <- unique(names(fit$obj$env$par)[fit$obj$env$random])
    tpl <- ns$gp_tape_order(fr, fit$estimates[names(fr$par_template)])
    nll <- ns$build_objective(fr)
    plain <- tryCatch(RTMB::MakeADFun(nll, tpl, random = rnd, map = fr$map, silent = TRUE), error = function(e) e)
    comp <- tryCatch(ns$gp_tape_scope(fr, RTMB::MakeADFun(nll, tpl, random = rnd, map = fr$map, silent = TRUE)), error = function(e) e)
    if (inherits(plain, "error") || inherits(comp, "error")) { cat("RC", nm, seed, "TAPE ERROR\n"); next }
    p0 <- plain$par
    for (ipt in 1:3) {
      set.seed(ipt + seed * 10)
      p <- p0 + if (ipt == 1L) 0 else rnorm(length(p0), sd = 0.15)
      g1 <- tryCatch(plain$gr(p), error = function(e) NA); g2 <- tryCatch(comp$gr(p), error = function(e) NA)
      f1 <- plain$fn(p); f2 <- comp$fn(p)
      fd <- vapply(seq_along(p), function(i) { h <- 1e-4; e <- replace(numeric(length(p)), i, h); (comp$fn(p + e) - comp$fn(p - e)) / (2 * h) }, 0)
      md <- max(abs(g1 - g2), na.rm = TRUE); rd <- md / max(1, max(abs(g1), na.rm = TRUE))
      fdd <- max(abs(g2 - fd) / pmax(1, abs(fd)), na.rm = TRUE)
      tot <- tot + 1L
      flag <- if (!is.finite(rd) || rd > 1e-4) "BAD" else ""
      if (nzchar(flag)) bad <- bad + 1L
      cat(sprintf("RC %-10s seed %d pt %d fn diff %.1e | gr plain-vs-compressed %.2e | compressed-vs-fd %.2e %s\n", nm, seed, ipt, f1 - f2, rd, fdd, flag))
    }
    flush.console()
  }
}
cat(sprintf("RCSUM %d readings, %d differ\n", tot, bad))
