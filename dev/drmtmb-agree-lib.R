# Shared data generators and the comparison engine for the drmTMB
# agreement work. Sourced by dev/drmtmb-agree.R and by exploratory
# scripts, so every number in dev/drmtmb-findings.md comes from one
# construction.
.libPaths(c("C:/Users/adf44/source/r/drmtmb-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)          # attached: frm(), bf(), gr(), families
  requireNamespace("drmTMB")
})
D <- asNamespace("drmTMB")
dbf <- drmTMB::bf

# Grouped design shared by the univariate models. Group-level
# covariate w is constant within group so sd(g) ~ w is admissible.
sim_grouped <- function(seed, ng = 40, nper = 10) {
  set.seed(seed)
  n <- ng * nper
  g <- factor(rep(seq_len(ng), each = nper))
  u <- rnorm(ng, 0, 0.6)[g]
  v <- rnorm(ng, 0, 0.3)[g]
  x <- rnorm(n); z <- rnorm(n)
  d <- data.frame(g = g, x = x, z = z, w = rnorm(ng)[g])
  d$y <- 1 + 0.5 * x + u + rnorm(n, 0, exp(-0.2 + 0.3 * z + v))
  m <- plogis(0.2 + 0.4 * x + u)
  d$yb <- rbeta(n, m * 20, (1 - m) * 20)
  d$yc <- rnbinom(n, mu = exp(1 + 0.3 * x + u), size = 3)
  d$yt <- 1 + 0.5 * x + u + 0.8 * rt(n, df = 5)
  lat <- 0.8 * x + u + rlogis(n)
  d$yo <- factor(cut(lat, c(-Inf, -1, 0.5, 2, Inf), labels = FALSE),
                 ordered = TRUE)
  d$ysn <- 1 + 0.5 * x + u + 0.8 * (abs(rnorm(n)) - sqrt(2 / pi))
  # y2's group effect is correlated with y's but not perfectly, so the
  # group correlation is interior. y2b shares y's group effect exactly,
  # which puts that correlation on its boundary on purpose.
  e1 <- d$y - 1 - 0.5 * x - u
  u2 <- 0.5 * u + rnorm(ng, 0, 0.4)[g]
  d$y2 <- 0.3 + 0.2 * x + u2 + 0.5 * e1 + rnorm(n, 0, 0.8)
  d$y2b <- 0.3 + 0.2 * x + 0.7 * u + 0.5 * e1 + rnorm(n, 0, 0.8)
  d
}

# Additive relationship matrix from a pedigree by the tabular method.
# Written here rather than taken from either package so that neither
# fit is compared against its own builder.
additive_A <- function(ped) {
  n <- nrow(ped)
  A <- matrix(0, n, n, dimnames = list(ped$id, ped$id))
  idx <- function(p) if (is.na(p)) NA_integer_ else match(p, ped$id)
  for (i in seq_len(n)) {
    s <- idx(ped$sire[i]); dm <- idx(ped$dam[i])
    for (j in seq_len(i - 1)) {
      a <- 0
      if (!is.na(s)) a <- a + A[j, s]
      if (!is.na(dm)) a <- a + A[j, dm]
      A[i, j] <- A[j, i] <- a / 2
    }
    A[i, i] <- 1 + if (!is.na(s) && !is.na(dm)) A[s, dm] / 2 else 0
  }
  A
}

sim_pedigree <- function(seed, nfound = 20, nper_gen = 40, ngen = 3,
                         nrec = 3) {
  set.seed(seed)
  ped <- data.frame(id = paste0("f", seq_len(nfound)), sire = NA_character_,
                    dam = NA_character_)
  parents <- ped$id
  for (k in seq_len(ngen)) {
    ids <- paste0("g", k, "_", seq_len(nper_gen))
    half <- length(parents) %/% 2
    sires <- sample(parents[seq_len(half)], nper_gen, TRUE)
    dams <- sample(parents[-seq_len(half)], nper_gen, TRUE)
    ped <- rbind(ped, data.frame(id = ids, sire = sires, dam = dams))
    parents <- ids
  }
  A <- additive_A(ped)
  a <- as.vector(t(chol(A)) %*% rnorm(nrow(A))) * 0.7
  d <- data.frame(id = factor(rep(ped$id, each = nrec), levels = ped$id),
                  x = rnorm(nrow(ped) * nrec))
  d$y <- 1 + 0.4 * d$x + a[as.integer(d$id)] + rnorm(nrow(d), 0, 0.8)
  list(data = d, A = A, ped = ped)
}

sim_meta <- function(seed, k = 60) {
  set.seed(seed)
  vi <- runif(k, 0.02, 0.3)
  x <- rnorm(k)
  data.frame(study = factor(seq_len(k)), x = x, vi = vi, sei = sqrt(vi),
             yi = 0.3 + 0.2 * x + rnorm(k, 0, 0.25) + rnorm(k, 0, sqrt(vi)))
}

# Evaluate a fitted Laplace/REML objective at an arbitrary outer point.
# Both packages expose a TMB ADFun whose fn() is the negative
# log-likelihood of the outer parameters.
negll_at <- function(obj, par) as.numeric(obj$fn(unname(par)))

# A map row states one parameter pair: f index into frmtmb's outer
# vector, d index into drmTMB's, fun (drm -> frm), inv (frm -> drm), dfun
# (derivative of fun, for mapping a standard error by the delta method)
# and a label. Every pair is stated so the mapping is auditable.
mrow <- function(f, d, fun, inv, dfun, label) {
  list(f = f, d = d, fun = fun, inv = inv, dfun = dfun, label = label)
}
lin <- function(f, d, a = 1, b = 0) {
  mrow(f, d, function(x) a * x + b, function(y) (y - b) / a,
       function(x) a + 0 * x,
       if (a == 1 && b == 0) "identity" else sprintf("f = %g*d%+g", a, b))
}
# drmTMB: nu = 2 + exp(d). brms and frmtmb: nu = 1 + exp(f).
nu_map <- function(f, d) {
  mrow(f, d, function(x) log1p(exp(x)), function(y) log(expm1(y)),
       function(x) exp(x) / (1 + exp(x)),
       "nu: 2+exp(d) = 1+exp(f)")
}
# drmTMB: rho = cc * tanh(d). frmtmb: rho = f / sqrt(1 + f^2). drmTMB
# caps random-effect correlations at cc = 0.999999, measured in
# dev/drmtmb-boundary.R; rho12 uses cc = 1, measured off the optimum.
rho_map <- function(f, d, cc = 1) {
  fun <- function(x) { r <- cc * tanh(x); r / sqrt(1 - r^2) }
  inv <- function(y) atanh(y / sqrt(1 + y^2) / cc)
  dfun <- function(x) {
    r <- cc * tanh(x); cc * (1 - tanh(x)^2) / (1 - r^2)^1.5
  }
  mrow(f, d, fun, inv, dfun,
       sprintf("rho: %s*tanh(d) = f/sqrt(1+f^2)", format(cc, digits = 7)))
}
rho_re_map <- function(f, d) rho_map(f, d, cc = 0.999999)
# NB2 sigma random-intercept SD: log shape = -2 log sigma, so the SD of
# the shape-scale effect is twice that of the sigma-scale effect.
log2x <- function(f, d) lin(f, d, 1, log(2))

compare_fits <- function(label, fd, ff, map, extra = NULL) {
  pd <- fd$opt$par
  pf <- ff$opt$par
  fi <- as.integer(vapply(map, `[[`, 0, "f"))
  di <- as.integer(vapply(map, `[[`, 0, "d"))
  stopifnot(length(map) == length(pd), length(map) == length(pf),
            setequal(fi, seq_along(pf)), setequal(di, seq_along(pd)))
  d2f <- pf; f2d <- pd
  for (m in map) {
    d2f[m$f] <- m$fun(pd[m$d]); f2d[m$d] <- m$inv(pf[m$f])
  }
  Vf <- vcov(ff, full = TRUE)
  stopifnot(nrow(Vf) == length(pf))
  sef <- sqrt(diag(Vf)); sed <- sqrt(diag(fd$sdr$cov.fixed))
  ll_d <- as.numeric(logLik(fd)); ll_f <- as.numeric(logLik(ff))
  Ff_at_d <- -negll_at(ff$obj, d2f); Fd_at_d <- -negll_at(fd$obj, pd)
  Fd_at_f <- -negll_at(fd$obj, f2d); Ff_at_f <- -negll_at(ff$obj, pf)
  # The optima test the maps only to second order, because the gradient
  # vanishes there. A point one standard error away tests them to first
  # order. The sign pattern is fixed rather than drawn, so that this
  # function leaves the RNG stream alone; the test file uses the same
  # pattern.
  pd_off <- pd + sed * rep(c(1, -1), length.out = length(pd))
  pf_off <- pf
  for (m in map) pf_off[m$f] <- m$fun(pd_off[m$d])
  off_gap <- (-negll_at(ff$obj, pf_off)) - (-negll_at(fd$obj, pd_off))
  est_d <- vapply(map, function(m) m$fun(pd[m$d]), 0)
  se_d <- vapply(map, function(m) abs(m$dfun(pd[m$d])) * sed[m$d], 0)
  est_z <- (pf[fi] - est_d) / sef[fi]
  se_ratio <- sef[fi] / se_d
  tab <- data.frame(
    par_frm = paste0(names(pf)[fi], "[", fi, "]"),
    par_drm = paste0(names(pd)[di], "[", di, "]"),
    map = vapply(map, `[[`, "", "label"),
    est_frm = pf[fi], est_drm_mapped = est_d, diff_over_se = est_z,
    se_frm = sef[fi], se_drm_mapped = se_d, se_ratio = se_ratio,
    row.names = NULL)
  summ <- data.frame(
    model = label, ll_frm = ll_f, ll_drm = ll_d, ll_diff = ll_f - ll_d,
    # Same point, two objectives: zero means the same likelihood.
    obj_gap_at_drm_opt = Ff_at_d - Fd_at_d,
    obj_gap_at_frm_opt = Ff_at_f - Fd_at_f,
    obj_gap_off_opt = off_gap,
    max_abs_diff_over_se = max(abs(est_z)),
    max_abs_log_se_ratio = max(abs(log(se_ratio))),
    grad_frm = max(abs(ff$obj$gr(pf))), grad_drm = max(abs(fd$obj$gr(pd))),
    conv_frm = ff$opt$convergence, conv_drm = fd$opt$convergence)
  list(summary = summ, table = tab, extra = extra)
}

fmt_block <- function(res) {
  op <- options(width = 200, digits = 10)
  on.exit(options(op))
  cat("\n==========", res$summary$model, "\n")
  print(res$summary, row.names = FALSE)
  print(res$table, row.names = FALSE)
  if (!is.null(res$extra)) print(res$extra)
}
