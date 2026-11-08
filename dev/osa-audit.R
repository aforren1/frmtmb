# Lane osa: every core family against residuals(type = "osa").
#
#   Rscript dev/osa-audit.R before|after
#
# For each family of the core registry: fit one fixed design with no
# random effect (data seed 9, n = 150), record the oneStepPredict() route
# residuals() takes (traced from the call itself), and compare the
# residuals with the fitted model's own PIT.
#
# Truth. With no random effect the one-step predictive of a row is the
# fitted row distribution, so the residual is the (randomized) quantile
# residual of the fitted CDF. The CDF is formed here from the family's
# DENSITY alone, on doubles, never from its cdf_log slot or TMB: a sum of
# the mass over the support for a discrete family, stats::integrate() of
# the density for a continuous one (normalized by the integral over the
# whole support, as TMB normalizes). For a discrete family the residual
# must lie in [qnorm(F(y-1)), qnorm(F(y))] ("outside" counts the rows that
# do not, with a slack of 1e-6 of the row's own interval), and with TMB's
# own uniforms (set.seed(123); runif(n)) it is qnorm(F(y) - U p(y))
# ("maxdiff"). For a continuous family "maxdiff" is |r - qnorm(F(y))|.
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "after"
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm:", arm, " frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), " RTMBdist",
    format(packageVersion("RTMBdist")), " RTMB",
    format(packageVersion("RTMB")), "\n")

route <- new.env()
trace("oneStepPredict", where = asNamespace("RTMB"), print = FALSE,
      tracer = quote({
        a <- list(...)
        assign("method", a$method %||% "(default)", envir = route)
        assign("discrete", if (is.null(a$discrete)) NA else a$discrete,
               envir = route)
        ds <- a$discreteSupport
        assign("support", if (is.null(ds)) "-" else
                 sprintf("%g:%g", min(ds), max(ds)), envir = route)
        rg <- a$range
        assign("range", if (is.null(rg)) "-" else
                 sprintf("[%g,%g]", rg[1], rg[2]), envir = route)
      }))
`%||%` <- function(a, b) if (is.null(a)) b else a

set.seed(9)
n <- 150
x <- rnorm(n)
size <- 10L
eta <- -0.3 + 0.8 * x
pb <- plogis(eta)
mcount <- exp(0.8 + 0.5 * x)
zi <- runif(n) < 0.25
ys <- list(
  real = 1 + 0.7 * x + rnorm(n),
  pos = rgamma(n, shape = 3, rate = 3 / exp(0.3 + 0.4 * x)),
  unit = rbeta(n, pb * 6, (1 - pb) * 6),
  count = rnbinom(n, size = 3, mu = mcount),
  trials = rbinom(n, size, pb),
  binary = rbinom(n, 1, pb),
  ord = as.integer(cut(eta + rlogis(n), c(-Inf, -1, 0, 1, Inf))),
  zcount = ifelse(zi, 0, rpois(n, mcount)),
  hcount = ifelse(zi, 0, 1 + rpois(n, mcount)),
  ztrials = ifelse(zi, 0, rbinom(n, size, pb)),
  zpos = ifelse(zi, 0, rgamma(n, 3, 3)),
  zunit = ifelse(zi, 0, rbeta(n, 2, 3)),
  circ = atan2(sin(0.5 * x + rnorm(n, 0, 0.5)), cos(0.5 * x +
                                                       rnorm(n, 0, 0.5))))
# support of each design: where the truth's sum or integral starts and ends
supp <- list(real = c(-Inf, Inf), pos = c(0, Inf), unit = c(0, 1),
             count = c(0, Inf), trials = c(0, size), binary = c(0, 1),
             ord = c(1, 4), zcount = c(0, Inf), hcount = c(0, Inf),
             ztrials = c(0, size), zpos = c(0, Inf), zunit = c(0, 1),
             circ = c(-pi, pi))
design <- c(
  gaussian = "real", student = "real", skew_normal = "real",
  exgaussian = "real", asym_laplace = "real", huber = "real",
  gen_extreme_value = "real",
  Gamma = "pos", lognormal = "pos", inverse.gaussian = "pos",
  exponential = "pos", weibull = "pos", shifted_lognormal = "pos",
  frechet = "pos", tweedie = "zpos",
  beta = "unit", xbeta = "unit", zero_inflated_beta = "zunit",
  zero_one_inflated_beta = "zunit",
  poisson = "count", negbinomial = "count", negbinomial2 = "count",
  nbinom1 = "count", compois = "count", geometric = "count",
  discrete_weibull = "count",
  binomial = "trials", beta_binomial = "trials", bernoulli = "binary",
  zero_inflated_poisson = "zcount", zero_inflated_negbinomial = "zcount",
  hurdle_poisson = "hcount", hurdle_negbinomial = "hcount",
  zero_inflated_binomial = "ztrials",
  zero_inflated_beta_binomial = "ztrials",
  hurdle_gamma = "zpos", hurdle_lognormal = "zpos",
  zero_inflated_asym_laplace = "real",
  cumulative = "ord", sratio = "ord", cratio = "ord", acat = "ord",
  hurdle_cumulative = "ord", von_mises = "circ",
  categorical = "ord", multinomial = "-", dirichlet_multinomial = "-",
  dirichlet = "-", dirichlet2 = "-", logistic_normal = "-")

U <- local({ set.seed(123); runif(n) })
reg <- frmtmb:::family_registry

# P(Y <= y) and p(y) of row i from the family's density on doubles
row_truth <- function(fam, dp, av, i, y, sp, disc) {
  sub <- function(v) if (length(v) == n) v[i] else v
  dpi <- lapply(dp, sub)
  avi <- lapply(av, sub)
  f <- function(t) exp(as.numeric(fam$lpdf(t, dpi, avi)))
  if (disc) {
    top <- if (is.finite(sp[2])) sp[2] else max(y + 2000, 5000)
    lo <- if (y > sp[1]) sum(f(seq(sp[1], y - 1))) else 0
    p <- f(y)
    return(c(lo = lo, hi = lo + p))
  }
  g <- function(t) vapply(t, f, 0)
  i1 <- stats::integrate(g, sp[1], y, rel.tol = 1e-10,
                         subdivisions = 1000L)$value
  i2 <- stats::integrate(g, y, sp[2], rel.tol = 1e-10,
                         subdivisions = 1000L)$value
  c(lo = i1 / (i1 + i2), hi = i1 / (i1 + i2))
}

rows <- list()
for (nm in names(design)) {
  dz <- design[[nm]]
  rec <- data.frame(family = nm, method = NA, discrete = NA, support = NA,
                    range = NA, nan = NA, outside = NA, maxdiff = NA,
                    secs = NA, note = "")
  if (dz == "-") {
    rec$note <- "matrix/simplex response: not probed"
    rows[[nm]] <- rec
    next
  }
  d <- data.frame(y = ys[[dz]], x = x, k = size)
  form <- if (dz %in% c("trials", "ztrials")) bf(y | trials(k) ~ x) else
    bf(y ~ x)
  famobj <- tryCatch(reg[[nm]](), error = function(e) NULL)
  if (nm == "categorical") d$y <- factor(d$y)
  if (dz == "ord" && nm != "categorical") d$y <- factor(d$y, ordered = TRUE)
  fit <- tryCatch(suppressWarnings(frm(form, family = famobj, data = d)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    rec$note <- paste("FIT:", substr(conditionMessage(fit), 1, 70))
    rows[[nm]] <- rec
    next
  }
  rm(list = ls(route), envir = route)
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(suppressWarnings(residuals(fit, type = "osa")[, "Estimate"]),
                error = function(e) e)
  rec$secs <- round(proc.time()[["elapsed"]] - t0, 2)
  for (k in c("method", "discrete", "support", "range")) {
    if (exists(k, envir = route, inherits = FALSE)) {
      rec[[k]] <- get(k, envir = route, inherits = FALSE)
    }
  }
  if (inherits(r, "error")) {
    rec$note <- paste("REFUSED:", substr(gsub("\n", " ",
                                             conditionMessage(r)), 1, 90))
    rows[[nm]] <- rec
    next
  }
  fin <- is.finite(r)
  rec$nan <- sum(!fin)
  fam <- fit$spec$responses[[1]]$family
  disc <- fam$type %in% c("discrete", "ordinal")
  dp <- lapply(frmtmb:::eval_dpars(fit)[[1]], function(v) as.numeric(v))
  av <- fit$frame$aterm_values[[1]]
  yv <- as.numeric(fit$frame$y[[1]])
  tr <- if (identical(fam$type, "ordinal")) {
    # an ordinal density reads thresholds the objective passes as extra
    # parameters; the fitted category probabilities are the same model
    P <- frmtmb:::ord_probs(fit, fit$spec$responses[[1]])
    cP <- t(apply(P, 1, cumsum))
    cbind(lo = ifelse(yv > 1, cP[cbind(seq_len(n), pmax(yv - 1, 1))], 0),
          hi = cP[cbind(seq_len(n), yv)])
  } else {
    tryCatch(t(vapply(seq_len(n), function(i) {
      row_truth(fam, dp, av, i, yv[i], supp[[dz]], disc)
    }, c(lo = 0, hi = 0))), error = function(e) e)
  }
  if (!inherits(tr, "error")) tr <- pmin(tr, 1)
  if (inherits(tr, "error")) {
    rec$note <- paste("truth:", substr(conditionMessage(tr), 1, 60))
    rows[[nm]] <- rec
    next
  }
  if (disc) {
    qlo <- qnorm(tr[, "lo"])
    qhi <- qnorm(tr[, "hi"])
    wid <- pmax(qhi - qlo, 0)
    wid[!is.finite(wid)] <- 1
    rec$outside <- sum(fin & (r < qlo - 1e-6 * wid | r > qhi + 1e-6 * wid))
    ex <- qnorm(tr[, "hi"] - U * (tr[, "hi"] - tr[, "lo"]))
  } else {
    ex <- qnorm(tr[, "hi"])
  }
  ok <- fin & is.finite(ex)
  rec$maxdiff <- if (any(ok)) signif(max(abs(r - ex)[ok]), 3) else NA
  rows[[nm]] <- rec
}
res <- do.call(rbind, rows)
rownames(res) <- NULL
cat("BEGIN osa-audit", arm, "\n")
op <- options(width = 250)
print(res, right = FALSE)
options(op)
cat("END osa-audit\n")
