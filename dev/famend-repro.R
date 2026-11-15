# Lane famend: the family items of the brief, before and after, with the
# guard-absent cases, on one build.
#   Rscript dev/famend-repro.R <lib or "base"> [case ...]
# "base" reads rellib-r21 alone; any other value puts that library first.
a <- commandArgs(trailingOnly = TRUE)
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] != "base") a[1], base, ulib))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
cases <- a[-1]
want <- function(k) !length(cases) || k %in% cases

cap <- function(expr) {
  w <- character(); m <- character()
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
    }),
    error = function(e) structure(conditionMessage(e), class = "err"))
  list(v = v, w = w, m = m)
}
showw <- function(r, n = 400) {
  for (x in r$w) cat("   W:", substr(gsub("\\s+", " ", x), 1, n), "\n")
  for (x in r$m) cat("   M:", substr(gsub("\\s+", " ", x), 1, n), "\n")
  if (inherits(r$v, "err")) cat("   ERR:", substr(r$v, 1, n), "\n")
}
brief <- function(label, r) {
  cat("\n==", label, "\n")
  showw(r)
  f <- r$v
  if (inherits(f, "err")) return(invisible(NULL))
  cat("   code", f$opt$convergence, " logLik",
      format(as.numeric(logLik(f)), digits = 10), "\n")
  se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
  print(data.frame(par = names(f$opt$par), est = signif(f$opt$par, 6),
                   se = signif(se[seq_along(f$opt$par)], 4)),
        row.names = FALSE)
  invisible(f)
}

if (want("frechet")) {
  # 3060: frechet at nu0 = 1.1, seed 31 of dev/famuni-rev-frechet11.R
  nu0 <- 1.1
  set.seed(9000 + 31)
  d <- data.frame(x = rnorm(300))
  d$y <- exp(0.5 + 0.3 * d$x) / gamma(1 - 1 / nu0) * rexp(300)^(-1 / nu0)
  brief("frechet: nu0 = 1.1, seed 9031 (3060)",
        cap(frm(y ~ x, family = frechet(), data = d)))
}
if (want("frechet_absent")) {
  # guard absent: nu0 = 1.1 seeds whose fit has an interior maximum, and
  # nu0 = 1.5 and 4
  for (s in c(1, 2, 3)) {
    nu0 <- 1.1
    set.seed(9000 + s)
    d <- data.frame(x = rnorm(300))
    d$y <- exp(0.5 + 0.3 * d$x) / gamma(1 - 1 / nu0) * rexp(300)^(-1 / nu0)
    brief(paste0("frechet absent: nu0 = 1.1, seed ", 9000 + s),
          cap(frm(y ~ x, family = frechet(), data = d)))
  }
}
if (want("gev")) {
  # 3075: xi in (-1, -0.5), and -0.3 (guard absent)
  rg <- function(n, mu, s, xi) mu + s * ((-log(runif(n)))^(-xi) - 1) / xi
  for (xi in c(-0.63, -0.7, -0.89, -0.3, -0.45)) {
    set.seed(4)
    dv <- data.frame(y = rg(300, 0, 1, xi))
    brief(paste0("gev: xi0 = ", xi, ", seed 4 (3075)"),
          cap(frm(y ~ 1, family = gen_extreme_value(), data = dv)))
  }
}
if (want("gevgrad")) {
  # 3067: an outside GEV row, xi 0.06, sigma 1: x where the gradient of
  # the density is not finite (dev/famuni-rev-numerics2.R, part B)
  famg <- gen_extreme_value()
  for (xi0 in c(0.06, 0.3, 1)) {
    xs <- seq(-1 / xi0 - 0.01, -1 / xi0 - 60, by = -0.01)
    p <- c(0, 0, log1p(xi0))
    bad <- logical(length(xs))
    badh <- logical(length(xs))
    for (j in seq_along(xs)) {
      y <- xs[j]
      tp <- RTMB::MakeTape(function(p) famg$lpdf(y, list(mu = p[1],
        sigma = exp(p[2]), xi = expm1(p[3]), .eta_sigma = p[2],
        .eta_xi = p[3]), list()), p)
      gr <- tp$jacobian(p)
      bad[j] <- any(!is.finite(gr))
      if (j %% 10 == 0) {
        hs <- tp$jacfun()$jacobian(p)
        badh[j] <- any(!is.finite(hs))
      }
    }
    cat(sprintf(paste0("gevgrad xi %.2f: x from %.2f to %.2f step 0.01: ",
                       "%d of %d with a non-finite gradient%s; Hessian ",
                       "non-finite at %d of %d probed\n"),
                xi0, max(xs), min(xs), sum(bad), length(xs),
                if (any(bad)) paste0(" (x in ", paste(range(xs[bad]),
                                                      collapse = " .. "),
                                    ")") else "",
                sum(badh), sum(seq_along(xs) %% 10 == 0)))
  }
  # the frechet's capped term on the same band: nu z near -700
  famf <- frechet()
  bad <- 0
  for (lz in seq(-0.5, -0.9, length.out = 401)) {
    y <- exp(lz * 1000 / 3)
    tp <- RTMB::MakeTape(function(p) famf$lpdf(y, list(mu = exp(p[1]),
      nu = 1 + exp(p[2]), .eta_mu = p[1], .eta_nu = p[2]), list()),
      c(0, log(2)))
    bad <- bad + any(!is.finite(tp$jacobian(c(0, log(2)))))
  }
  cat("gevgrad frechet: rows with log y from", -0.5 * 1000 / 3, "to",
      -0.9 * 1000 / 3, ":", bad, "of 401 with a non-finite gradient\n")
}
if (want("exg")) {
  # 3003: exgaussian, group sd 30 times sigma (dev/optstart-exg.R)
  for (s in 1:10) {
    set.seed(s)
    n <- 400
    g <- factor(rep(1:20, length.out = n))
    x <- rnorm(n)
    bg <- rnorm(20, 0, 30)
    mu <- 0.3 + 0.5 * x + bg[g]
    y <- ns$family_registry$exgaussian()$sim(
      list(mu = mu, sigma = rep(1, n), beta = rep(1.5, n)), list(), n)
    d <- data.frame(y = y, x = x, g = g)
    r <- cap(frm(y ~ x + (1 | g), data = d, family = exgaussian()))
    f <- r$v
    if (inherits(f, "err")) {
      cat(sprintf("exg seed %d: ERROR %s\n", s, substr(f, 1, 120)))
    } else {
      cat(sprintf(paste0("exg seed %d: code %d logLik %.6f sigma %.4f ",
                         "beta %.4f sd_g %.3f warnings %d\n"), s,
                  f$opt$convergence, as.numeric(logLik(f)),
                  exp(f$estimates$betad[1]), exp(f$estimates$betad[2]),
                  exp(f$estimates$theta[1]), length(r$w)))
    }
  }
}
if (want("gprior")) {
  # 3560: gamma(2, 1) on a log-scale b of shape ~ 0 + fac
  # (dev/fuzzfix-rev-r2-edges.R)
  set.seed(20261021)
  n <- 80
  d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                  fac = factor(rep(c("a", "b"), 40)),
                  w = runif(n, 0.5, 2), g = factor(rep(1:20, each = 4)))
  eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
  d$y <- eta + rnorm(n, 0, 0.8)
  d$yg <- rgamma(n, shape = 3, rate = 3 / exp(eta))
  fg <- bf(yg | weights(w) ~ x + (1 | obs), shape ~ 0 + fac)
  brief("gprior: Gamma shape ~ 0 + fac, gamma(2, 1) on b (3560)",
        cap(frm(fg, family = Gamma("log"), data = d,
                prior = set_prior("gamma(2, 1)", class = "b",
                                  dpar = "shape"))))
  brief("gprior: the same with lb = 0",
        cap(frm(fg, family = Gamma("log"), data = d,
                prior = set_prior("gamma(2, 1)", class = "b",
                                  dpar = "shape", lb = 0))))
  brief("gprior: the same, normal(0, 1) (guard absent)",
        cap(frm(fg, family = Gamma("log"), data = d,
                prior = set_prior("normal(0, 1)", class = "b",
                                  dpar = "shape"))))
  # no weights, no (1 | obs): a positive-support prior whose mode is
  # inside, on a coefficient whose estimate is positive
  fg2 <- bf(yg ~ x, shape ~ 0 + fac)
  brief("gprior: no weights, gamma(2, 1) on b",
        cap(frm(fg2, family = Gamma("log"), data = d,
                prior = set_prior("gamma(2, 1)", class = "b",
                                  dpar = "shape"))))
  # a coefficient whose likelihood maximum is negative
  d$yn <- rgamma(n, shape = 0.4, rate = 0.4 / exp(eta))
  brief("gprior: shape 0.4 data, gamma(2, 1) on b",
        cap(frm(bf(yn ~ x, shape ~ 0 + fac), family = Gamma("log"),
                data = d, prior = set_prior("gamma(2, 1)", class = "b",
                                            dpar = "shape"))))
}
if (want("nbmirror")) {
  # 3644: negbinomial vs negbinomial2 at the Poisson limit
  # (dev/nb2dm-rev-mirror.R)
  for (seed in c(4, 6, 1)) {
    set.seed(seed)
    n <- 300
    d <- data.frame(x = rnorm(n))
    d$y <- rpois(n, exp(1 + 0.3 * d$x))
    r1 <- cap(frm(y ~ x, family = negbinomial(), data = d))
    r2 <- cap(frm(y ~ x, family = brmsfamily("negbinomial2"), data = d))
    for (k in 1:2) {
      r <- list(r1, r2)[[k]]
      f <- r$v
      cat(sprintf("nbmirror seed %d %s: code %s logLik %.10f last par %.4f\n",
                  seed, c("nb ", "nb2")[k], f$opt$convergence,
                  as.numeric(logLik(f)), f$opt$par[length(f$opt$par)]))
      showw(r, 200)
    }
  }
}
if (want("csmix")) {
  # 2750: the cs() ordinal-mixture wall (dev/optima-csmix-wall.R)
  mk <- function(seed, n = 300) {
    set.seed(seed)
    d <- data.frame(x = rnorm(n), z = rnorm(n))
    cls <- rbinom(n, 1, 0.4)
    lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
    d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
    d
  }
  for (s in 1:20) {
    d <- mk(s)
    r <- cap(frm(bf(y ~ x + cs(z)), family = mixture(cumulative(), sratio()),
                 data = d))
    f <- r$v
    if (inherits(f, "err")) {
      cat(sprintf("csmix seed %d: ERROR %s\n", s, substr(f, 1, 120)))
      next
    }
    pl <- f$estimates
    tr <- pl$tau_raw1
    tau <- cumsum(c(tr[1], exp(tr[-1])))
    eta <- pl$beta[1] * d$x
    A <- outer(-eta, tau, "+") - outer(d$z, pl$bcs3)
    lo <- cbind(-Inf, A)
    up <- cbind(A, Inf)
    gap <- up[cbind(seq_len(nrow(d)), d$y)] - lo[cbind(seq_len(nrow(d)), d$y)]
    ow <- r$w[!grepl("Category specific effects for this family", r$w,
                     fixed = TRUE)]
    named <- any(grepl("cross|closes|closed", c(ow, r$m)))
    cat(sprintf(paste0("csmix seed %d: code %d logLik %.4f min gap %.3g ",
                       "warnings %d named %s\n"), s, f$opt$convergence,
                as.numeric(logLik(f)), min(gap), length(ow), named))
    for (x in unique(c(ow, r$m))) {
      cat("    ", substr(gsub("\\s+", " ", x), 1, 160), "\n")
    }
  }
}
cat("\ndone\n")
