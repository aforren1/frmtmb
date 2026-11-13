# Lane boundfam: the defects of the brief, and the guard-absent cases,
# on one build.
#   Rscript dev/boundfam-repro.R <lib or "base"> [case ...]
# "base" reads rellib-r20 alone; any other value puts that library first.
a <- commandArgs(trailingOnly = TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
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
showw <- function(r, n = 600) {
  for (x in r$w) cat("   W:", substr(gsub("\\s+", " ", x), 1, n), "\n")
  for (x in r$m) cat("   M:", substr(gsub("\\s+", " ", x), 1, n), "\n")
  if (inherits(r$v, "err")) cat("   ERR:", substr(r$v, 1, n), "\n")
}
report <- function(label, r) {
  cat("\n==", label, "\n")
  showw(r)
  f <- r$v
  if (inherits(f, "err")) return(invisible(NULL))
  cat("   code", f$opt$convergence, " logLik",
      format(as.numeric(logLik(f)), digits = 10), "\n")
  s <- cap(sqrt(diag(vcov(f, full = TRUE))))
  showw(s)
  p <- f$opt$par
  nm <- ns$outer_par_names(f)
  se <- if (!inherits(s$v, "err")) s$v else rep(NA, length(p))
  lost <- tryCatch(ns$sdr_of(f)$se_lost, error = function(e) NULL)
  why <- if (length(lost)) unname(lost[nm]) else rep(NA, length(nm))
  why[is.na(why)] <- ""
  print(data.frame(par = nm, est = signif(p, 6), se = signif(se, 4),
                   lost = why), row.names = FALSE)
  invisible(f)
}

if (want("pois0")) {
  # 3333: no MLE, every count zero
  set.seed(2)
  de <- data.frame(x = rnorm(50), y = 0L)
  report("pois0: poisson, all-zero response (3333)",
         cap(frm(y ~ x, family = poisson(), data = de)))
  print(cap(suppressWarnings(glm(y ~ x, family = poisson(), data = de)))$w)
  g <- cap(glm(y ~ x, family = poisson(), data = de))
  cat("   glm warnings:", g$w, "\n")
  print(summary(g$v)$coefficients[, 1:2])
}
if (want("pois1")) {
  # guard absent: one nonzero count at an interior x has a maximum
  set.seed(2)
  de <- data.frame(x = rnorm(50), y = 0L)
  de$y[order(de$x)[25]] <- 1L
  report("pois1: poisson, one count at the median x (guard absent)",
         cap(frm(y ~ x, family = poisson(), data = de)))
  g <- cap(glm(y ~ x, family = poisson(), data = de))
  cat("   glm warnings:", g$w, "\n")
  print(summary(g$v)$coefficients[, 1:2])
  # one count at the largest x: the slope runs away (no MLE)
  de$y <- 0L
  de$y[which.max(de$x)] <- 1L
  report("pois1max: poisson, one count at the largest x (no MLE)",
         cap(frm(y ~ x, family = poisson(), data = de)))
  g <- cap(glm(y ~ x, family = poisson(), data = de))
  cat("   glm warnings:", g$w, "\n")
  print(summary(g$v)$coefficients[, 1:2])
}
if (want("cens")) {
  # 3333: every row right-censored
  for (s in 1:6) {
    set.seed(s)
    dc <- data.frame(x = rnorm(60))
    dc$y <- 1 + dc$x + rnorm(60)
    dc$cen <- "right"
    report(paste0("cens seed ", s, ": every row right-censored (3333)"),
           cap(frm(y | cens(cen) ~ x, data = dc)))
    if (s %in% 2:3) {
      sv <- cap(survival::survreg(survival::Surv(y, rep(0, 60)) ~ x,
                                  data = dc, dist = "gaussian"))
      cat("   survreg warnings:", sv$w, "\n")
      if (!inherits(sv$v, "err")) {
        print(summary(sv$v)$table)
      } else cat("   survreg ERR", sv$v, "\n")
    }
  }
}
if (want("maling")) {
  # 2710: Malingering_2, component-2 precision runs to the binomial limit
  d <- data.frame(k = c(45, 45, 44, 45, 44, 45, 45, 45, 45, 45, 30, 20, 6,
                        44, 44, 27, 25, 17, 14, 27, 35, 30), n = 45L)
  report("maling: mixture(beta_binomial, beta_binomial) (2710)",
         cap(frm(bf(k | trials(n) ~ 1),
                 family = mixture(beta_binomial, beta_binomial), data = d,
                 start = list(beta = qlogis(c(0.55, 0.95)),
                              betad = c(log(c(15, 15)), 0)))))
}
if (want("frechet")) {
  # 3022: frechet at nu0 = 1.1, seed 31 of dev/famuni-rev-frechet11.R
  nu0 <- 1.1
  set.seed(9000 + 31)
  d <- data.frame(x = rnorm(300))
  d$y <- exp(0.5 + 0.3 * d$x) / gamma(1 - 1 / nu0) * rexp(300)^(-1 / nu0)
  report("frechet: nu0 = 1.1, seed 9031 (3022)",
         cap(frm(y ~ x, family = frechet(), data = d)))
}
if (want("gev")) {
  # 3037: xi in (-1, -0.5)
  rg <- function(n, mu, s, xi) mu + s * ((-log(runif(n)))^(-xi) - 1) / xi
  for (xi in c(-0.63, -0.7, -0.89, -0.3)) {
    set.seed(4)
    dv <- data.frame(y = rg(300, 0, 1, xi))
    report(paste0("gev: xi0 = ", xi, ", seed 4 (3037)"),
           cap(frm(y ~ 1, family = gen_extreme_value(), data = dv)))
  }
}
if (want("wobs")) {
  # 3482: weights summing to exactly n with an observation-level effect
  for (s in 1:3) {
    set.seed(s)
    n <- 60
    d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)))
    d$y <- 1 + d$x + rnorm(n, 0, 0.5) + rnorm(n, 0, 0.5)
    w <- runif(n, 0.5, 2)
    d$w <- w / mean(w)
    report(paste0("wobs gaussian seed ", s, ": sum(w) = n, (1 | obs) (3482)"),
           cap(frm(y | weights(w) ~ x + (1 | obs), data = d)))
    d$yg <- rgamma(n, shape = 5, rate = 5 / exp(0.5 + 0.3 * d$x +
                                                  rnorm(n, 0, 0.3)))
    report(paste0("wobs Gamma seed ", s, ": sum(w) = n, (1 | obs) (3482)"),
           cap(frm(yg | weights(w) ~ x + (1 | obs), family = Gamma("log"),
                   data = d)))
  }
}
if (want("eam")) {
  # 2591 = 2769: frmtmb.eam test-defects.R:69's mixture, weight at 6e-11
  suppressPackageStartupMessages(library(frmtmb.eam))
  set.seed(3)
  dat <- ddm_simulate(300, mu = 0.9, bs = 1.4, ndt = 0.30)
  k <- sample(300, 18)
  dat$rt[k] <- stats::runif(18, 0.12, 2.5)
  report("eam: mixture(wiener(max_ndt = 0.4), lognormal()) (2591)",
         cap(frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
                 family = mixture(wiener(max_ndt = 0.4,
                                         allow_unreachable = TRUE),
                                  lognormal()),
                 data = dat)))
}
if (want("negbin")) {
  # a negbinomial fit to poisson data: the shape has no finite maximum
  for (s in 1:3) {
    set.seed(s)
    d <- data.frame(x = rnorm(200))
    d$y <- rpois(200, exp(1 + 0.4 * d$x))
    report(paste0("negbin on poisson data, seed ", s),
           cap(frm(y ~ x, family = negbinomial(), data = d)))
    m <- cap(MASS::glm.nb(y ~ x, data = d))
    cat("   glm.nb warnings:", m$w, " theta", m$v$theta, "SE",
        m$v$SE.theta, "\n")
  }
}
if (want("absent")) {
  # family parameters near (not at) a boundary, with a real maximum
  set.seed(11)
  n <- 400
  d <- data.frame(x = rnorm(n))
  d$y <- MASS::rnegbin(n, exp(1 + 0.5 * d$x), theta = 30)
  report("absent nb: negbinomial shape 30, n = 400",
         cap(frm(y ~ x, family = negbinomial(), data = d)))
  d$k <- rbinom(n, 20, plogis(0.3 + 0.4 * d$x))
  report("absent bb: beta_binomial on binomial data, n = 400",
         cap(frm(k | trials(20) ~ x, family = beta_binomial(), data = d)))
  z <- rbinom(n, 1, 0.08)
  d$ym <- ifelse(z == 1, rnorm(n, 4, 1), rnorm(n, 0, 1))
  report("absent mix: mixture(gaussian, gaussian), weight 0.08",
         cap(frm(bf(ym ~ 1, mu2 ~ 1), family = mixture(gaussian, gaussian),
                 data = d)))
  d$yl <- exp(1 + 0.2 * d$x + rnorm(n, 0, 0.02))
  report("absent lognormal: sigma 0.02",
         cap(frm(yl ~ x, family = lognormal(), data = d)))
  d$yb <- rbinom(n, 1, plogis(-4 + 0.5 * d$x))
  report("absent bern: rare bernoulli",
         cap(frm(yb ~ x, family = bernoulli(), data = d)))
}
cat("\ndone\n")
