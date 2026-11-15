# Lane selimit2: the items of the brief, on one build.
#   Rscript dev/selimit2-probe.R <lib or "base"> [case ...]
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
showw <- function(r, n = 700) {
  for (x in r$w) cat("   W:", substr(gsub("[[:space:]]+", " ", x), 1, n), "\n")
  for (x in r$m) cat("   M:", substr(gsub("[[:space:]]+", " ", x), 1, n), "\n")
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
fitted_report <- function(label, fit, newdata = NULL) {
  r <- cap(if (is.null(newdata)) fitted(fit) else
    fitted(fit, newdata = newdata))
  if (inherits(r$v, "err")) {
    cat(sprintf("   %-40s ERROR %s\n", label, r$v))
    return(invisible())
  }
  se <- r$v[, "Est.Error"]
  cat(sprintf("   %-40s rows %4d NaN %4d est %s se %s warn %d %s\n", label,
              length(se), sum(!is.finite(se)),
              format(signif(r$v[1, "Estimate"], 4)),
              format(signif(se[1], 4)), length(r$w),
              if (length(r$w)) substr(r$w[1], 1, 60) else ""))
}

if (want("pois1max")) {
  # 4800: one count at the largest x, a joint runaway
  set.seed(2)
  de <- data.frame(x = rnorm(50), y = 0L)
  de$y[which.max(de$x)] <- 1L
  report("pois1max: poisson, one count at the largest x",
         cap(frm(y ~ x, family = poisson(), data = de)))
  # guard absent: the one count at the median x has a maximum
  de$y <- 0L
  de$y[order(de$x)[25]] <- 1L
  report("pois1: poisson, one count at the median x (absent)",
         cap(frm(y ~ x, family = poisson(), data = de)))
}
if (want("grp0")) {
  # 4800, the coding dependence
  set.seed(5)
  d <- data.frame(g = factor(rep(c("a", "b"), each = 50)))
  d$y <- ifelse(d$g == "a", 0L, rpois(100, 3))
  report("grp0 poisson y ~ g", cap(frm(y ~ g, family = poisson(), data = d)))
  report("grp0 poisson y ~ 0 + g",
         cap(frm(y ~ 0 + g, family = poisson(), data = d)))
  report("grp0 zip y ~ g",
         cap(frm(y ~ g, family = zero_inflated_poisson(), data = d)))
  report("grp0 negbinomial y ~ g",
         cap(frm(y ~ g, family = negbinomial(), data = d)))
  # guard absent: group a has a single count
  d$y[1] <- 1L
  report("grp0 absent: group a one count, y ~ g",
         cap(frm(y ~ g, family = poisson(), data = d)))
}
if (want("zif")) {
  # 4800 (treatment-coded zi ~ f) and 4814 (zi ~ 0 + f)
  set.seed(21)
  n <- 300
  d <- data.frame(f = factor(rep(c("a", "b"), each = n / 2)), x = rnorm(n))
  d$y1 <- ifelse(d$f == "a", 0L, rpois(n, 4))
  report("zi ~ 0 + f, a all zero",
         cap(frm(bf(y1 ~ 1, zi ~ 0 + f), family = zero_inflated_poisson(),
                 data = d)))
  report("zi ~ f, a all zero",
         cap(frm(bf(y1 ~ 1, zi ~ f), family = zero_inflated_poisson(),
                 data = d)))
  d$z1 <- as.integer(d$y1 == 0)
  report("bernoulli analogue z ~ 0 + f",
         cap(frm(z1 ~ 0 + f, family = bernoulli(), data = d)))
  # hu: group a all zero, hurdle
  report("hu ~ f, a all zero",
         cap(frm(bf(y1 ~ 1, hu ~ f), family = hurdle_poisson(), data = d)))
  # guard absent: zi ~ 0 + f where group b has no excess zeros (zi_fb to
  # -Inf, the family limit on b's rows) and group a has some
  d$y2 <- ifelse(d$f == "a" & runif(n) < 0.4, 0L, rpois(n, 4))
  report("zi ~ 0 + f, a 40% zeros (absent)",
         cap(frm(bf(y2 ~ 1, zi ~ 0 + f), family = zero_inflated_poisson(),
                 data = d)))
}
if (want("zifq")) {
  # a joint end that is the family's limit: level a (the reference) has
  # no excess zeros, level b has 40 percent, so zi_(Intercept) runs to
  # -Inf and zi_fb to +Inf together; group b's zi stays determined
  for (s in 1:3) {
    set.seed(30 + s)
    n <- 300
    d <- data.frame(f = factor(rep(c("a", "b"), each = n / 2)))
    d$y <- ifelse(d$f == "b" & runif(n) < 0.4, 0L, rpois(n, 4))
    f <- report(paste("zi ~ f, a poisson, b inflated, seed", 30 + s),
                cap(frm(bf(y ~ 1, zi ~ f), family = zero_inflated_poisson(),
                        data = d)))
    if (!is.null(f)) {
      fitted_report("fitted, data", f)
      fitted_report("fitted, level a", f, data.frame(f = "a"))
      fitted_report("fitted, level b", f, data.frame(f = "b"))
    }
  }
}
if (want("zix")) {
  # 5433 and zi ~ x seed 5's "flat"
  zib <- function(seed) {
    set.seed(seed)
    n <- 200
    d <- data.frame(x = stats::rnorm(n), N = sample(5:15, n, TRUE))
    d$y <- stats::rbinom(n, d$N, stats::plogis(0.3 + 0.5 * d$x))
    d
  }
  for (s in c(5, 7)) {
    d <- zib(s)
    f <- report(paste0("zib seed ", s, " zi ~ x"),
                cap(frm(bf(y | trials(N) ~ x, zi ~ x),
                        family = zero_inflated_binomial(), data = d)))
    fitted_report("fitted, data", f)
    for (x in c(-100, -30, -10, -5, 0, 3, 10)) {
      fitted_report(sprintf("fitted newdata x = %g", x), f,
                    data.frame(x = x, N = 10))
    }
    if (s == 7) {
      # how well a zi slope of -3 (and others) fits: the profile of zi_x
      p <- f$opt$par
      nm <- ns$outer_par_names(f)
      j <- match("zi_x", nm)
      for (b in c(-7, -3, -1, -0.32, 1, 3, 7)) {
        q <- p
        q[j] <- b
        cat(sprintf("   zi_x = %5.2f  obj change %.3g\n", b,
                    f$obj$fn(q) - f$obj$fn(p)))
      }
    }
  }
}
if (want("nb17")) {
  # 4824
  set.seed(17)
  d <- data.frame(x = rnorm(200))
  d$y <- rpois(200, exp(1 + 0.4 * d$x))
  f <- report("nb17 negbinomial on poisson data",
              cap(frm(y ~ x, family = negbinomial(), data = d)))
  cat("   max|grad|", max(abs(f$obj$gr(f$opt$par))), "\n")
  fitted_report("fitted", f)
}
if (want("asym")) {
  # 5423, 3389, 3575
  set.seed(5)
  g <- factor(rep(1:10, each = 10))
  x <- stats::rnorm(100)
  dd <- data.frame(g = g, x = x, y = 1 + 0.5 * x + stats::rnorm(100))
  f <- report("asym_laplace y ~ x", cap(frm(y ~ x, data = dd,
                                           family = asym_laplace())))
  fitted_report("fitted", f)
  f <- report("asym_laplace y ~ x, quantile 0.5 fixed",
              cap(frm(bf(y ~ x, quantile = 0.5), data = dd,
                      family = asym_laplace())))
  fitted_report("fitted", f)
}
