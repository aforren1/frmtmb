# rel084-kt copy of dev/reviews/alapse-review-rq.R: sources dev/rel084-kt/alapse-review-common.R (rellib-r23). No other change.
# Review of lane alapse, check 1: is ald_polish()'s location exactly the
# weighted quantile regression (rq br and fn), on weights, ties, a
# non-unique solution, a rank-deficient X, an offset; is the "proved"
# point a local maximum (random perturbations of every outer parameter).
#   Rscript dev/reviews/alapse-review-rq.R <base|lane>
lib_arg <- commandArgs(TRUE)[1]
source("dev/rel084-kt/alapse-review-common.R")

one <- function(id, fo, d, X, y, w = rep(1, length(y)), off = 0,
                fam = asym_laplace(), probe = TRUE) {
  cat("==", id, "\n")
  r <- cap(frm(fo, family = fam, data = d))
  say("frm", r)
  if (inherits(r$v, "error")) return(invisible())
  f <- r$v
  cert <- f$opt$ald_cert
  cat(sprintf("  code %d | cert ok %s | why %s | ll %.10f\n",
              f$opt$convergence, format(isTRUE(cert$ok)),
              cert$why %||% "-", as.numeric(logLik(f))))
  if (isTRUE(cert$ok)) {
    cat(sprintf("  cert: nlminb code %s gain %.3g margin %.3g gap %.3g moved %.3g rounds %s\n",
                format(cert$code), cert$gain, cert$margin, cert$gap,
                cert$moved, format(cert$rounds)))
  }
  dp <- ns$eval_dpars(f)[[1L]]
  s <- rep_len(dp$sigma, length(y))
  p <- rep_len(dp$quantile, length(y))
  b <- ns$fixef_by_dpar(f)$mu
  cat("  frm b:", format(signif(b, 10)), "\n")
  yo <- y - off
  obj <- function(bb) sum(w / s * (yo - drop(X %*% bb)) *
                            (p - ((yo - drop(X %*% bb)) < 0)))
  of <- obj(b[colnames(X)] %||% b)
  if (length(unique(p)) == 1L) {
    for (m in c("br", "fn")) {
      rr <- tryCatch(suppressWarnings(
        quantreg::rq.wfit(X, yo, tau = p[1], weights = w / s, method = m)),
        error = function(e) NULL)
      if (is.null(rr)) { cat("  rq", m, "failed\n"); next }
      bq <- rr$coefficients
      cat(sprintf("  rq %s: obj frm %.12g rq %.12g rel(frm-rq) %.3g | max|db| %.3g\n",
                  m, of, obj(bq), (of - obj(bq)) / abs(obj(bq)),
                  max(abs(b - bq), na.rm = TRUE)))
    }
    nu <- tryCatch(withCallingHandlers(
      quantreg::rq.wfit(X, yo, tau = p[1], weights = w / s, method = "br"),
      warning = function(z) { cat("  rq br says:", conditionMessage(z), "\n")
        invokeRestart("muffleWarning") }), error = function(e) NULL)
  }
  if (probe && identical(lib_arg, "lane")) {
    lp <- local_probe(f)
    cat("  local probe (max gain of -obj by step):",
        paste(names(lp), format(lp, digits = 3), sep = "=", collapse = " "),
        "\n")
  }
  invisible(f)
}

# 1 weights, fixed quantile 0.3 and estimated
set.seed(11)
n <- 150
d <- data.frame(x = rnorm(n), w = runif(n, 0.2, 3))
d$y <- brms::rasym_laplace(n, 1 + 0.5 * d$x, 1, 0.3)
X <- cbind(`(Intercept)` = 1, x = d$x)
one("weights-fixq", bf(y | weights(w) ~ x, quantile = 0.3), d, X, d$y, d$w)
one("weights-estq", y | weights(w) ~ x, d, X, d$y, d$w)

# 2 ties: integer y, discrete x, median: rq's solution is not unique
for (s in 1:3) {
  set.seed(20 + s)
  n <- 60
  d <- data.frame(x = sample(0:2, n, TRUE))
  d$y <- rpois(n, 3 + d$x)
  X <- cbind(`(Intercept)` = 1, x = d$x)
  one(paste0("ties-fixq0.5-s", s), bf(y ~ x, quantile = 0.5), d, X, d$y)
  one(paste0("ties-estq-s", s), y ~ x, d, X, d$y)
}

# 3 intercept only, even n, median: any value between the middle two
set.seed(31)
d <- data.frame(y = round(rnorm(10), 3))
X <- cbind(`(Intercept)` = rep(1, 10))
f3 <- one("even-median", bf(y ~ 1, quantile = 0.5), d, X, d$y)
cat("  middle two:", sort(d$y)[5:6], "\n")

# 4 non-unique face with a row-varying scale: discrete y and x, sigma ~ g
for (s in 1:6) {
  set.seed(40 + s)
  n <- 80
  d <- data.frame(x = sample(0:1, n, TRUE), g = factor(sample(c("a", "b"), n,
                                                              TRUE)))
  d$y <- rpois(n, 2 + d$x) * ifelse(d$g == "a", 1, 2)
  X <- cbind(`(Intercept)` = 1, x = d$x)
  one(paste0("ties-sigmag-s", s), bf(y ~ x, sigma ~ g, quantile = 0.5), d,
      X, d$y)
  one(paste0("ties-sigmag-quantg-s", s), bf(y ~ x, sigma ~ g, quantile ~ g),
      d, X, d$y, probe = TRUE)
}

# 5 rank-deficient X
set.seed(51)
n <- 100
d <- data.frame(x1 = rnorm(n))
d$x2 <- 2 * d$x1
d$y <- 1 + d$x1 + rnorm(n)
X <- cbind(`(Intercept)` = 1, x1 = d$x1, x2 = d$x2)
r5 <- cap(frm(bf(y ~ x1 + x2, quantile = 0.5), family = asym_laplace(),
              data = d))
cat("== rank-deficient\n")
say("frm", r5)
if (!inherits(r5$v, "error")) {
  f <- r5$v
  cat("  code", f$opt$convergence, "| cert", format(isTRUE(f$opt$ald_cert$ok)),
      "|", f$opt$ald_cert$why %||% "-", "| ll", as.numeric(logLik(f)), "\n")
  print(ns$fixef_by_dpar(f)$mu)
  rr <- quantreg::rq(y ~ x1, tau = 0.5, data = d)
  s0 <- ns$eval_dpars(f)[[1]]$sigma[1]
  cat("  ALD ll at rq(y ~ x1) with fit sigma:",
      ald_ll(d$y, fitted(rr), s0, 0.5), "\n")
}

# 6 offset
set.seed(61)
n <- 120
d <- data.frame(x = rnorm(n), o = runif(n, -2, 2))
d$y <- 1 + 0.5 * d$x + d$o + rnorm(n)
X <- cbind(`(Intercept)` = 1, x = d$x)
one("offset-fixq", bf(y ~ x + offset(o), quantile = 0.7), d, X, d$y,
    off = d$o)

# 7 multivariate, ALD and gaussian responses
set.seed(71)
n <- 100
d <- data.frame(x = rnorm(n))
d$y1 <- 1 + d$x + rnorm(n)
d$y2 <- 2 - d$x + rnorm(n)
r7 <- cap(frm(bf(y1 ~ x, family = asym_laplace()) +
                bf(y2 ~ x, family = gaussian()), data = d))
cat("== mv\n")
say("frm", r7)
if (!inherits(r7$v, "error")) {
  f <- r7$v
  cat("  code", f$opt$convergence, "| cert", format(isTRUE(f$opt$ald_cert$ok)),
      "|", f$opt$ald_cert$why %||% "-", "\n")
  print(ns$sdr_of(f)$se_lost)
}
cat("DONE\n")
