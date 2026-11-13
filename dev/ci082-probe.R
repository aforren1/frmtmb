# Lane ci082: the three BLAS-dependent CI failures of 0.82.0, one case
# per process.
#   Rscript dev/ci082-probe.R <lib|base> <case> [k]
# lib: the lane library before rellib-r21; base: rellib-r21 only.
# case: xcs (test-covstruct-seam.R:126), mix (test-mixture-starts.R:333),
# nb (test-se-limit.R:98). k perturbs the data by 2e-15 * k (ulps).
args <- commandArgs(TRUE)
arm <- args[1]
case <- args[2]
k <- if (length(args) >= 3) as.numeric(args[3]) else 0
base <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/wt-ci082-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
cat("R:", R.home(), " CORETYPE:", Sys.getenv("OPENBLAS_CORETYPE"), "\n")
cat("La:", La_library(), "\n")

conds <- function(expr) {
  w <- list()
  m <- list()
  v <- withCallingHandlers(expr, warning = function(x) {
    w[[length(w) + 1L]] <<- x
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m[[length(m) + 1L]] <<- x
    invokeRestart("muffleMessage")
  })
  list(value = v, warnings = w, messages = m)
}
show <- function(r) {
  for (x in r$warnings) cat("  W [", class(x)[1], "]", conditionMessage(x),
                            "\n")
  for (x in r$messages) cat("  M [", class(x)[1], "]",
                            substr(conditionMessage(x), 1, 300), "\n")
}

if (case %in% c("xcs", "xcsk")) {
  ks <- if (case == "xcsk") 0:19 else k
} else ks <- numeric()
# (a for loop over nothing sets its variable to NULL, so k is not reused)
for (kx in ks) if (case %in% c("xcs", "xcsk")) {
  k <- kx
  xcs_C <- function(d) matrix(0.5, d, d) + diag(0.5, d)
  frmtmb_register_covstruct(
    "xcs",
    lpdf = function(b, theta, blk) {
      d <- blk[["dim"]]
      Sigma <- exp(2 * theta[1]) * xcs_C(d)
      dim(b) <- c(d, length(b) %/% d)
      sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
    },
    vcov = function(theta, blk) {
      V <- exp(2 * theta[1]) * xcs_C(blk[["dim"]])
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    npar = 1L, sd_idx = 1L)
  set.seed(11L)
  G <- 30L
  d <- expand.grid(t = factor(c("a", "b", "c")), rep = 1:3,
                   g = factor(sprintf("g%02d", seq_len(G))))
  u <- t(crossprod(chol(0.8^2 * xcs_C(3L)), matrix(rnorm(3 * G), 3)))
  d$x <- rnorm(nrow(d))
  d$y <- 0.4 + 0.7 * d$x + u[cbind(as.integer(d$g), as.integer(d$t))] +
    rnorm(nrow(d), 0, 0.6)
  d$y <- d$y * (1 + 2e-15 * k)
  fit <- frm(bf(y ~ x + xcs(0 + t | g)), family = gaussian(), data = d)
  cat(sprintf("theta %.17g\n", fit$estimates[["theta"]][1]))
  tr <- frmtmb:::varcorr_trans_rows(fit)
  print(format(tr, digits = 17))
  r <- conds(confint_varcorr(fit))
  show(r)
  cv <- r$value
  cat("width cor:", format(cv$upr[4:6] - cv$lwr[4:6], digits = 17), "\n")
  cat("RESULT xcs k", k, "se_cor_zero",
      paste(tr$se_t[4:6] == 0, collapse = ","), "na",
      paste(is.na(cv$upr[4:6]), collapse = ","), "\n")
}

if (case == "mixk") {
  # k = 0..K-1 in one process: does each arm warn of a step function?
  oms_data <- function(seed, n = 300, k = 0) {
    set.seed(seed)
    d <- data.frame(x = rnorm(n), z = rnorm(n))
    cls <- stats::rbinom(n, 1, 0.4)
    lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
    d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
    d$x <- d$x * (1 + 2e-15 * k)
    d
  }
  fam <- mixture(cumulative(), cumulative())
  # args: [3] first k, [4] how many
  for (kk in k + seq_len(if (length(args) >= 4) as.integer(args[4]) else
    10) - 1L) {
    d <- oms_data(4, k = kk)
    d$xs <- d$x * 1e-4
    for (f in list(bf(y ~ x + z), bf(y ~ xs + z))) {
      r <- conds(frm(f, family = fam, data = d))
      fit <- r$value
      ms <- fit$opt$mixture_starts
      step <- any(vapply(r$warnings, inherits, NA,
                         "frmtmb_mixture_step_start"))
      cat(sprintf(paste("RESULT mixk k %d scaled %s code %s logLik %.6f",
                        "reached %s gain %.4g above %s tier %s step %s",
                        "nwarn %d\n"),
                  kk, isTRUE(fit$autoscaled), fit$opt$convergence,
                  as.numeric(logLik(fit)), ms[["reached"]], ms[["gain"]],
                  ms[["degenerate_above"]], ms[["tier"]], step,
                  length(r$warnings)))
      for (x in r$warnings) {
        if (!inherits(x, "frmtmb_mixture_step_start")) {
          cat("  W", substr(conditionMessage(x), 1, 200), "\n")
        }
      }
    }
  }
}

if (case == "mix") {
  oms_data <- function(seed, n = 300, k = 0) {
    set.seed(seed)
    d <- data.frame(x = rnorm(n), z = rnorm(n))
    cls <- stats::rbinom(n, 1, 0.4)
    lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
    d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
    d$x <- d$x * (1 + 2e-15 * k)
    d
  }
  fam <- mixture(cumulative(), cumulative())
  d <- oms_data(4, k = k)
  d$xs <- d$x * 1e-4
  for (f in list(bf(y ~ x + z), bf(y ~ xs + z))) {
    r <- conds(frm(f, family = fam, data = d))
    fit <- r$value
    ms <- fit$opt$mixture_starts
    step <- any(vapply(r$warnings, inherits, NA,
                       "frmtmb_mixture_step_start"))
    cat(sprintf("formula %s autoscaled %s code %s logLik %.10f\n",
                deparse(f$formula), isTRUE(fit$autoscaled),
                fit$opt$convergence, as.numeric(logLik(fit))))
    cat("  starts:", paste(names(ms), unlist(lapply(ms, function(v) {
      if (length(v) == 1L) format(v, digits = 10) else
        paste(format(v, digits = 6), collapse = "/")
    })), sep = "=", collapse = " "), "\n")
    cl <- frmtmb:::mixture_ord_start_class(fit, "y")
    cat("  class step:", cl$step, " boundary:", cl$boundary, "\n")
    show(r)
    cat("RESULT mix k", k, "scaled", isTRUE(fit$autoscaled), "step", step,
        "\n")
  }
}

if (case == "nlgp") {
  # test-se-limit.R "a family limit on a nonlinear fit still triggers the
  # restart", timed
  n <- 120L
  set.seed(81100 + 1000 * 4 + 100 * 2 + n + 10 * 6)
  x <- sort(round(stats::runif(n, 0, 10), 3))
  z <- round(stats::rnorm(n), 3)
  w <- round(stats::runif(n, 0.2, 5), 3)
  f <- 2 * sin(3 * x) + x
  fs <- f / max(abs(f))
  y <- (2 + fs) * w^0.5 + 0.1 * stats::rt(n, 3)
  d <- data.frame(x, z, w, y)
  t0 <- proc.time()
  r <- conds(frm(bf(y ~ a * w^b, a ~ gp(x), b ~ 1, nl = TRUE),
                 data = d, family = student(),
                 control = frmtmb_control(verbose = 1)))
  print(proc.time() - t0)
  show(r)
  cat("RESULT nlgp logLik", as.numeric(logLik(r$value)), "\n")
}

if (case == "nbone") {
  # one perturbation k of the negbinomial case, looked at closely
  set.seed(2)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- stats::rpois(200, exp(1 + 0.4 * d$x))
  d$x <- d$x * (1 + 2e-15 * k)
  r <- conds(frm(y ~ x, family = negbinomial(), data = d))
  fit <- r$value
  show(r)
  cat(sprintf("code %s (%s) par %s\n", fit$opt$convergence, fit$opt$message,
              paste(sprintf("%.10f", fit$opt$par), collapse = " ")))
  g <- fit$obj$gr(fit$opt$par)
  cat("gradient:", format(g, digits = 6), "\n")
  cat("par_units:", format(fit$par_units, digits = 6), "\n")
  cat("se_lost:", paste(names(sdr_of(fit)$se_lost), sdr_of(fit)$se_lost),
      "\n")
  g0 <- glm(y ~ x, family = poisson(), data = d)
  cat(sprintf("poisson glm logLik %.12f  fit logLik %.12f\n",
              as.numeric(logLik(g0)), as.numeric(logLik(fit))))
  p <- fit$opt$par
  for (s in c(14, 16, 18, 20, 22, 23, 24, p[3], 25, 26, 28, 30, 35)) {
    q <- p
    q[3] <- s
    q2 <- c(coef(g0), s)
    cat(sprintf("log shape %.4f  -nll(beta_fit) %.12f  -nll(beta_glm) %.12f  grad3 %.4g\n",
                s, -fit$obj$fn(q), -fit$obj$fn(q2), fit$obj$gr(q)[3]))
  }
}

if (case == "nbk") {
  # k = 0..19 in one process: the ignore fit's code and warnings
  set.seed(2)
  d0 <- data.frame(x = stats::rnorm(200))
  d0$y <- stats::rpois(200, exp(1 + 0.4 * d0$x))
  for (kk in 0:19) {
    d <- d0
    d$x <- d$x * (1 + 2e-15 * kk)
    for (cs in c("warn", "ignore")) {
      r <- conds(frm(y ~ x, family = negbinomial(), data = d,
                     control = frmtmb_control(check_se = cs)))
      fit <- r$value
      cat(sprintf("RESULT nbk k %d %s code %s (%s) shape %.6f nwarn %d nmsg %d\n",
                  kk, cs, fit$opt$convergence, fit$opt$message,
                  fit$opt$par[3], length(r$warnings), length(r$messages)))
      for (x in r$warnings) cat("  W", substr(conditionMessage(x), 1, 200),
                                "\n")
    }
  }
}

if (case == "nb") {
  set.seed(2)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- stats::rpois(200, exp(1 + 0.4 * d$x))
  d$x <- d$x * (1 + 2e-15 * k)
  for (cs in c("warn", "ignore")) {
    r <- conds(frm(y ~ x, family = negbinomial(), data = d,
                   control = frmtmb_control(check_se = cs)))
    fit <- r$value
    cat(sprintf("check_se %s code %s msg %s logLik %.12f\n", cs,
                fit$opt$convergence, fit$opt$message,
                as.numeric(logLik(fit))))
    cat("  par:", format(fit$opt$par, digits = 12), "\n")
    show(r)
    cat("RESULT nb k", k, cs, "nwarn", length(r$warnings), "nmsg",
        length(r$messages), "\n")
  }
}
