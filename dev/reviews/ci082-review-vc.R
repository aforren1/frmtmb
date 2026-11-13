# Review of lane ci082: does the "fixed row" rule of varcorr_trans_rows()
# call an ESTIMATED row fixed (zero-width interval, no warning)?
#   Rscript dev/reviews/ci082-review-vc.R <lib|base> <case>
# lib: lane library before rellib-r21; base: rellib-r21 only.
args <- commandArgs(TRUE)
arm <- args[1]
case <- args[2]
base <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/wt-ci082-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")),
    format(packageVersion("frmtmb")), "\n")
cat("La:", La_library(), " CORETYPE:", Sys.getenv("OPENBLAS_CORETYPE"), "\n")

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
report <- function(tag, r) {
  fit <- r$value
  cat("==", tag, "code", fit$opt$convergence, fit$opt$message, "\n")
  cat("theta:", format(fit$estimates[["theta"]], digits = 10), "\n")
  for (x in r$warnings) cat("  fit W:", substr(conditionMessage(x), 1, 160),
                            "\n")
  for (x in r$messages) cat("  fit M:", substr(conditionMessage(x), 1, 160),
                            "\n")
  tr <- frmtmb:::varcorr_trans_rows(fit)
  cv <- conds(confint_varcorr(fit))
  ci <- cv$value
  out <- data.frame(term = tr$term, type = tr$type,
                    est_t = signif(tr$est_t, 8), se_t = signif(tr$se_t, 4),
                    estimate = signif(ci$estimate, 10),
                    width = signif(ci$upr - ci$lwr, 4))
  print(out, row.names = FALSE)
  for (x in cv$warnings) cat("  confint W:", substr(conditionMessage(x), 1,
                                                    200), "\n")
  cat(sprintf("RESULT %s %s rows %d se0 %d na %d tiny(se<1e-8,>0) %d nwarn %d\n",
              arm, tag, nrow(tr), sum(tr$se_t == 0, na.rm = TRUE),
              sum(is.na(tr$se_t)),
              sum(tr$se_t > 0 & tr$se_t < 1e-8, na.rm = TRUE),
              length(cv$warnings)))
  invisible(tr)
}

if (case == "cslow") {
  # an exchangeable block whose correlation runs to its lower limit
  # -1/(d-1) = -0.5: each group's three effects sum to zero
  for (s in 1:3) {
    set.seed(s)
    G <- 40L
    d <- expand.grid(t = factor(c("a", "b", "c")), rep = 1:4,
                     g = factor(sprintf("g%02d", seq_len(G))))
    z <- matrix(rnorm(3 * G), 3)
    u <- t(sweep(z, 2, colMeans(z)))
    d$y <- 1 + 1.2 * u[cbind(as.integer(d$g), as.integer(d$t))] +
      rnorm(nrow(d), 0, 0.3)
    for (st in c("homcs", "cs")) {
      f <- stats::as.formula(paste0("y ~ 1 + ", st, "(0 + t | g)"))
      report(paste0(st, " seed ", s), conds(frm(f, family = gaussian(),
                                               data = d)))
    }
  }
}

if (case == "ar1") {
  # a weak AR(1) over ten levels: the far lags are rho^k, which the
  # jacobian moves by less than 16 eps
  for (s in 1:3) {
    set.seed(s)
    G <- 60L
    K <- 10L
    rho <- 0.05
    S <- rho^abs(outer(1:K, 1:K, "-"))
    u <- matrix(rnorm(G * K), G) %*% chol(S)
    d <- expand.grid(t = factor(sprintf("t%02d", 1:K)), rep = 1:2,
                     g = factor(sprintf("g%02d", seq_len(G))))
    d$y <- u[cbind(as.integer(d$g), as.integer(d$t))] +
      rnorm(nrow(d), 0, 0.5)
    report(paste0("ar1 seed ", s), conds(frm(y ~ 1 + ar1(0 + t | g),
                                             family = gaussian(), data = d)))
  }
}

if (case == "ou") {
  # far-apart times: exp(-rate * D) is tiny but estimated
  for (s in 1:2) {
    set.seed(s)
    G <- 50L
    times <- c(0, 0.5, 1, 2, 8, 20, 40)
    Sig <- exp(-1.0 * abs(outer(times, times, "-")))
    U <- matrix(rnorm(G * length(times)), G) %*% chol(Sig)
    d <- data.frame(y = as.vector(t(U)) + rnorm(G * length(times), 0, 0.4),
                    g = factor(rep(seq_len(G), each = length(times))),
                    tim = num_factor(rep(times, G)))
    report(paste0("ou seed ", s), conds(frm(y ~ 1 + ou(tim + 0 | g),
                                            family = gaussian(), data = d)))
  }
}

if (case == "ouind") {
  # independent data: the OU rate runs to +Inf (independence, the edge
  # of the space), and the far correlations underflow to exactly 0
  for (s in 1:3) {
    set.seed(s)
    G <- 50L
    times <- c(0, 0.5, 1, 2, 8, 20, 40)
    # two replicates per group and time: b iid over (g, t), no
    # correlation across times
    nt <- length(times)
    b <- rnorm(G * nt)
    d <- data.frame(y = rep(b, each = 2L) + rnorm(2L * G * nt, 0, 0.5),
                    g = factor(rep(seq_len(G), each = 2L * nt)),
                    tim = num_factor(rep(rep(times, each = 2L), G)))
    report(paste0("ouind seed ", s), conds(frm(y ~ 1 + ou(tim + 0 | g),
                                               family = gaussian(),
                                               data = d)))
  }
}

if (case == "csstart") {
  # an exchangeable block started (as a caller may) far out on its
  # logistic correlation, where the objective is flat: the stop is at
  # rho = -1/(d-1), a singular correlation matrix, the edge of the space
  set.seed(1)
  G <- 40L
  d <- expand.grid(t = factor(c("a", "b", "c")), rep = 1:4,
                   g = factor(sprintf("g%02d", seq_len(G))))
  z <- matrix(rnorm(3 * G), 3)
  u <- t(sweep(z, 2, colMeans(z)))
  d$y <- 1 + 1.2 * u[cbind(as.integer(d$g), as.integer(d$t))] +
    rnorm(nrow(d), 0, 0.3)
  for (t2 in c(-20, -26, -28, -30, -34)) {
    report(paste0("homcs start ", t2),
           conds(frm(y ~ 1 + homcs(0 + t | g), family = gaussian(),
                     data = d, start = list(theta = c(0, t2)))))
  }
}

if (case %in% c("fix999", "fix9")) {
  # a FIXED correlation where atanh is steep: is "fixed" still decided
  # against rounding?
  rr <- if (case == "fix999") 0.999 else 0.99
  xc <- function(dd) matrix(rr, dd, dd) + diag(1 - rr, dd)
  frmtmb_register_covstruct(
    "xcs",
    lpdf = function(b, theta, blk) {
      dd <- blk[["dim"]]
      Sigma <- exp(2 * theta[1]) * xc(dd)
      dim(b) <- c(dd, length(b) %/% dd)
      sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
    },
    vcov = function(theta, blk) {
      V <- exp(2 * theta[1]) * xc(blk[["dim"]])
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    npar = 1L, sd_idx = 1L)
  set.seed(11L)
  G <- 30L
  d0 <- expand.grid(t = factor(c("a", "b", "c")), rep = 1:3,
                    g = factor(sprintf("g%02d", seq_len(G))))
  u <- t(crossprod(chol(0.8^2 * xc(3L)), matrix(rnorm(3 * G), 3)))
  d0$x <- rnorm(nrow(d0))
  d0$y <- 0.4 + 0.7 * d0$x + u[cbind(as.integer(d0$g), as.integer(d0$t))] +
    rnorm(nrow(d0), 0, 0.6)
  for (k in 0:19) {
    d <- d0
    d$y <- d0$y * (1 + 2e-15 * k)
    r <- conds(frm(bf(y ~ x + xcs(0 + t | g)), family = gaussian(),
                   data = d))
    tr <- frmtmb:::varcorr_trans_rows(r$value)
    cv <- conds(confint_varcorr(r$value))
    cr <- tr$type == "cor"
    w <- cv$value$upr[cr] - cv$value$lwr[cr]
    cat(sprintf("RESULT %s %s k %d se_cor %s width %s nwarn %d\n", arm, case,
                k, paste(signif(tr$se_t[cr], 3), collapse = ","),
                paste(signif(w, 3), collapse = ","), length(cv$warnings)))
  }
}
