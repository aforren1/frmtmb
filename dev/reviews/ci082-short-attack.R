# Short check of lane ci082's fix round: attacks on varcorr_fixed_rows().
#   Rscript dev/reviews/ci082-short-attack.R <lib|base> <case>
# lib: lane library before rellib-r21; base: rellib-r21 only (0.82.0).
# Cases:
#   usclamp  us() d = 2 at and past the 0.9999 clamp: the probe directly
#            (lib only), and a fit whose slopes are 2x the intercepts
#   oufar    ou() on independent data with one far time (lag >= 960),
#            where the far correlations underflow at every probe point
#   fixclamp a registered structure FIXING its correlation at 0.99995
#   invpath  a registered structure fixing an AR(1) correlation matrix
#            (0.995, d = 6) with 6 estimated sds, its vcov read back
#            through two matrix inverses
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
ns <- asNamespace("frmtmb")
has_probe <- exists("varcorr_fixed_rows", envir = ns, inherits = FALSE)

conds <- function(expr) {
  w <- list()
  v <- withCallingHandlers(expr, warning = function(x) {
    w[[length(w) + 1L]] <<- x
    invokeRestart("muffleWarning")
  }, message = function(x) invokeRestart("muffleMessage"))
  list(value = v, warnings = w)
}
report <- function(tag, fit) {
  cat("==", tag, "code", fit$opt$convergence, fit$opt$message, "\n")
  cat("theta:", format(fit$estimates[["theta"]], digits = 8), "\n")
  tr <- ns$varcorr_trans_rows(fit)
  cv <- conds(confint_varcorr(fit))
  ci <- cv$value
  fx <- if (has_probe) {
    unlist(lapply(fit$frame[["re_blocks"]], function(bk) {
      ns$varcorr_fixed_rows(bk, fit$estimates[["theta"]][bk[["theta_idx"]]])
    }))
  } else rep(NA, nrow(tr))
  out <- data.frame(term = tr$term, type = tr$type,
                    est = signif(ci$estimate, 8), se_t = signif(tr$se_t, 4),
                    width = signif(ci$upr - ci$lwr, 4), fixed = fx)
  print(out, row.names = FALSE)
  for (x in cv$warnings) cat("  confint W:", substr(conditionMessage(x), 1,
                                                    160), "\n")
  cr <- tr$type == "cor"
  cat(sprintf(paste("RESULT %s %s cor rows %d fixed %d se0 %d na %d",
                    "pos %d zero-width %d est0 %d est0-na %d nwarn %d\n"),
              arm, tag, sum(cr), sum(fx[cr] %in% TRUE),
              sum(tr$se_t[cr] == 0, na.rm = TRUE), sum(is.na(tr$se_t[cr])),
              sum(tr$se_t[cr] > 0, na.rm = TRUE),
              sum((ci$upr - ci$lwr)[cr] == 0, na.rm = TRUE),
              sum(ci$estimate[cr] == 0),
              sum(is.na(tr$se_t[cr]) & ci$estimate[cr] == 0),
              length(cv$warnings)))
}

if (case == "usclamp") {
  set.seed(3)
  G <- 60L
  d <- data.frame(g = factor(rep(seq_len(G), each = 8L)),
                  x = rnorm(8L * G))
  b0 <- rnorm(G)
  d$y <- 1 + b0[d$g] + 2 * b0[d$g] * d$x + rnorm(nrow(d), 0, 0.3)
  fit <- conds(frm(y ~ x + (1 + x | g), family = gaussian(), data = d))$value
  report("us perfectly correlated", fit)
  if (has_probe) {
    bk <- fit$frame[["re_blocks"]][[1]]
    for (L in c(70.7, 1e3, 1e5, 1e7, 6.7e7, 1e8, 1e9)) {
      t0 <- c(0, 0, L)
      r <- L / sqrt(1 + L^2)
      cat(sprintf("PROBE us L %.3g cor %.17g fixed %s\n", L, r,
                  paste(ns$varcorr_fixed_rows(bk, t0), collapse = ",")))
    }
  }
}

if (case == "oufar") {
  for (s in 1:3) {
    set.seed(s)
    G <- 50L
    times <- c(0, 0.5, 1, 2, 8, 20, 40, 1000)
    nt <- length(times)
    b <- rnorm(G * nt)
    d <- data.frame(y = rep(b, each = 2L) + rnorm(2L * G * nt, 0, 0.5),
                    g = factor(rep(seq_len(G), each = 2L * nt)),
                    tim = num_factor(rep(rep(times, each = 2L), G)))
    fit <- conds(frm(y ~ 1 + ou(tim + 0 | g), family = gaussian(),
                     data = d))$value
    report(paste0("oufar seed ", s), fit)
  }
}

reg_fixed <- function(name, R) {
  frmtmb_register_covstruct(
    name,
    lpdf = function(b, theta, blk) {
      dd <- blk[["dim"]]
      s <- exp(theta[seq_len(dd)])
      Sigma <- R * (RTMB::matrix(s, ncol = 1) %*% RTMB::matrix(s, nrow = 1))
      dim(b) <- c(dd, length(b) %/% dd)
      sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
    },
    vcov = function(theta, blk) {
      dd <- blk[["dim"]]
      s <- exp(theta[seq_len(dd)])
      # the covariance read back through its precision, as a structure
      # parameterized by a precision matrix would
      Q <- diag(1 / s, dd) %*% solve(R) %*% diag(1 / s, dd)
      V <- solve(Q)
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    npar = nrow(R), sd_idx = seq_len(nrow(R)))
}

if (case == "fixclamp") {
  rr <- 0.99995
  R <- matrix(rr, 2, 2) + diag(1 - rr, 2)
  frmtmb_register_covstruct(
    "xcl",
    lpdf = function(b, theta, blk) {
      Sigma <- exp(2 * theta[1]) * R
      dim(b) <- c(2L, length(b) %/% 2L)
      sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
    },
    vcov = function(theta, blk) {
      V <- exp(2 * theta[1]) * R
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    npar = 1L, sd_idx = 1L)
  set.seed(11L)
  G <- 30L
  d <- expand.grid(t = factor(c("a", "b")), rep = 1:3,
                   g = factor(sprintf("g%02d", seq_len(G))))
  u <- t(crossprod(chol(0.8^2 * R), matrix(rnorm(2 * G), 2)))
  d$y <- 0.4 + u[cbind(as.integer(d$g), as.integer(d$t))] +
    rnorm(nrow(d), 0, 0.6)
  fit <- conds(frm(y ~ 1 + xcl(0 + t | g), family = gaussian(),
                   data = d))$value
  report("fixed at 0.99995", fit)
}

if (case %in% c("invpath", "invgauss")) {
  # invgauss: a gaussian kernel exp(-lag^2 / 2) over 10 levels, whose far
  # correlations (down to 2.6e-18) the double inverse returns with an
  # absolute error near eps, so a relative one near 1
  dd <- if (case == "invpath") 6L else 10L
  lag <- abs(outer(seq_len(dd), seq_len(dd), "-"))
  R <- if (case == "invpath") 0.995^lag else exp(-lag^2 / 2)
  cat("kappa(R):", format(kappa(R, exact = TRUE), digits = 4), "\n")
  reg_fixed("xinv", R)
  sds <- exp(seq(-3, 7, length.out = dd))
  set.seed(4)
  G <- 80L
  tl <- sprintf("t%d", seq_len(dd))
  d0 <- expand.grid(t = factor(tl, levels = tl), rep = 1:2,
                    g = factor(sprintf("g%02d", seq_len(G))))
  U <- matrix(rnorm(G * dd), G) %*% chol(R * (sds %o% sds))
  d0$y <- U[cbind(as.integer(d0$g), as.integer(d0$t))] +
    rnorm(nrow(d0), 0, 0.05)
  for (k in 0:9) {
    d <- d0
    d$y <- d0$y * (1 + 2e-15 * k)
    fit <- conds(frm(y ~ 1 + xinv(0 + t | g), family = gaussian(),
                     data = d))$value
    tr <- ns$varcorr_trans_rows(fit)
    cv <- conds(confint_varcorr(fit))
    cr <- tr$type == "cor"
    w <- (cv$value$upr - cv$value$lwr)[cr]
    bk <- fit$frame[["re_blocks"]][[1]]
    t0 <- fit$estimates[["theta"]][bk[["theta_idx"]]]
    fx <- if (has_probe) ns$varcorr_fixed_rows(bk, t0)[-seq_len(dd)] else NA
    # largest relative move of a correlation row over the probe points
    rows <- function(tt) {
      V <- ns$cs_entry("xinv")$vcov(tt, bk)
      C <- stats::cov2cor(V)
      atanh(C[lower.tri(C)])
    }
    g0 <- rows(t0)
    mv <- 0
    for (i in seq_along(t0)) for (s in c(-1, 1)) {
      tt <- t0
      tt[i] <- tt[i] + s
      gp <- rows(tt)
      mv <- max(mv, abs(gp - g0) / pmax(abs(g0), abs(gp)))
    }
    cat(sprintf(paste("RESULT %s %s k %d code %s cor %d fixed %d se0 %d",
                      "na %d pos %d zero-width %d maxwidth %.3g",
                      "max-rel-move %.3g nwarn %d\n"),
                arm, case, k, fit$opt$convergence, sum(cr), sum(fx %in% TRUE),
                sum(tr$se_t[cr] == 0, na.rm = TRUE), sum(is.na(tr$se_t[cr])),
                sum(tr$se_t[cr] > 0, na.rm = TRUE), sum(w == 0, na.rm = TRUE),
                max(w, na.rm = TRUE), mv, length(cv$warnings)))
    for (x in cv$warnings) cat("  confint W:", substr(conditionMessage(x), 1,
                                                      160), "\n")
  }
}
