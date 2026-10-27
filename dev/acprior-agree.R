# The posterior of an autocorrelation parameter near its bound under
# frm_sample()'s default priors, before and after this lane, against
# brms 2.23.0 under its own defaults, on the same data.
#
# Usage: Rscript dev/acprior-agree.R <arm> where arm is "before"
# (rellib-r13), "after" (the lane library) or "brms". Each arm saves
# its draws to dev/acprior-out/agree-<arm>.rds; "summary" reads all
# three and prints the comparison (log: dev/acprior-out/agree.txt).
#
# Data (seed 20261027):
#   ar:   y = 1 + 0.5 x + e, e an AR(1) with phi = 0.9 within each of
#         8 groups of 6 time points; fitted with ar(t, g, cov = TRUE).
#   cosy: y = 1 + 0.5 x + u_g + e, cosy correlation 0.03 (u_g ~
#         N(0, 0.03), e ~ N(0, 0.97)), 30 groups of 4; cosy(t, g).
# brms parameterizes cov = TRUE ar() by the innovation sd (sigma_i) and
# frmtmb by the marginal one (sigma_m = sigma_i / sqrt(1 - phi^2)), so
# the same student_t(3, 0, 2.5) default on sigma is a different joint
# prior. brms's ar draws are therefore importance-weighted to frmtmb's
# model: w = pt(sigma_i / s) / (pt(sigma_i) s), s = sqrt(1 - phi^2),
# which is the density ratio after the change of variables sigma_i ->
# sigma_m. cosy uses the marginal sigma on both sides; no weights.
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
out <- file.path(wt, "dev", "acprior-out")
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
lib <- c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-acprior-lib",
         "C:/Users/adf44/source/r/rellib-r13", ulib)
.libPaths(lib)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
Sys.setenv(FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))

make_data <- function() {
  set.seed(20261027)
  ng <- 8; k <- 6
  da <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                   x = rnorm(ng * k))
  e <- numeric(nrow(da))
  for (i in seq_len(nrow(da))) {
    e[i] <- if (da$t[i] == 1) rnorm(1) else 0.9 * e[i - 1] +
      rnorm(1, 0, sqrt(1 - 0.81))
  }
  da$y <- 1 + 0.5 * da$x + e
  ng2 <- 30; k2 <- 4
  dc <- data.frame(t = rep(seq_len(k2), ng2),
                   g = rep(seq_len(ng2), each = k2), x = rnorm(ng2 * k2))
  dc$y <- 1 + 0.5 * dc$x + rnorm(ng2, 0, sqrt(0.03))[dc$g] +
    rnorm(nrow(dc), 0, sqrt(0.97))
  list(ar = da, cosy = dc)
}
dat <- make_data()
fo <- list(ar = y ~ x + ar(time = t, gr = g, cov = TRUE),
           cosy = y ~ x + cosy(time = t, gr = g))
CH <- 4; IT <- 4000; SEED <- 3

if (arm %in% c("before", "after")) {
  suppressPackageStartupMessages({
    library(frmtmb)
    library(frmtmb.sample)
  })
  cat("frmtmb", format(packageVersion("frmtmb")), "from",
      dirname(find.package("frmtmb")), "; frmtmb.sample",
      format(packageVersion("frmtmb.sample")), "from",
      dirname(find.package("frmtmb.sample")), "\n")
  res <- list()
  for (m in names(fo)) {
    fit <- frm(fo[[m]], data = dat[[m]])
    msgs <- character(0)
    ds <- withCallingHandlers(
      frm_sample(fit, chains = CH, iter = IT, seed = SEED, refresh = 0),
      message = function(c) {
        msgs <<- c(msgs, conditionMessage(c))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        msgs <<- c(msgs, paste("WARNING:", conditionMessage(w)))
        invokeRestart("muffleWarning")
      })
    cat("\n--", m, "messages:\n", paste(msgs, collapse = "\n"), "\n")
    mm <- as.matrix(ds)
    col <- if (m == "ar") "ar[1]" else "cosy"
    res[[m]] <- list(par = mm[, col], chains = CH,
                     ml = unname(frmtmb:::autocor_natural(
                       fit$estimates$thetaac, fit$frame$autocor[[1]])))
  }
  saveRDS(res, file.path(out, paste0("agree-", arm, ".rds")))
} else if (identical(arm, "brms")) {
  suppressPackageStartupMessages(library(brms))
  cat("brms", format(packageVersion("brms")), "\n")
  res <- list()
  for (m in names(fo)) {
    code <- stancode(fo[[m]], data = dat[[m]])
    sdat <- standata(fo[[m]], data = dat[[m]])
    f <- tempfile(fileext = ".stan")
    writeLines(c(code, paste0("// rstan ", packageVersion("rstan"))), f)
    key <- unname(tools::md5sum(f))
    path <- file.path(Sys.getenv("FRMTMB_STAN_CACHE"), paste0(key, ".rds"))
    mod <- if (file.exists(path)) readRDS(path) else {
      mm <- rstan::stan_model(model_code = code, save_dso = TRUE)
      saveRDS(mm, path)
      mm
    }
    sf <- rstan::sampling(mod, data = sdat, chains = CH, iter = IT,
                          seed = SEED, refresh = 0)
    par <- if (m == "ar") "ar[1]" else "cosy"
    a <- as.matrix(sf, pars = c(par, "sigma"))
    print(rstan::summary(sf, pars = c(par, "sigma"))$summary)
    res[[m]] <- list(par = a[, 1], sigma = a[, 2], chains = CH)
  }
  saveRDS(res, file.path(out, "agree-brms.rds"))
} else {
  .libPaths(ulib)
  rd <- function(a) readRDS(file.path(out, paste0("agree-", a, ".rds")))
  B <- rd("before"); A <- rd("after"); R <- rd("brms")
  pt <- function(s) dt(s / 2.5, 3) / 2.5
  wq <- function(x, w, p) {
    o <- order(x); cw <- cumsum(w[o]) / sum(w)
    vapply(p, function(pp) x[o][which(cw >= pp)[1]], 0)
  }
  ess <- function(x, ch) {
    posterior::ess_bulk(matrix(x, ncol = ch))
  }
  for (m in c("ar", "cosy")) {
    cat("\n=====", m, "\n")
    cat(sprintf("ML estimate %.4f\n", A[[m]]$ml))
    rb <- R[[m]]
    w <- rep(1, length(rb$par))
    if (m == "ar") {
      s <- sqrt(1 - rb$par^2)
      w <- pt(rb$sigma / s) / (pt(rb$sigma) * s)
    }
    wn <- w / sum(w)
    ess_w <- 1 / sum(wn^2)
    row <- function(lab, x, wts = NULL, essv) {
      if (is.null(wts)) wts <- rep(1, length(x))
      mu <- sum(wts * x) / sum(wts)
      sdv <- sqrt(sum(wts * (x - mu)^2) / sum(wts))
      q <- wq(x, wts, c(0.05, 0.5, 0.95))
      cat(sprintf("%-26s mean %.4f sd %.4f q05 %.4f q50 %.4f q95 %.4f  ESS %.0f  mcse %.4f\n",
                  lab, mu, sdv, q[1], q[2], q[3], essv, sdv / sqrt(essv)))
      c(mu = mu, se = sdv / sqrt(essv))
    }
    b0 <- row("frm_sample before", B[[m]]$par, essv = ess(B[[m]]$par, 4))
    a0 <- row("frm_sample after", A[[m]]$par, essv = ess(A[[m]]$par, 4))
    r0 <- row("brms", rb$par, essv = ess(rb$par, 4))
    if (m == "ar") {
      r1 <- row("brms, weighted to frmtmb", rb$par, w,
                essv = min(ess_w, ess(rb$par, 4)))
    } else r1 <- r0
    z <- function(a, b) (a[["mu"]] - b[["mu"]]) / sqrt(a[["se"]]^2 + b[["se"]]^2)
    cat(sprintf("z(after - brms ref) %.2f; z(before - brms ref) %.2f\n",
                z(a0, r1), z(b0, r1)))
    if (m == "cosy") {
      cat(sprintf("share of draws below 0: before %.4f after %.4f brms %.4f\n",
                  mean(B[[m]]$par < 0), mean(A[[m]]$par < 0),
                  mean(rb$par < 0)))
    }
  }
}
