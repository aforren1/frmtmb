# Punch round 1, m9: the log density of a multivariate model with one
# ordinal-mixture response, against brms 2.23.0's compiled program at
# matched parameter values (not by fitting both). Response y is
# mixture(cumulative(), sratio()) with theta1 ~ z, response w is a
# gaussian; brms needs set_rescor(FALSE) for the pair. Every brms prior
# is flat, so its log density with adjust_transform = FALSE is the log
# likelihood plus the constant of its dirichlet(1) on a simplex theta
# (none here: theta1 has a predictor). Points: frmtmb's optimum and
# three perturbed vectors (seed 99).
# Usage: Rscript dev/ordmix-p1-mv.R
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(brms)
  library(rstan)
})
cache_dir <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/stan-cache"
stan_mod <- function(code) {
  f <- tempfile(fileext = ".stan")
  writeLines(c(code, paste0("// rstan ", packageVersion("rstan"))), f)
  key <- unname(tools::md5sum(f))
  path <- file.path(cache_dir, paste0(key, ".rds"))
  if (file.exists(path)) {
    m <- try(readRDS(path), silent = TRUE)
    if (!inherits(m, "try-error")) return(m)
  }
  m <- rstan::stan_model(model_code = code, save_dso = TRUE)
  saveRDS(m, path)
  m
}
set.seed(20261005 + 77)
n <- 400
x <- rnorm(n)
z <- rnorm(n)
cls <- rbinom(n, 1, plogis(-0.4 + 0.8 * z))
lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
w <- 0.5 + 0.7 * x + rnorm(n, 0, 1.3)
d <- data.frame(y, w, x, z)

ff <- frmtmb::bf(y ~ x, theta1 ~ z) + frmtmb::bf(w ~ x)
fit <- frm(ff, family = list(frmtmb::mixture(frmtmb::cumulative(),
                                              frmtmb::sratio()),
                             gaussian()),
           data = d, control = frmtmb_control(grad_tol = 1e-8))
cat("frmtmb logLik", format(as.numeric(logLik(fit)), digits = 15), "\n")
fb <- brms::bf(y ~ x, theta1 ~ z,
               family = brms::mixture(brms::cumulative(), brms::sratio())) +
  brms::bf(w ~ x, family = gaussian()) + brms::set_rescor(FALSE)
pr <- brms::get_prior(fb, data = d)
pr$prior <- ""
code <- as.character(brms::stancode(fb, data = d, prior = pr))
sdat <- brms::standata(fb, data = d, prior = pr)
mod <- stan_mod(code)
sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
pn <- sf@model_pars
cat("brms parameters:", paste(setdiff(pn, "lp__"), collapse = " "), "\n")

## translation, by name --------------------------------------------------
resp <- names(fit$spec$responses)
cat("frmtmb responses:", paste(resp, collapse = " "), "\n")
lpk <- function(r, dp) {
  fit$frame$linpreds[[frmtmb:::linpred_key(fit$spec$responses[[r]]$resp_name,
                                            dp)]]
}
coef_of <- function(est, r, dp) {
  lp <- lpk(r, dp)
  v <- est[[lp$par]][lp$idx]
  names(v) <- colnames(lp$X)[seq_along(v)]
  v
}
ry <- resp[1]
rw <- resp[2]
famy <- fit$spec$responses[[ry]]$family
mx <- famy$mix$ord
# brms centers each design: its Intercept is the uncentered intercept
# plus colMeans(X) %*% b, and an ordinal threshold minus that sum
cm <- function(sx) {
  m <- colMeans(sdat[[paste0("X", sx)]])
  m[setdiff(names(m), "Intercept")]
}
translate <- function(est) {
  out <- list()
  for (k in 1:2) {
    dp <- paste0("mu", k)
    b <- coef_of(est, ry, dp)
    out[[paste0("b_", dp, "_y")]] <- array(unname(b[["x"]]), 1)
    tm <- names(fit$frame$extra_map[[ry]] %||% list())
    raw <- est[[fit$frame$extra_map[[ry]][[mx$tau_names[k]]] %||%
                  mx$tau_names[k]]]
    tau <- frmtmb:::ord_threshold_values(mx$views[[k]], raw)
    out[[paste0("Intercept_", dp, "_y")]] <- tau - sum(cm(paste0("_", dp, "_y")) *
                                                 b[["x"]])
  }
  bt <- coef_of(est, ry, "theta1")
  out[["b_theta1_y"]] <- array(unname(bt[["z"]]), 1)
  out[["Intercept_theta1_y"]] <- unname(bt[["(Intercept)"]]) +
    sum(cm("_theta1_y") * bt[["z"]])
  bw <- coef_of(est, rw, "mu")
  out[["b_w"]] <- array(unname(bw[["x"]]), 1)
  out[["Intercept_w"]] <- unname(bw[["(Intercept)"]]) +
    sum(cm("_w") * bw[["x"]])
  out[["sigma_w"]] <- exp(unname(coef_of(est, rw, "sigma")[["(Intercept)"]]))
  out
}
`%||%` <- function(a, b) if (is.null(a)) b else a
ulp <- function(a, b) abs(a - b) / (.Machine$double.eps * max(abs(a), abs(b)))
lp_at <- function(p) {
  est <- fit$obj$env$parList(p)
  up <- rstan::unconstrain_pars(sf, translate(est))
  list(brms = rstan::log_prob(sf, up, adjust_transform = FALSE),
       frm = -fit$obj$fn(p), up = up)
}
r0 <- lp_at(fit$opt$par)
g0 <- rstan::grad_log_prob(sf, r0$up, adjust_transform = FALSE)
cat(sprintf("LP point=opt brms=%.15g frm=%.15g ulp=%.1f\n", r0$brms, r0$frm,
            ulp(r0$brms, r0$frm)))
cat(sprintf("GRAD max|grad brms at frmtmb opt|=%.3g npar=%d\n",
            max(abs(g0)), length(g0)))
set.seed(99)
for (j in 1:3) {
  p <- fit$opt$par + rnorm(length(fit$opt$par), 0, 0.3)
  r <- lp_at(p)
  cat(sprintf("LP point=perturb%d brms=%.15g frm=%.15g ulp=%.1f\n", j,
              r$brms, r$frm, ulp(r$brms, r$frm)))
}
