# Reviewer of lane optstart, item R1: brms_missings' fit_imp2 and
# fit_imp3 against a closed-form marginal likelihood written here.
#   Rscript dev/optstart-rev-r1ref.R <lib>      (cwd: the lane worktree)
# Model: chl ~ N(a0 + a1 age, sc); bmi | chl ~ N(b0 + b1 age + (b2 + b3
# age) chl, sb). bmi is linear in chl, so each row's (chl, bmi) is
# bivariate normal and the latent missing values integrate out exactly;
# fit_imp3 adds chl_obs ~ N(chl, se) on the observed rows. No frmtmb code
# is involved in the reference: base R densities, optim(), and
# numDeriv::hessian() of the reference objective (function values only).
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
source("dev/brms-port/shim.R")
data("nhanes", package = "mice")
d <- nhanes
set.seed(port_seed("brms_missings.10.1"))
d$se <- rexp(nrow(d), 2)
ref_nll <- function(th, me) {
  b <- th[1:4]; sb <- exp(th[5]); a <- th[6:7]; sc <- exp(th[8])
  age <- d$age
  mc <- a[1] + a[2] * age
  cc <- b[3] + b[4] * age
  mb <- b[1] + b[2] * age + cc * mc
  ve <- if (me) d$se^2 else rep(0, nrow(d))
  vC <- sc^2 + ve                   # var of the observed chl
  vB <- sb^2 + cc^2 * sc^2          # var of bmi, chl integrated out
  cv <- cc * sc^2                   # their covariance
  oc <- !is.na(d$chl); ob <- !is.na(d$bmi)
  ll <- 0
  i <- oc & ob
  if (any(i)) {
    # bmi | chl_obs, then chl_obs
    m <- mb[i] + cv[i] / vC[i] * (d$chl[i] - mc[i])
    v <- vB[i] - cv[i]^2 / vC[i]
    ll <- ll + sum(dnorm(d$chl[i], mc[i], sqrt(vC[i]), log = TRUE)) +
      sum(dnorm(d$bmi[i], m, sqrt(v), log = TRUE))
  }
  i <- oc & !ob
  ll <- ll + sum(dnorm(d$chl[i], mc[i], sqrt(vC[i]), log = TRUE))
  i <- !oc & ob
  ll <- ll + sum(dnorm(d$bmi[i], mb[i], sqrt(vB[i]), log = TRUE))
  -ll
}
forms <- list(
  fit_imp2 = bf(bmi | mi() ~ age * mi(chl)) + bf(chl | mi() ~ age) +
    set_rescor(FALSE),
  fit_imp3 = bf(bmi | mi() ~ age * mi(chl)) + bf(chl | mi(se) ~ age) +
    set_rescor(FALSE))
for (nm in names(forms)) {
  me <- nm == "fit_imp3"
  w <- character()
  f <- withCallingHandlers(suppressMessages(frm(forms[[nm]], data = d)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  tp <- par_template(f)
  cat(sprintf("\n%s: code %d (%s), %d warnings, logLik %.10f\n", nm,
              f$opt$convergence, f$opt$message, length(w),
              as.numeric(logLik(f))))
  if (length(w)) cat("  warning:", substr(w, 1, 150), sep = "\n  ")
  print(names(tp$beta)); print(names(tp$betad))
  # reference parameter order: bmi b (4), log sb, chl a (2), log sc
  est <- c(tp$beta[1:4], tp$betad[1], tp$beta[5:6], tp$betad[2])
  best <- NULL
  for (k in 1:4) {
    st <- if (k == 1) est else est + rnorm(8, 0, 0.05)
    o <- optim(st, ref_nll, me = me, method = "BFGS",
               control = list(maxit = 5000, reltol = 1e-15))
    o <- optim(o$par, ref_nll, me = me, method = "Nelder-Mead",
               control = list(maxit = 20000, reltol = 1e-15))
    o <- optim(o$par, ref_nll, me = me, method = "BFGS",
               control = list(maxit = 5000, reltol = 1e-15))
    if (is.null(best) || o$value < best$value) best <- o
  }
  cat(sprintf("reference: logLik at frmtmb's estimate %.10f, max %.10f\n",
              -ref_nll(est, me), -best$value))
  Hr <- numDeriv::hessian(function(t) ref_nll(t, me), best$par)
  se_ref <- sqrt(diag(solve(Hr)))
  # frmtmb's route: the outer Hessian it inverts for vcov()
  se_v <- suppressWarnings(sqrt(diag(vcov(f))))
  tn <- gsub("[()]", "", names(tp$beta))
  se_v <- se_v[tn]
  Hf <- frmtmb:::fit_outer_hessian(f)$H
  se_h <- sqrt(diag(solve(Hf)))
  nb <- length(tp$beta)
  se_h <- c(se_h[1:4], se_h[nb + 1], se_h[5:6], se_h[nb + 2])
  se_vv <- c(se_v[1:4], NA, se_v[5:6], NA)
  out <- data.frame(par = c(names(tp$beta)[1:4], names(tp$betad)[1],
                            names(tp$beta)[5:6], names(tp$betad)[2]),
                    est = signif(est, 7), ref_est = signif(best$par, 7),
                    se_vcov = signif(se_vv, 6), se_outerH = signif(se_h, 6),
                    se_ref = signif(se_ref, 6),
                    rel_vcov = signif(se_vv / se_ref - 1, 3),
                    rel_outerH = signif(se_h / se_ref - 1, 3))
  print(out, row.names = FALSE)
  cat(sprintf("max |rel| vcov %.3g, outer Hessian %.3g\n",
              max(abs(out$rel_vcov), na.rm = TRUE),
              max(abs(out$rel_outerH))))
}
