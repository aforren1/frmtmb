# Probe C: the RUNG-1 test. The probe B1 model, expressed through
# frm() + custom_family() + vint().
#
# Questions:
#   - does the custom-family contract permit a non-factorizing
#     likelihood with per-row dpars?
#   - does check_custom_family() pass?
#   - does the fit equal probe B1 exactly?
#   - which post-fit methods stay meaningful?
#
# Run: FRMTMB_LIB=<scratch lib> Rscript dev/hmm/probeC1-frm-custom-family.R

lib <- Sys.getenv("FRMTMB_LIB", unset = "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
  library(RTMB)
  library(frmtmb)
})
source("dev/hmm/hmm-common.R")
source("dev/hmm/hmm-family.R")

cat("frmtmb", format(packageVersion("frmtmb")), "\n\n")

ref <- readRDS("dev/hmm/probeB1.rds")
dat <- ref$dat
K <- ref$K

cat("== probe C1: rung 1, HMM through frm() ==\n\n")

## ---- 1. check_custom_family ------------------------------------------

fam <- hmm2_family()
n <- nrow(dat)
test_dpars <- list(mu = rep(0, n), mu2 = rep(3, n),
                   sigma1 = rep(0.7, n), sigma2 = rep(0.7, n),
                   tr12 = rep(-2, n), tr22 = rep(1.5, n))
test_aterms <- list(vint1 = dat$g, vint2 = dat$t)
res <- tryCatch(check_custom_family(fam, y = dat$y, dpars = test_dpars,
                                    aterms = test_aterms),
                error = function(e) conditionMessage(e))
cat("1. check_custom_family: ",
    if (isTRUE(res)) "PASS" else paste("FAIL -", res), "\n\n")

## ---- 2. the fit ------------------------------------------------------

form <- bf(y | vint(g, t) ~ 1, mu2 ~ 1, sigma1 ~ 1, sigma2 ~ 1,
           tr12 ~ x, tr22 ~ x)
t_fit <- system.time(
  fit <- try(frm(form + fam, data = dat), silent = TRUE)
)[["elapsed"]]
if (inherits(fit, "try-error")) {
  cat("2. frm() FAILED:\n", attr(fit, "condition")$message, "\n")
  quit(status = 0)
}
cat(sprintf("2. frm() fit in %.2f s\n", t_fit))
cat("   logLik      :", format(as.numeric(logLik(fit)), digits = 12),
    "\n")
cat("   probe B1 ML :", format(ref$ll, digits = 12), "\n")
cat("   |diff|      :",
    format(abs(as.numeric(logLik(fit)) - ref$ll), digits = 3), "\n")
cat("   df          :", attr(logLik(fit), "df"), " AIC ",
    format(AIC(fit), digits = 10), "\n\n")

cat("3. coefficients vs the hand-rolled MakeADFun optimum\n")
fx <- fixef(fit)
pe <- ref$par
cmp <- rbind(
  frm = c(mu1 = fx$mu[[1]], mu2 = fx$mu2[[1]],
          lsigma1 = fx$sigma1[[1]], lsigma2 = fx$sigma2[[1]],
          tr12_int = fx$tr12[[1]], tr12_x = fx$tr12[[2]],
          tr22_int = fx$tr22[[1]], tr22_x = fx$tr22[[2]]),
  rtmb = c(pe[1], pe[2], pe[3], pe[4], pe[5], pe[7], pe[6], pe[8])
)
print(round(cmp, 8))
cat("   max |diff| :", format(max(abs(cmp[1, ] - cmp[2, ])), digits = 3),
    "\n")
ep <- fit$estimates[["hmm_ldel"]]
cat("   hmm_ldel   :", format(ep, digits = 8), " (rtmb ",
    format(pe[9], digits = 8), ")\n\n")

## ---- 4. sdreport / confint -------------------------------------------

cat("4. inference surface\n")
sm <- try(summary(fit), silent = TRUE)
cat("   summary()  :", if (inherits(sm, "try-error")) "ERROR" else "ok",
    "\n")
ci <- try(confint(fit, method = "wald"), silent = TRUE)
if (inherits(ci, "try-error")) {
  cat("   confint()  : ERROR -", attr(ci, "condition")$message, "\n")
} else {
  cat("   confint(wald) rows:", nrow(ci), "\n")
  print(round(as.data.frame(ci), 5))
}
cat("\n")

## ---- 5. post-fit methods ---------------------------------------------

cat("5. post-fit method survey (family has no post/sim)\n")
probe_method <- function(label, expr) {
  v <- tryCatch(expr, error = function(e) e, warning = function(w) w)
  if (inherits(v, "error")) {
    cat(sprintf("   %-28s REFUSED: %s\n", label,
                substr(conditionMessage(v), 1, 90)))
  } else if (inherits(v, "warning")) {
    cat(sprintf("   %-28s WARNS: %s\n", label,
                substr(conditionMessage(v), 1, 90)))
  } else {
    cat(sprintf("   %-28s returns %s, first 3: %s\n", label,
                paste(class(v), collapse = "/"),
                paste(format(utils::head(as.numeric(unlist(v)), 3),
                             digits = 4), collapse = " ")))
  }
  invisible(v)
}
probe_method("fitted()", fitted(fit))
probe_method("residuals()", residuals(fit))
probe_method("predict(type='link')", predict(fit, type = "link"))
probe_method("predict(type='response')", predict(fit, type = "response"))
probe_method("predict(dpar='tr12')", predict(fit, dpar = "tr12"))
probe_method("simulate()", simulate(fit, nsim = 1, seed = 1))
probe_method("ranef()", ranef(fit))
probe_method("VarCorr()", VarCorr(fit))
cat("\n")

## ---- 6. a family WITH a post-processing block ------------------------

# The honest question is not whether fitted() runs but whether the
# number means anything. E[y_t] under the HMM is
# sum_k P(S_t = k | theta) mu_k with the STATE-OCCUPANCY probability,
# which no per-row mean_fn can see. Supplying dpars$mu makes fitted()
# run and report state 1's mean at every row.
fam2 <- hmm2_family(
  post = list(mean_fn = function(dpars, aterms) dpars$mu,
              var_fn = function(dpars, aterms) dpars$sigma1^2)
)
fit2 <- frm(form + fam2, data = dat)
cat("6. with mean_fn = dpars$mu\n")
cat("   logLik identical  :",
    isTRUE(all.equal(as.numeric(logLik(fit2)),
                     as.numeric(logLik(fit)), tolerance = 1e-10)), "\n")
ft <- fitted(fit2)
cat("   fitted() range    :", format(range(ft), digits = 5),
    " (state 1 mean; y range", format(range(dat$y), digits = 4), ")\n")
cat("   cor(fitted, y)    :", format(cor(ft, dat$y), digits = 4),
    " <- nonsense, as expected\n\n")

## ---- 7. forward-backward post-processing off the fit -----------------

# What a real HMM user wants instead: smoothed state probabilities.
# Computable from the fit numerically, no tape involved.
fb <- function(lpmat, Gof, rows, delta) {
  Tl <- length(rows)
  Kk <- ncol(lpmat)
  la <- matrix(0, Tl, Kk)
  lb <- matrix(0, Tl, Kk)
  a <- delta * exp(lpmat[rows[1], ]); s <- sum(a); lsc <- log(s)
  la[1, ] <- log(a / s) + lsc
  for (k in 2:Tl) {
    a <- as.vector((a / s) %*% Gof(rows[k - 1L])) * exp(lpmat[rows[k], ])
    s <- sum(a); lsc <- lsc + log(s)
    la[k, ] <- log(a / s) + lsc
  }
  lb[Tl, ] <- 0
  for (k in (Tl - 1):1) {
    G <- Gof(rows[k])
    v <- G %*% (exp(lpmat[rows[k + 1L], ] + lb[k + 1L, ]))
    m <- max(log(v))
    lb[k, ] <- log(v) - m + m
  }
  lg <- la + lb
  ex <- exp(lg - apply(lg, 1, max))
  ex / rowSums(ex)
}

pe2 <- ref$par
B_e <- matrix(0, K, K);  B_e[1, 2] <- pe2[5];  B_e[2, 2] <- pe2[6]
Bx_e <- matrix(0, K, K); Bx_e[1, 2] <- pe2[7]; Bx_e[2, 2] <- pe2[8]
lpm <- lpmat_gauss(dat$y, pe2[1:2], exp(pe2[3:4]))
Gof <- function(r) tpm_tv_num(B_e, Bx_e, dat$x[r], K)
rbg <- split(seq_len(nrow(dat)), dat$g)
sp <- do.call(rbind, lapply(rbg, function(r)
  fb(lpm, Gof, r, softmax0(pe2[9]))))
ord <- unlist(rbg)
post_p <- matrix(0, nrow(dat), K)
post_p[ord, ] <- sp
viterbi_state <- max.col(post_p)
cat("7. forward-backward off the fit (numeric, post hoc)\n")
cat("   local-decoding accuracy vs the simulated states:",
    format(mean(viterbi_state == dat$state), digits = 4), "\n")
cat("   E[y] = sum_k P(S=k) mu_k, cor with y:",
    format(cor(post_p %*% pe2[1:2], dat$y), digits = 4), "\n")
cat("   cost: O(n K^2) numeric, ",
    format(system.time(for (i in 1:20) lapply(rbg, function(r)
      fb(lpm, Gof, r, softmax0(pe2[9]))))[["elapsed"]] / 20 * 1000,
      digits = 3), "ms for n =", nrow(dat), "\n\n")

saveRDS(list(ll = as.numeric(logLik(fit)), coef = cmp),
        "dev/hmm/probeC1.rds")
cat("saved dev/hmm/probeC1.rds\n")
