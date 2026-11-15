# Review of lane arspat2, attack 1 (last part): fits that do not take the
# new route. Estimates, logLik, SEs, messages and confint_varcorr() saved
# per build; run base (REVLIB="") then lane, and the lane run compares.
source("dev/reviews/arspat2-review-common.R")
tag <- if (nzchar(lib)) "lane" else "base"
n <- 60
M <- 0.5^abs(outer(1:n, 1:n, "-"))
set.seed(4)
d <- data.frame(x = rnorm(n), h = factor(rep(1:6, each = 10)),
                t = rep(1:10, 6))
d$y <- 1 + 0.5 * d$x + rnorm(6)[d$h] +
  as.vector(t(chol(M)) %*% rnorm(n)) * sqrt(4 / rchisq(1, 4))
dm <- d
dm$y[c(3, 11, 25, 40, 52)] <- NA
W <- rook(7)
W <- W / rowSums(W)
set.seed(7)
ds <- data.frame(x = rnorm(49), h = factor(rep(1:7, length.out = 49)))
ds$y <- as.vector(solve(diag(49) - 0.5 * W, 1 + 0.5 * ds$x + rnorm(49)))
dsm <- ds
dsm$y[c(4, 9, 20, 33)] <- NA
pr <- set_prior("gamma(2, 0.1)", class = "nu")
cases <- list(
  gauss_fcor_mi = function() frm(y | mi() ~ x + fcor(M), data = dm,
                                 data2 = list(M = M)),
  gauss_fcor_mi_re = function() frm(y | mi() ~ x + (1 | h) + fcor(M),
                                    data = dm, data2 = list(M = M)),
  stud_fcor_mi_const = function() frm(y | mi() ~ x + fcor(M), data = dm,
                                      data2 = list(M = M), family = student(),
                                      prior = pr),
  stud_fcor_re_nomi = function() frm(y ~ x + (1 | h) + fcor(M), data = d,
                                     data2 = list(M = M), family = student(),
                                     prior = pr),
  gauss_fcor_reml = function() frm(y ~ x + (1 | h) + fcor(M), data = d,
                                   data2 = list(M = M), REML = TRUE),
  stud_sar_mi_const = function() frm(y | mi() ~ x + sar(W), data = dsm,
                                     data2 = list(W = W), family = student(),
                                     prior = pr),
  gauss_sar_re = function() frm(y ~ x + (1 | h) + sar(W, type = "error"),
                                data = ds, data2 = list(W = W)),
  gauss_ar_mi_re = function() frm(y | mi() ~ x + (1 | h) + ar(t, h),
                                  data = dm),
  gauss_ar2_re = function() frm(y ~ x + (1 | h) + ar(t, h, p = 2), data = d)
)
out <- list()
for (nm in names(cases)) {
  msgs <- character(0)
  fit <- withCallingHandlers(cases[[nm]](),
    warning = function(w) {
      msgs <<- c(msgs, paste("W:", conditionMessage(w)))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      msgs <<- c(msgs, paste("M:", conditionMessage(m)))
      invokeRestart("muffleMessage")
    })
  cv <- tryCatch(quiet(confint_varcorr(fit)), error = function(e) NULL)
  s <- quiet(summary(fit))
  out[[nm]] <- list(par = fit$opt$par, ll = as.numeric(logLik(fit)),
                    se = sqrt(diag(vcov(fit))), msgs = msgs, cv = cv,
                    cor_pars = s$cor_pars, conv = fit$opt$convergence)
}
saveRDS(out, sprintf("dev/arspat2-out/ident-%s-fix.rds", tag))
if (tag == "lane") {
  b <- readRDS("dev/reviews/arspat2-review-out/ident-base.rds")
  for (nm in names(out)) {
    a <- out[[nm]]; o <- b[[nm]]
    cat(sprintf("%-20s par %s ll %s (%.3e) se %s msgs %s cv %s cor_pars %s\n", nm,
                identical(unname(a$par), unname(o$par)),
                identical(a$ll, o$ll), a$ll - o$ll,
                identical(unname(a$se), unname(o$se)),
                identical(a$msgs, o$msgs),
                identical(a$cv, o$cv), identical(a$cor_pars, o$cor_pars)))
    if (!identical(a$cv, o$cv)) {
      print(o$cv); print(a$cv)
    }
    if (!identical(a$msgs, o$msgs)) {
      cat("  base msgs:", o$msgs, sep = "\n   "); cat("  lane msgs:", a$msgs, sep = "\n   ")
    }
    if (!identical(unname(a$par), unname(o$par))) {
      cat("  max |par diff|", max(abs(a$par - o$par)), "\n")
    }
  }
}
