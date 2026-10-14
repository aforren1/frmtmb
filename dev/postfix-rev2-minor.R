# Reviewer of lane postfix, round 2: the minor items of fix round 1.
.libPaths(c("C:/Users/adf44/source/r/wt-postfix-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
cnd <- function(expr) {
  out <- character()
  tryCatch(withCallingHandlers(expr,
    warning = function(w) {
      out <<- c(out, paste("[warning]", conditionMessage(w)))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage")),
    error = function(e) out <<- c(out, paste("[error]", conditionMessage(e))))
  out
}
set.seed(1)
dd <- data.frame(x = rnorm(60))
dd$o <- cut(dd$x + rlogis(60), c(-Inf, -1, 0, 1, Inf), ordered_result = TRUE)
for (e in list(quote(frm(bf(o ~ x) + cratio(), data = dd, threshold = "equidistant")),
               quote(frm(bf(o ~ x) + cratio("logit"), data = dd, threshold = "equidistant")),
               quote(frm(bf(o ~ x, family = acat()), data = dd, threshold = "equidistant")),
               quote(frm(bf(o ~ x), family = sratio, data = dd, threshold = "equidistant")),
               quote(frm(bf(o ~ x), family = "acat", data = dd, threshold = "equidistant")),
               quote(frm(bf(o ~ x) + gaussian(), data = dd, threshold = "equidistant")))) {
  cat(deparse1(e), "\n  ", substr(cnd(eval(e)), 1, 400), "\n")
}
## yhash
source("dev/postfix-rev-data.R")
q <- function(e) suppressWarnings(suppressMessages(e))
sp <- q(frm_sample(frm(bf(y ~ x), family = poisson(), data = dp),
                   chains = 2, iter = 600, seed = 1, refresh = 0))
s0 <- q(frm_sample(frm(bf(y ~ 1), family = poisson(), data = dp),
                   chains = 2, iter = 600, seed = 1, refresh = 0))
dp2 <- dp; dp2$y <- rev(dp2$y)
s2 <- q(frm_sample(frm(bf(y ~ x), family = poisson(), data = dp2),
                   chains = 2, iter = 600, seed = 1, refresh = 0))
cat("loo(sp, s2) different responses:\n ",
    unique(grep("Pareto", cnd(loo(sp, s2)), value = TRUE, invert = TRUE)), "\n")
cat("loo(sp, s0) same response:\n ",
    unique(grep("Pareto", cnd(loo(sp, s0)), value = TRUE, invert = TRUE)), "\n")
cat("waic(sp, s2):\n ",
    unique(grep("p_waic", cnd(waic(sp, s2)), value = TRUE, invert = TRUE)), "\n")
cat("loo_compare(sp, s2):\n ",
    unique(grep("Pareto", cnd(loo_compare(sp, s2)), value = TRUE, invert = TRUE)), "\n")
yh <- attr(q(loo(sp)), "yhash")
cat("yhash attr of one loo:", if (is.null(yh)) "NULL" else class(yh), "\n")
## the same loo objects through brms's compare_ic, which reads yhash
two <- q(loo(sp, s2))
cat("brms compare_ic on them:\n ", unique(cnd(brms:::compare_ic(x = two$loos))), "\n")
## transform warning
cat("ce transform:\n ", unique(cnd(conditional_effects(sp, "x", method = "predict",
                                         transform = sqrt, seed = 1))), "\n")
cat("posterior_predict transform:\n ",
    unique(cnd(posterior_predict(sp, transform = sqrt))), "\n")
cat("ce without transform:", length(cnd(conditional_effects(sp, "x", method = "predict",
                                                   seed = 1))), "conditions\n")
