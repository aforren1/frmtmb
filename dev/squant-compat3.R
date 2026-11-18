# Lane squant, fix round 1: fits behind three compat rows of
# quantile_elf: ar() through the latent residual, an nl formula, and a
# mixture of two ELFs (seed 7). Run: Rscript dev/squant-compat3.R
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
out <- function(lab, expr) {
  ws <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, substr(conditionMessage(w), 1, 100))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  if (inherits(r, "error")) {
    cat(lab, ": ERROR", substr(conditionMessage(r), 1, 200), "\n")
  } else {
    cat(lab, ": conv", r$opt$convergence, "logLik", as.numeric(logLik(r)),
        "warnings", length(ws), paste(ws, collapse = " | "), "\n")
    print(fixef(r)[, 1:2])
  }
}
set.seed(7)
n <- 200
d <- data.frame(x = runif(n), t = seq_len(n), g = factor(rep(1:10, each = 20)))
d$y <- 1 + 2 * d$x + as.numeric(stats::filter(rnorm(n), 0.5, "recursive"))
out("ar", frm(bf(y ~ x + ar(time = t, gr = g), quantile = 0.5),
              family = quantile_elf(), data = d))
d$yn <- 3 * exp(-0.8 * d$x) + rnorm(n, 0, 0.3)
out("nl", frm(bf(yn ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE,
                 quantile = 0.5), family = quantile_elf(), data = d,
              prior = c(prior(normal(3, 2), nlpar = "a"),
                        prior(normal(1, 1), nlpar = "b"))))
z <- rbinom(n, 1, 0.4)
d$ym <- ifelse(z == 1, -2 + rnorm(n, 0, 0.5), 2 + rnorm(n, 0, 0.5))
out("mixture", frm(bf(ym ~ 1, mu1 ~ 1, mu2 ~ 1, quantile1 = 0.5,
                      quantile2 = 0.5),
                   family = mixture(quantile_elf(), quantile_elf()),
                   data = d))
