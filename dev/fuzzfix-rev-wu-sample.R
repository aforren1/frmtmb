# Reviewer: frmtmb.sample on the weighted observation-level models
# (data of dev/fuzzfix-rev-wu.R, seed 20261021). The formula route skips
# check_weights_unbounded() (dry_run = "objective").
#   Rscript dev/fuzzfix-rev-wu-sample.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "| sample", format(packageVersion("frmtmb.sample")), "\n")
set.seed(20261021)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n), xs = rnorm(n),
                w = runif(n, 0.5, 2), g = factor(rep(1:20, each = 4)))
eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$y <- eta + rnorm(n, 0, exp(0.3 * d$xs))
d$w7 <- d$w * 0.7
run <- function(lbl, expr) {
  ws <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-26s ERROR %s\n", lbl, substr(conditionMessage(r), 1, 200)))
    return(invisible())
  }
  s <- tryCatch(summary(r), error = function(e) NULL)
  dm <- as.matrix(r)
  cn <- grep("sigma|betad", colnames(dm), value = TRUE)[1]
  q <- stats::quantile(dm[, cn], c(0.05, 0.5, 0.95))
  rh <- tryCatch(max(posterior::summarise_draws(r)$rhat, na.rm = TRUE),
                 error = function(e) NA)
  cat(sprintf("%-26s ok %s q05/50/95 %.3g %.3g %.3g max rhat %.3f warnings %d %s\n",
              lbl, cn, q[1], q[2], q[3], rh, length(ws),
              if (length(ws)) substr(ws[1], 1, 80) else ""))
}
G <- gaussian()
fo <- bf(y | weights(w) ~ x + (1 | obs))
run("sample sum(w) > n", frm_sample(fo, data = d, family = G, chains = 2,
                                     iter = 1000, seed = 1, refresh = 0))
run("sample sum(w) < n", frm_sample(bf(y | weights(w7) ~ x + (1 | obs)),
                                     data = d, family = G, chains = 2,
                                     iter = 1000, seed = 1, refresh = 0))
run("fit route, ML fit first", frm_sample(frm(fo, family = G, data = d),
                                           chains = 2, iter = 1000,
                                           seed = 1, refresh = 0))
