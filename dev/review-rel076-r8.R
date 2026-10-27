# Review, short check of R8: posterior_predict() accepts and refuses
# the same calls on 0.76.0 as on 0.75.0 (structured and unstructured
# models), and log_lik(newdata =)'s refusal names posterior_predict()
# only where it answers.
#   Rscript dev/review-rel076-r8.R <lib> > dev/review-rel076-out/r8-<lib>.txt
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb", format(packageVersion("frmtmb")), "sample",
    format(packageVersion("frmtmb.sample")), "\n")
set.seed(7)
d <- expand.grid(t = 1:5, g = 1:12)
d$x <- rnorm(nrow(d))
d$k <- sample(0:3, nrow(d), TRUE)
d$y <- 1 + 0.5 * d$x + 0.3 * d$k + rnorm(12)[d$g] + rnorm(nrow(d))
d$y2 <- 1 - 0.2 * d$x + rnorm(nrow(d))
nd <- d[d$g %in% 1:2, ]
models <- list(
  plain = y ~ x,
  ar_cond = y ~ x + ar(t, g),
  ma_cond = y ~ x + ma(t, g),
  arma_cond = y ~ x + arma(t, g),
  ar_cov = y ~ x + ar(t, g, cov = TRUE),
  ma_cov = y ~ x + ma(t, g, cov = TRUE),
  cosy = y ~ x + cosy(t, g),
  unstr = y ~ x + unstr(t, g),
  mo = y ~ x + mo(k),
  mo_ar_cond = y ~ x + mo(k) + ar(t, g),
  re = y ~ x + (1 | g),
  mv = bf(y ~ x + ar(t, g, cov = TRUE)) + bf(y2 ~ x) + set_rescor(FALSE),
  mv_cond = bf(y ~ x + ar(t, g)) + bf(y2 ~ x) + set_rescor(FALSE))
out <- function(lab, expr) {
  r <- tryCatch({v <- suppressWarnings(expr)
    paste("OK", paste(dim(v), collapse = "x"))},
    error = function(e) paste("ERROR:", substr(gsub("\n", " ",
                                                    conditionMessage(e)),
                                               1, 140)))
  cat(sprintf("%-11s %-34s %s\n", nm, lab, r))
}
for (nm in names(models)) {
  ds <- tryCatch(suppressWarnings(suppressMessages(frm_sample(
    models[[nm]], data = d, chains = 1, iter = 200, seed = 3,
    refresh = 0))), error = function(e) e)
  if (inherits(ds, "error")) {
    cat(sprintf("%-11s frm_sample ERROR: %s\n", nm, conditionMessage(ds)))
    next
  }
  resps <- if (grepl("^mv", nm)) c("y", "y2") else list(NULL)
  for (r in resps) {
    tag <- if (is.null(r)) "" else paste0(" resp=", r)
    out(paste0("pp()", tag), posterior_predict(ds, resp = r, ndraws = 20))
    out(paste0("pp(newdata)", tag), posterior_predict(ds, newdata = nd,
                                                      resp = r, ndraws = 20))
    out(paste0("pp(re_formula = NA)", tag),
        posterior_predict(ds, re_formula = NA, resp = r, ndraws = 20))
    out(paste0("pp(newdata, re_formula = NA)", tag),
        posterior_predict(ds, newdata = nd, re_formula = NA, resp = r,
                          ndraws = 20))
  }
  e <- tryCatch(log_lik(ds, newdata = nd), error = function(e) e)
  msg <- if (inherits(e, "error")) conditionMessage(e) else "RAN"
  cat(sprintf("%-11s %-34s %s\n", nm, "log_lik(newdata) names",
              if (grepl("refused for this model too", msg, fixed = TRUE))
                "epred only (structured branch)" else
                  if (grepl("or posterior_predict(newdata =)", msg,
                            fixed = TRUE)) "epred and predict" else
                    substr(msg, 1, 100)))
}
