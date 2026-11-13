# Review of lane sampfront: subset newdata on frmtmb.spatial's sar() and
# car(), and on log_lik() of laplace draws.
#   Rscript dev/review-sampfront-subset2.R [base]
args <- commandArgs(TRUE)
base <- length(args) && identical(args[1], "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.sample); library(frmtmb.spatial)
})
cat("arm", if (base) "base" else "lane", find.package("frmtmb.sample"), "\n")
full_draws <- function(fit, n = 20L, seed = 1L) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit),
            class = "frmtmb_draws")
}
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
show <- function(lab, expr) {
  cat("\n-- ", lab, "\n", sep = "")
  r <- tryCatch(suppressWarnings(expr), error = function(e) e)
  if (inherits(r, "error")) {
    m <- gsub("\\s+", " ", conditionMessage(r))
    cat("  REFUSED:", if (grepl("This newdata", m))
      sub(".*(This newdata)", "\\1", m) else substr(m, 1, 200), "\n")
  } else {
    cat("  ACCEPTED, dim", dim(r), "\n")
  }
  invisible(r)
}
set.seed(3)
n <- 30
W <- matrix(0, n, n)
for (i in seq_len(n - 1)) W[i, i + 1] <- W[i + 1, i] <- 1
d <- data.frame(x = rnorm(n), loc = factor(seq_len(n)))
d$y <- 0.5 * d$x + rnorm(n)
fs <- tryCatch(q(bf(y ~ x + sar(W, type = "lag")), data = d,
                 data2 = list(W = W)), error = function(e) e)
if (inherits(fs, "error")) {
  cat("sar fit failed:", conditionMessage(fs), "\n")
} else {
  cat("sar frame$autocor:", length(fs$frame$autocor), "\n")
  show("sar(lag): log_lik(newdata = d[1:10, ])",
       log_lik(full_draws(fs), newdata = d[1:10, ]))
}
fse <- tryCatch(q(bf(y ~ x + sar(W, type = "error")), data = d,
                  data2 = list(W = W)), error = function(e) e)
if (!inherits(fse, "error")) {
  show("sar(error): log_lik(newdata = d[1:10, ])",
       log_lik(full_draws(fse), newdata = d[1:10, ]))
}
fc <- tryCatch(q(bf(y ~ x + car(W, gr = loc)), data = d,
                 data2 = list(W = W)), error = function(e) e)
if (inherits(fc, "error")) {
  cat("car fit failed:", conditionMessage(fc), "\n")
} else {
  show("car(W, gr = loc): log_lik(newdata = d[1:10, ])",
       log_lik(full_draws(fc), newdata = d[1:10, ]))
}
cat("\nDONE\n")
