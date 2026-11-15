# Re-check of lane priors2 after fix round 2, item 3: the held warning is
# raised again where the box run's own result is returned. A unit run of
# open_end_hold(), then real fits whose box run starts inside and crosses
# mu < 0 rows on its way (identity links, intercept near 0).
# Rscript priors2-short-r2c.R <lane|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
if (arm == "lane") {
  b <- list(lower = 0, upper = Inf, open_lower = TRUE, open_upper = FALSE)
  for (refit in c(FALSE, TRUE)) {
    h <- frmtmb:::open_end_hold(b, 1L, "nlminb")
    withCallingHandlers(
      warning("NA/NaN function evaluation", call. = FALSE),
      warning = h$catch)
    got <- character(0)
    withCallingHandlers(h$release(list(open_refit = refit)),
                        warning = function(w) {
                          got <<- c(got, paste0(class(w)[1], ": ",
                                                conditionMessage(w)))
                          invokeRestart("muffleWarning")
                        })
    cat("unit: open_refit =", refit, "-> re-raised:",
        if (length(got)) got else "none", "\n")
  }
  h0 <- frmtmb:::open_end_hold(list(lower = 0, upper = Inf), 1L, "nlminb")
  r <- tryCatch(withCallingHandlers(
    warning("NA/NaN function evaluation", call. = FALSE),
    warning = h0$catch), warning = function(w) "passed through")
  cat("unit: no open end -> the warning", if (identical(r, "passed through"))
    "passes through" else "is held", "\n")
}
show <- function(lab, expr) {
  ws <- character(0)
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, paste0(class(w)[1], ": ", substr(conditionMessage(w), 1, 70)))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  if (inherits(r, "error")) {
    cat(lab, ": ERROR ", substr(conditionMessage(r), 1, 120), "\n", sep = "")
  } else {
    cat(lab, ": code ", r$opt$convergence, " hex ",
        paste(sprintf("%a", c(r$opt$objective, r$opt$par)), collapse = ","),
        "\n", sep = "")
  }
  for (w in ws) cat("    [warning]", w, "\n")
}
sp <- set_prior
for (s in 1:6) {
  set.seed(100 + s)
  n <- 120
  d <- data.frame(x = runif(n, 0, 3))
  d$ye <- rexp(n, 1 / (0.03 + 0.8 * d$x))
  d$yp <- rpois(n, 0.02 + 1.2 * d$x)
  for (pr in c("lognormal(0, 1)", "gamma(2, 2)")) {
    show(sprintf("seed %d exponential(identity) %s on x", 100 + s, pr),
         frm(ye ~ x, data = d, family = exponential(link = "identity"),
             prior = sp(pr, class = "b", coef = "x", lb = 0),
             start = list(beta = c(1, 1))))
    show(sprintf("seed %d poisson(identity) %s on x", 100 + s, pr),
         frm(yp ~ x, data = d, family = poisson(link = "identity"),
             prior = sp(pr, class = "b", coef = "x", lb = 0),
             start = list(beta = c(1, 1))))
  }
}
