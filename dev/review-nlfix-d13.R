# Fix round 2: bobyqa under frm_allfit(start_from_mle = FALSE), both builds; the saddle table's
# reading from my own run of dev/review-nlfix-saddle.R.
source("dev/review-nlfix-common.R")
set.seed(9202); n <- 300; d <- data.frame(x = round(rnorm(n), 3)); d$y <- 1.5 * exp(0.6 * d$x) + rnorm(n, 0, 0.5)
f <- suppressWarnings(frm(bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE), data = d, start = list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 0.2))))
x <- tt(frm_allfit(f, start_from_mle = FALSE)); print(if (iserr(x)) x$msg else x)
g <- suppressWarnings(frm(bf(y ~ a + b * x, a + b ~ 1, nl = TRUE), data = d))
x <- tt(frm_allfit(g, start_from_mle = FALSE)); print(if (iserr(x)) x$msg else x)
if (arm == "after" && dir.exists("dev/review-nlfix-saddle-out")) {
  fs <- list.files("dev/review-nlfix-saddle-out", "txt$", full.names = TRUE)
  r <- do.call(rbind, lapply(fs, function(p) read.table(p, header = TRUE, sep = "\t", stringsAsFactors = FALSE)))
  print(names(r)); cat("rows", nrow(r), "\n")
  r$gap <- r$brms_best - r$default
  cat("at or above brms's best (within 1e-3):", sum(r$gap <= 1e-3, na.rm = TRUE), "of", sum(is.finite(r$gap)), "\n")
  print(r[which(r$gap > 1e-3), c("family", "body", "default", "brms_best", "gap")], row.names = FALSE)
}
