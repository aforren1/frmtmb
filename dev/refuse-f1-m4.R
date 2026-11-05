# M4: the review's 24 beta fits (dev/review-refuse-misc.R, seed 10200): the
# unread route, then a read mi() on the same draws, with a gaussian reader.
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(10200)
n <- 200
x <- rnorm(n)
res <- list()
for (e0 in c(2, 3, 4, 5, -2, -3, -4, -5)) for (phi in c(0.3, 2, 8)) {
  m <- plogis(e0 + 0.5 * x)
  y <- pmin(pmax(rbeta(n, m * phi, (1 - m) * phi), 1e-10), 1 - 1e-10)
  d <- data.frame(x = x, y = y)
  # the reader's own noise makes y's coefficient identified at the mean's
  # end, where y hardly moves on its own scale
  d$z <- 0.5 + 2 * qlogis(y) / 10 + rnorm(n, 0, 0.5)
  d$y[seq(5, 195, 10)] <- NA
  a <- tryCatch(suppressWarnings(frm(bf(y | mi() ~ x), family = Beta(),
                                     data = d)), error = function(e) e)
  b <- tryCatch(suppressWarnings(frm(bf(y ~ x), family = Beta(),
                                     data = d[!is.na(d$y), ])),
                error = function(e) e)
  r <- tryCatch(suppressWarnings(frm(
    mvbf(bf(z ~ mi(y)), bf(y | mi() ~ x, family = Beta()), rescor = FALSE),
    data = d)), error = function(e) e)
  lab <- function(f) {
    if (inherits(f, "error")) {
      paste("ERROR", substr(conditionMessage(f), 1, 50))
    } else {
      sprintf("logLik %.3f conv %d", logLik(f), f$opt$convergence)
    }
  }
  cat(sprintf("eta0=%+d phi=%-4g unread [%s] cc [%s] read [%s]\n", e0, phi,
              lab(a), lab(b), lab(r)))
  res[[length(res) + 1L]] <- c(u = !inherits(a, "error"),
                               c = !inherits(b, "error"),
                               r = !inherits(r, "error"))
}
cat("fits without an error, of", length(res), ": ",
    paste(names(res[[1]]), rowSums(do.call(cbind, res))), "\n")
