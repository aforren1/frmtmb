# Lane hollow: what a fit's $data (the model frame) carries, and whether
# a model frame rebuilt from the raw data reproduces it.
#   Rscript dev/hollow-frame-probe.R > dev/hollow-out/frame-probe.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
d6 <- as.data.frame(brms:::brmsfit_example6$data)
attr(d6, "terms") <- NULL
fit6 <- suppressWarnings(frm(bf(volume ~ Trt + gp(Age, by = Trt, gr = TRUE)) +
  gaussian() + bf(count ~ Trt + Age) + poisson() + set_rescor(FALSE),
  data = d6))
show <- function(f, label) {
  cat("\n==", label, "==\n")
  dd <- f$data
  cat("class:", class(dd), " dim:", dim(dd), "\n")
  cat("names:", names(dd), "\n")
  cat("attributes:", setdiff(names(attributes(dd)),
                             c("names", "row.names", "class")), "\n")
  cat("identical(model.frame(fit), fit$data):",
      identical(model.frame(f), f$data), "\n")
  tt <- attr(dd, "terms")
  if (!is.null(tt)) {
    re <- tryCatch(stats::model.frame(tt, data = f$frame$raw_data %||%
                                        NULL), error = function(e) e)
  }
  invisible(NULL)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
show(fit6, "fixture 6 (mv, gp)")
cat("frame elements:", names(fit6$frame), "\n")

set.seed(1)
dd <- data.frame(x = rnorm(50), Age = runif(50, 1, 3),
                 g = factor(rep(1:5, 10)), w = runif(50, 0.5, 2))
dd$y <- rpois(50, exp(0.3 + 0.4 * dd$x + log(dd$Age)))
dd$extra <- letters[rep(1:5, 10)]
fo <- frm(bf(y | weights(w) ~ poly(x, 2) + offset(log(Age)) + (1 | g)),
          family = poisson(), data = dd)
show(fo, "offset, poly, weights")
str(fo$data, give.attr = FALSE)
cat("refit data columns:",
    names(frmtmb:::fit_refit_data(fo, "data")), "\n")
cat("PROBE END\n")
