# fitted() standard errors of the two fits of test-predict-draw-edge.R's
# "a predictor on a zero inflation or hurdle at its end is narrowed"
# (zero_inflated_binomial with zi ~ x at zi's end; hurdle_poisson with
# hu ~ x where no y is 0), on one library: the warning, the count of NaN
# Est.Error, and the SE-check reasons. The merged build warns where the
# lanes' own builds do not (dev/round-20261111.md, conflict 2).
#
#   Rscript dev/rel082-fitted-held.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
# test-predict-draw-edge.R's pde_zib_data(0)
pde_zib_data <- function(zi, n = 200, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), N = 10L)
  p <- stats::plogis(-0.5 + 0.8 * d$x)
  d$y <- ifelse(stats::runif(n) < zi, 0L, stats::rbinom(n, d$N, p))
  d
}
src <- readLines("tests/testthat/test-predict-draw-edge.R")
i <- grep("^pde_zib_data <- function", src)
if (length(i)) {
  j <- i
  while (!grepl("^}", src[j])) j <- j + 1L
  eval(parse(text = src[i:j]))
}
d <- pde_zib_data(0)
f1 <- frm(bf(y | trials(N) ~ x, zi ~ x), family = zero_inflated_binomial(),
          data = d)
set.seed(12)
dh <- data.frame(x = stats::rnorm(300))
dh$y <- stats::rpois(300, exp(2.5 + 0.2 * dh$x))
f4 <- frm(bf(y ~ x, hu ~ x), family = hurdle_poisson(), data = dh)
for (nm in c("zi ~ x", "hu ~ x")) {
  f <- if (nm == "zi ~ x") f1 else f4
  w <- character()
  fe <- withCallingHandlers(fitted(f), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  lost <- f$sdr$se_lost
  cat(nm, ": Est.Error NaN", sum(is.nan(fe[, "Est.Error"])), "of", nrow(fe),
      "; median Est.Error", format(stats::median(fe[, "Est.Error"],
                                                 na.rm = TRUE), digits = 6),
      "; warnings", length(w), if (length(w)) paste0("(", substr(w[1], 1, 60), ")"),
      "; se_lost", if (length(lost)) paste(names(lost), lost, sep = "=",
                                            collapse = " ") else "none", "\n")
}
cat("DONE\n")
