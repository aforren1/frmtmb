# Reviewer round 2: the nugget decision replicated. The lane's
# dev/gpby-p1-nugget.R fits ONE data set per design. This fits 15 seeds
# of three designs at one nugget and records convergence, pdHess, the
# largest gradient, logLik and se.fit past the data, so the arms can be
# paired by seed (dev/gpby-rev2-nugget-summ.R).
# Usage: Rscript dev/gpby-rev2-nugget.R <nugget>
nug <- as.numeric(commandArgs(TRUE)[1])
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
utils::assignInNamespace("gp_nugget", nug, "frmtmb")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-gpby"
out <- file.path(wt, sprintf("dev/gpby-rev2-nugget-%s.tsv", format(nug)))
cat("design\tseed\tnugget\tconv\tpdHess\tmax_grad\tlogLik\tse_past\tse_in\n",
    file = out)
mk <- function(design, s) {
  set.seed(s)
  if (design == "x60") {
    d <- data.frame(x = round(runif(60, 0, 6), 1))
    d$y <- sin(d$x) + rnorm(60, 0, 0.3)
  } else if (design == "x200") {
    d <- data.frame(x = round(runif(200, 0, 6), 2))
    d$y <- sin(d$x) + rnorm(200, 0, 0.3)
  } else {
    d <- data.frame(x = round(runif(60, 0, 6), 1),
                    f = factor(rep(c("a", "b", "c"), 20)))
    d$y <- ifelse(d$f == "a", sin(d$x), ifelse(d$f == "b", cos(d$x),
                  0.2 * d$x)) + rnorm(60, 0, 0.3)
  }
  d
}
for (design in c("x60", "x200", "byf")) {
  for (s in 1:15) {
    d <- mk(design, s)
    fo <- if (design == "byf") bf(y ~ gp(x, by = f)) else bf(y ~ gp(x))
    r <- tryCatch({
      fit <- suppressWarnings(frm(fo, data = d))
      dg <- diagnose(fit, quiet = TRUE)
      nd <- data.frame(x = c(7, 2.55))
      if (design == "byf") nd$f <- factor("a", levels = c("a", "b", "c"))
      se <- suppressWarnings(as.numeric(frm_linpred(fit, newdata = nd,
                                                    se.fit = TRUE)$se.fit))
      c(as.integer(dg$convergence), as.integer(isTRUE(dg$pdHess)),
        dg$max_grad, as.numeric(logLik(fit)), se)
    }, error = function(e) rep(NA_real_, 6))
    cat(sprintf("%s\t%d\t%s\t%s\n", design, s, format(nug),
                paste(r, collapse = "\t")), file = out, append = TRUE)
  }
}
cat("DONE\n", file = out, append = TRUE)
