# Lane reslot, item 3: what prior = "flat" resolves to, one model per
# class. Usage: Rscript dev/reslot-flat-probe.R before|after. Seed 1.
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
dd <- re_data(1L)
dd$ord <- cut(dd$y, 4, labels = FALSE)
dd$yp <- abs(dd$y) + 0.1
show <- function(label, f, ...) {
  fit <- suppressWarnings(suppressMessages(frm(f, data = dd, ...)))
  r <- suppressWarnings(suppressMessages(
    frmtmb.sample:::sample_resolve_priors(fit, "flat")))
  cat("--", label, "\n")
  if (is.null(r$ri)) cat("  no entries\n") else show_entries(r$ri)
  if (!is.null(r$effective)) {
    print(as.data.frame(r$effective)[, c("prior", "class", "dpar", "lb",
                                         "ub")], row.names = FALSE)
  }
}
show("(1 + x + z | g)", bf(y ~ x + (1 + x + z | g)))
show("ar(cov = TRUE)", bf(y ~ x + ar(t, g, cov = TRUE)))
show("rescor", bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE))
show("Gamma shape", bf(yp ~ x), family = Gamma("log"))
show("student nu", bf(y ~ x), family = student())
show("cumulative equidistant", bf(ord ~ x),
     family = cumulative(threshold = "equidistant"))
show("gp", bf(y ~ gp(x)))
