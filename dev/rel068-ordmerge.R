# The interplay of lanes fixes and ordmix on the merged tree: the
# per-threshold prior rows, a per-threshold prior on a mixture
# component, the incl_thres path and the confint() labels of an ordinal
# mixture, beside brms 2.23.0 where brms has an answer.
#
#   Rscript dev/rel068-ordmerge.R [lib]   (default rellib-r6)
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r6"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
set.seed(20261005)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n))
cls <- rbinom(n, 1, plogis(-0.4 + 0.5 * d$z))
lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
show <- function(tag, expr) {
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  [warning]", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("==", tag, "\n")
  print(r)
  invisible(r)
}
rows <- function(p) {
  p <- as.data.frame(p)
  keep <- p$class %in% c("Intercept", "delta")
  sort(paste(p$class, p$coef, p$group, p$dpar, sep = "|")[keep])
}
for (ord in c("none", "mu")) {
  fam <- mixture(cumulative(), sratio(), order = ord)
  show(paste("default_prior rows, order =", ord),
       rows(default_prior(bf(y ~ x), data = d, family = fam)))
  if (requireNamespace("brms", quietly = TRUE)) {
    bp <- suppressMessages(brms::default_prior(
      brms::bf(y ~ x), data = d,
      family = brms::mixture(brms::cumulative(), brms::sratio(),
                             order = ord)))
    show(paste("brms rows, order =", ord), rows(bp))
  }
}
f1 <- show("fit order none", frm(bf(y ~ x),
                                 family = mixture(cumulative(), sratio()),
                                 data = d))
show("confint rownames, order none", rownames(confint(f1)))
show("vcov(full) rownames, order none", rownames(vcov(f1, full = TRUE)))
show("ord_thres_linpred on a mixture",
     frmtmb:::ord_thres_linpred(f1))
f2 <- show("fit order mu", frm(bf(y ~ x),
                               family = mixture(cumulative(), sratio(),
                                                order = "mu"), data = d))
show("confint rownames, order mu", rownames(confint(f2)))
show("fixef rownames, order mu", rownames(fixef(f2)))
pr <- set_prior("normal(-1, 0.5)", class = "Intercept", coef = "2",
                dpar = "mu2")
f3 <- show("fit with a per-threshold prior on mu2",
           frm(bf(y ~ x), family = mixture(cumulative(), sratio()),
               data = d, prior = pr))
show("prior_summary", prior_summary(f3))
show("compat cells", sapply(c("disc", "equidistant", "sum_to_zero"),
                            function(s) frm_compat(s, "mixture")$status))
