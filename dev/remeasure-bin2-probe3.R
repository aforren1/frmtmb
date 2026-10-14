# Fix round 1: the ordinal densities of distributions:238, frmtmb's
# family objects (through dev/brmsport-bin2-adapters.R's
# brms_adapt_dord()) against brms's own closed-form helpers
# (tests/testthat/helpers/inv_link_ordinal_ch.R), per family and link,
# with brms's unsorted random thresholds and with sorted ones.
#
#   Rscript dev/remeasure-bin2-probe3.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
`%||%` <- function(a, b) if (is.null(a)) b else a
Sys.setenv(FRMTMB_PORT_ROOT = getwd())
source("dev/brmsport-bin2-adapters.R")
source(file.path(brms_adapt_suite_dir, "helpers", "inv_link_ordinal_ch.R"))
cat("frmtmb", format(packageVersion("frmtmb")), "\n")
set.seed(238)
ndraws <- 5; ncat <- 4
thres_u <- matrix(rnorm(ndraws * (ncat - 1)), nrow = ndraws)
thres_s <- t(apply(thres_u, 1, sort))
eta <- rnorm(ndraws)
ch <- list(cumulative = function(te, l) inv_link_cumulative_ch(te, l),
           sratio = function(te, l) inv_link_sratio_ch(te, l),
           cratio = function(et, l) inv_link_cratio_ch(et, l),
           acat = function(et, l) inv_link_acat_ch(et, l))
for (th_name in c("unsorted", "sorted")) {
  th <- if (th_name == "sorted") thres_s else thres_u
  te <- sweep(th, 1, eta)
  et <- sweep(-th, 1, eta, FUN = "+")
  cat("\n##", th_name, "thresholds\n")
  for (f in names(ch)) {
    for (l in c("logit", "probit", "cauchit", "cloglog")) {
      d <- brms_adapt_dord(f, seq_len(ncat), eta, th, link = l)
      ref <- ch[[f]](if (f %in% c("cumulative", "sratio")) te else et, l)
      ref <- matrix(ref, nrow = ndraws)
      cat(sprintf("%-10s %-7s frmtmb NaN %2d of %2d; max abs diff %s; ",
                  f, l, sum(is.nan(d)), length(d),
                  format(max(abs(d - ref), na.rm = TRUE), digits = 3)),
          "reference negative:", sum(ref < 0), "\n")
    }
  }
}

cat("\n## cumulative logit, unsorted: where frmtmb is NaN and the reference",
    "is not negative\n")
te <- sweep(thres_u, 1, eta)
d <- brms_adapt_dord("cumulative", seq_len(ncat), eta, thres_u)
ref <- matrix(inv_link_cumulative_ch(te, "logit"), nrow = ndraws)
for (s in seq_len(ndraws)) {
  cat(sprintf("draw %d thresholds %s | frmtmb %s | brms %s\n", s,
              paste(sprintf("%.3f", thres_u[s, ]), collapse = " "),
              paste(sprintf("%.3f", d[s, ]), collapse = " "),
              paste(sprintf("%.3f", ref[s, ]), collapse = " ")))
}
