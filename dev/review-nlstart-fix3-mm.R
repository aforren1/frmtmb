# Review of lane nlstart, third short check: a first fit that is short and
# that the probe does not improve, at covariate scale 1 and 1e-6. The
# Michaelis-Menten rows of review-nlstart-fresh.R (bernoulli, binomial),
# with xp and K rescaled by s. Records warnings by class and the messages.
# Usage: REV_ARM=after|before Rscript dev/review-nlstart-fix3-mm.R [ref]
do_ref <- isTRUE(commandArgs(trailingOnly = TRUE)[1] == "ref")
source("dev/review-nlstart-common.R")
src <- readLines("dev/review-nlstart-fresh.R")
eval(parse(text = src[grep("^F <- list", src):(grep("^cat\\(\"ROWHEAD", src) - 1L)]))
cat("ROWHEAD fam seed s ll conv nwarn wclass multi ref\n")
for (fam in c("bernoulli", "binomial")) for (seed in c(11, 22, 33, 44, 55)) {
  d0 <- gen(fam, "mm", seed * 1000 + match("mm", names(B)))
  lhs <- if (fam == "binomial") "y | trials(nt)" else "y"
  fm <- bf(as.formula(paste(lhs, "~ Vm * xs / (K + xs)")), Vm + K ~ 1, nl = TRUE)
  for (s in c(1, 1e-6)) {
    d <- d0
    d$xs <- d$xp * s
    r <- tt(frm(fm, family = F[[fam]], data = d))
    best <- NA
    if (do_ref) {
      set.seed(seed * 13 + 7)
      sts <- c(lapply(1:16, function(i) rnorm(2) * c(2, 2 * s) * if (i > 8) 3 else 1),
               list(c(3, s)))
      best <- max(vapply(sts, function(st) llv(tt(frm(fm, family = F[[fam]],
        data = d, se = FALSE, start = list(beta = c("Vm_(Intercept)" = st[1],
                                                    "K_(Intercept)" = st[2]))))),
        0), na.rm = TRUE)
    }
    cat(sprintf("ROW %s %d %g %.6f %s %d %s %d %.6f\n", fam, seed, s, llv(r),
                if (iserr(r)) NA else r$opt$convergence, length(attr(r, "w")),
                if (length(attr(r, "wc"))) paste(unique(attr(r, "wc")), collapse = "+") else "-",
                as.integer(any(attr(r, "mc") == "frmtmb_nl_multistart")), best))
    flush.console()
  }
}
cat("DONE\n")
