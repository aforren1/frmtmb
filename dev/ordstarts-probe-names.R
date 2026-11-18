# Lane ordstarts fix round: what the outer parameters of an ordinal
# mixture with random effects are called, and how each component's
# thresholds map from their raw values, so a start can permute the
# components. One seed per design; prints only.
#   Rscript dev/ordstarts-probe-names.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
q <- function(e) suppressWarnings(suppressMessages(e))
d <- mk(3, n = 300, mix = TRUE)
d$g <- factor(rep_len(1:20, 300))
for (f in list(bf(y ~ x + z + (1 | g)), bf(y ~ x + z + (1 + x | g)))) {
  for (fam in list(mixture(cumulative(), sratio()),
                   mixture(cumulative(), sratio(), acat()))) {
    A <- q(frm(f, family = fam, data = d,
               control = frmtmb_control(starts = 0)))
    nm <- ns$outer_par_map(A)$names
    cat("\n", deparse(f$formula), A$spec$responses[[1]]$family$family, "\n")
    print(data.frame(nm = nm, raw = names(A$opt$par),
                     par = signif(A$opt$par, 3)))
    cat("estimates:", names(A$estimates), "\n")
    for (b in A$frame$re_blocks) {
      cat("re block:", paste(names(b), collapse = ","), "\n")
      cat("  dpar", format(b$dpar %||% b$lp %||% NA), "terms",
          paste(b$cnms %||% b$terms %||% "", collapse = ","), "\n")
    }
    resp <- names(A$spec$responses)[1]
    K <- A$spec$responses[[resp]]$family$mix$ord$K
    for (k in seq_len(K)) {
      lp <- A$frame$linpreds[[ns$linpred_key(resp, paste0("mu", k))]]
      ob <- ns$ord_lp_block(A$frame, A$spec, lp)
      raw <- as.numeric(A$estimates[[ob$comp]])
      v <- as.numeric(ns$ord_threshold_values(ob$fam, raw))
      inv <- c(v[1], log(diff(v)))
      cat("  comp", k, ob$comp, "ordered", ob$ordered, "raw", signif(raw, 4),
          "values", signif(v, 4), "| identity", isTRUE(all.equal(raw, v)),
          "| log-diff inverse", isTRUE(all.equal(raw, inv)), "\n")
    }
  }
}
