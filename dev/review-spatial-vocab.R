# Reviewer: brms 2.23.0's vocabulary, recomputed from the installed brms
# and its source, against the lane's fixed data.
.libPaths(c("C:/Users/adf44/source/r/wt-spatial-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("brms", as.character(packageVersion("brms")), "\n")
ns <- asNamespace("brms")
voc <- frmtmb:::brms_vocabulary
fam <- sub("^\\.family_", "", grep("^\\.family_", ls(ns, all.names = TRUE),
                                   value = TRUE))
fam <- setdiff(fam, "info")
cat("families in brms:", length(fam), " lane:", length(voc$families), "\n")
cat("  brms only:", setdiff(fam, voc$families), "\n")
cat("  lane only:", setdiff(voc$families, fam), "\n")
# which family names brmsfamily() itself accepts (aliases too)
alias_try <- c("Gamma", "inverse.gaussian", "zero_inflated_negbinomial",
               "hurdle_cumulative", "student_t", "lognormal", "von_mises",
               "categorical", "skew_normal", "exponential", "geometric",
               "logistic", "gen_extreme_value", "multinomial", "beta",
               "ordbeta", "zero_one_inflated_beta", "xbeta")
for (a in alias_try) {
  r <- tryCatch({brms::brmsfamily(a); "ok"}, error = function(e) "refused")
  cat("  brmsfamily(", a, "):", r, "\n")
}
exps <- getNamespaceExports(ns)
resp <- sub("^resp_", "", grep("^resp_", exps, value = TRUE))
cat("resp_* exports:", length(resp), " lane aterms:", length(voc$aterms), "\n")
cat("  brms only:", setdiff(resp, voc$aterms), "\n")
cat("  lane only:", setdiff(voc$aterms, resp), "\n")
# addition terms brms's parser knows, from its source
src <- list.files("C:/Users/adf44/source/r/frmtmb-wt-spatial/dev/brms-suite",
                  pattern = "formula-ad.R$", recursive = TRUE, full.names = TRUE)
cat("formula-ad.R:", src, "\n")
if (length(src)) {
  txt <- readLines(src[1])
  i <- grep("all_aterms|ad_terms|\\.ad_funs|aterms <-|resp_", txt)
  cat(head(txt[i], 25), sep = "\n")
}
cands <- c("s", "t2", "te", "ti", "gp", "mo", "me", "mi", "cs", "cse",
           "car", "sar", "fcor", "ar", "ma", "arma", "cosy", "unstr", "gr",
           "mm", "mmc", "sbhaz", "lf", "nlf", "acformula", "offset", "sp",
           "cor_ar", "cor_arma", "cor_car", "cor_sar", "cor_fixed",
           "cor_cosy", "cor_bsts", "cor_arr", "cor_icar", "cor_ma")
cat("candidate specials exported by brms:\n")
for (s in cands) {
  cat(sprintf("  %-10s exported=%s lane=%s\n", s, s %in% exps,
              s %in% voc$specials))
}
cat("specials in lane, not exported by brms:",
    setdiff(voc$specials, exps), "\n")
