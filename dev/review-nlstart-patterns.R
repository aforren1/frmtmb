# Review of lane nlstart: the 4-pattern cap. For short rows of the fresh
# designs (review-nlstart-fresh.R), the default fit with the shipped
# patterns, with all 2^k patterns, and with each single pattern (the
# probe's nl_sign_patterns() replaced in this process only).
# Usage: REV_ARM=after Rscript dev/review-nlstart-patterns.R
source("dev/review-nlstart-common.R")
src <- readLines("dev/review-nlstart-fresh.R")
# the generator and the bodies of the fresh designs, without the loop
stop_at <- grep("^cat\\(\"ROWHEAD", src) - 1L
start_at <- grep("^F <- list", src)
eval(parse(text = src[start_at:stop_at]))
ns <- asNamespace("frmtmb")
orig <- get("nl_sign_patterns", ns)
set_pat <- function(fun) {
  unlockBinding("nl_sign_patterns", ns)
  assign("nl_sign_patterns", fun, envir = ns)
  lockBinding("nl_sign_patterns", ns)
}
all_pat <- function(k) {
  lapply(seq_len(2^k), function(i) {
    ifelse(bitwAnd(i - 1L, 2^(seq_len(k) - 1L)) > 0, -1, 1)
  })
}
rows <- list(c("bernoulli", "two_sat_rev", 33), c("bernoulli", "two_sat_rev", 22),
             c("poisson", "two_sat_rev", 11), c("binomial", "two_sat_rev", 33),
             c("gamma", "two_sat", 11), c("gamma", "two_sat", 44),
             c("cumulative", "bigx", 44), c("bernoulli", "mm", 55))
for (r in rows) {
  fam <- r[1]; body <- r[2]; seed <- as.numeric(r[3])
  d <- gen(fam, body, seed * 1000 + match(body, names(B)))
  fm <- form(fam, body)
  set_pat(orig)
  f0 <- tt(frm(fm, family = F[[fam]], data = d))
  set_pat(all_pat)
  fa <- tt(frm(fm, family = F[[fam]], data = d))
  k <- length(B[[body]]$p)
  cat(sprintf("== %s %s %d: shipped %.4f (multi %s) | all %d patterns %.4f (multi %s)\n",
              fam, body, seed, llv(f0), any(attr(f0, "mc") == "frmtmb_nl_multistart"),
              2^k, llv(fa), any(attr(fa, "mc") == "frmtmb_nl_multistart")))
  shipped <- vapply(orig(k), paste, "", collapse = "")
  for (p in all_pat(k)) {
    set_pat(function(kk) list(p))
    fp <- tt(frm(fm, family = F[[fam]], data = d))
    lab <- paste(ifelse(p > 0, "+", "-"), collapse = "")
    cat(sprintf("   pattern %s %s: %.4f\n", lab,
                if (paste(p, collapse = "") %in% shipped) "shipped" else "  -    ",
                llv(fp)))
  }
  flush.console()
}
set_pat(orig)
