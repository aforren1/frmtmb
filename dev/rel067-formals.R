.libPaths(c("C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
ns <- asNamespace("frmtmb"); nss <- asNamespace("frmtmb.sample")
bns <- asNamespace("brms")
cmp <- function(lab, a, b) {
  fa <- setdiff(names(formals(a)), "..."); fb <- setdiff(names(formals(b)), "...")
  cat(sprintf("%s\n  shared: %s\n  brms-only: %s\n  frmtmb-only: %s\n", lab,
              paste(intersect(fa, fb), collapse = ", "),
              paste(setdiff(fb, fa), collapse = ", "),
              paste(setdiff(fa, fb), collapse = ", ")))
}
for (f in c("acat", "cratio", "cumulative", "sratio"))
  cmp(f, get(f, ns), get(f, bns))
cmp("nsamples draws", get("nsamples.frmtmb_draws", nss), get("nsamples.brmsfit", bns))
cmp("posterior_samples draws", get("posterior_samples.frmtmb_draws", nss),
    get("posterior_samples.brmsfit", bns))
cmp("parnames fit", get("parnames.frmtmb_fit", ns), get("parnames.brmsfit", bns))
cmp("conditional_effects fit", get("conditional_effects.frmtmb_fit", ns),
    get("conditional_effects.brmsfit", bns))
cmp("fitted fit", get("fitted.frmtmb_fit", ns), get("fitted.brmsfit", bns))
extra <- commandArgs(trailingOnly = TRUE)
for (f in extra) {
  a <- if (exists(f, ns)) get(f, ns) else get(f, nss)
  b <- get(f, bns)
  cmp(f, a, b)
}
