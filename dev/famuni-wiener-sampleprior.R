# What frmtmb.sample's defaults and its disclosure say about a wiener
# model's bs, ndt and bias, against brms 2.23.0's default_prior() rows
# (dev/famuni-out/wiener-brms.txt). Seed 20261013.
# Output: dev/famuni-out/wiener-sampleprior.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all(
  "C:/Users/adf44/source/r/frmtmb-wt-famuni/extensions/frmtmb.eam",
  quiet = TRUE))
suppressMessages(library(frmtmb.sample))
set.seed(20261013)
n <- 80
d <- data.frame(x = rnorm(n))
sim <- RWiener::rwiener(n, alpha = 1.4, tau = 0.25, beta = 0.5, delta = 0.8)
d$rt <- sim$q
d$d <- as.integer(sim$resp == "upper")
ns <- asNamespace("frmtmb.sample")
for (fo in list(bf(rt | dec(d) ~ x),
                bf(rt | dec(d) ~ x, bs ~ x, ndt ~ 1, bias ~ 1))) {
  fit <- frm(fo, family = wiener(), data = d)
  cat("\n== ", deparse(fo$formula), "+",
      paste(names(fo$pforms), collapse = ","), "\n")
  pl <- ns$default_priors_for(fit)
  cat("default_priors_for():\n")
  for (s in unclass(pl)) {
    cat("  class", s$class, "dpar", s$dpar, ":", s$dist$kind,
        paste(unlist(s$dist[-1L]), collapse = ", "), "\n")
  }
  cat("default_prior_notes():\n")
  nt <- ns$default_prior_notes(fit)
  if (!length(nt)) cat("  (none)\n") else cat(paste0("  ", nt), sep = "\n")
}
