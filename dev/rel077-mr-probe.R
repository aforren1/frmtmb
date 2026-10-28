# mointerop x reslot: brms's (flat) on the simplex of a two-value mo()
# (a simplex of one weight). Rscript dev/rel077-mr-probe.R [lib]
a <- commandArgs(TRUE)
lib <- if (length(a)) a[1] else "C:/Users/adf44/source/r/rellib-r15"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), "\n")
set.seed(3)
dm <- data.frame(m2 = rep(1:2, 30), k = sample(1:4, 60, TRUE))
dm$y <- 0.5 * dm$m2 + rnorm(60)
fl <- list(bf(y ~ mo(k)))
if (!grepl("reslot|r14", lib)) fl <- c(list(bf(y ~ mo(m2))), fl)
for (f in fl) {
  fit <- frm(f, data = dm, dry_run = "objective")
  co <- if (grepl("m2", deparse(f$formula))) "mom21" else "mok1"
  defs <- frmtmb.sample:::default_priors_for(fit)
  for (s in unclass(defs)) if (s$class == "simo") str(s$dist)
  for (pr in list(set_prior("", class = "simo", coef = co),
                  set_prior("dirichlet(2)", class = "simo", coef = co))) {
    r <- tryCatch(suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, pr)),
                  error = function(e) conditionMessage(e))
    cat(deparse(f$formula), "coef", pr[[1]]$coef, ":",
        if (is.character(r)) r else "resolved", "\n")
  }
  r <- tryCatch(frm(f, data = dm, prior = set_prior("", class = "simo", coef = co)),
                error = function(e) conditionMessage(e))
  cat("  frm():", if (is.character(r)) r else "fits", "\n")
}
if (requireNamespace("brms", quietly = TRUE)) {
  sc <- brms::stancode(brms::bf(y ~ mo(m2)), data = dm,
    prior = brms::set_prior("", class = "simo", coef = "mom21"))
  cat("brms:", grep("simo", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n  ")
}
if (requireNamespace("brms", quietly = TRUE)) {
  sd <- brms::standata(brms::bf(y ~ mo(k)), data = dm,
    prior = brms::set_prior("", class = "simo", coef = "mok1"))
  cat("brms con_simo_1 under \"\":", sd$con_simo_1, "\n")
}
