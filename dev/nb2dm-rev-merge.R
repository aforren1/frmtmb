# Reviewer: a trial merge of this lane with lane cdfs (snapshot of its
# worktree, 2026-10-20 ~12:00), built in the session scratchpad from
# 4854e425 plus both diffs, loaded with pkgload::load_all() (nothing
# installed). What does negbinomial2 inherit, and what fails?
#   Rscript dev/nb2dm-rev-merge.R [test file ...]
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
tree <- paste0("C:/Users/adf44/AppData/Local/Temp/3/claude/",
               "c--Users-adf44-source-r-frmtmb/",
               "7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/",
               "nb2dm-rev-merge/tree")
suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = FALSE))
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R",
       echo = FALSE, local = TRUE)
cat("frmtmb from", getNamespaceInfo("frmtmb", "path"), "\n")
a <- commandArgs(trailingOnly = TRUE)
if (length(a)) {
  suppressMessages(library(testthat))
  for (f in a) {
    r <- as.data.frame(test_file(file.path(tree, "tests/testthat", f),
                                 package = "frmtmb",
                                 env = testthat::test_env("frmtmb"),
                                 reporter = "silent",
                                 stop_on_failure = FALSE))
    cat("RESULT merge ", f, " pass=", sum(r$passed), " fail=",
        sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
        "\n", sep = "")
    bad <- r[r$failed > 0 | r$error, "test"]
    if (length(bad)) cat("FAILING BLOCK:", bad, sep = "\n  ")
    cat("\n")
  }
  quit(save = "no")
}
nb <- negbinomial()
nb2 <- brmsfamily("negbinomial2")
probe("slots", list(nb_lcdf = !is.null(nb$lcdf), nb2_lcdf = !is.null(nb2$lcdf),
                    nb2_lccdf = !is.null(nb2$lccdf),
                    nb2_trunc_mean = !is.null(nb2$post$trunc_mean_fn),
                    nb2_aterms = nb2$accepts_aterms))
probe("lcdf at shape = 1 / sigma vs pnbinom",
      c(nb2 = nb2$lcdf(4, list(mu = 2, sigma = 0.25), list()),
        ref = pnbinom(4, size = 4, mu = 2, log.p = TRUE)))
set.seed(5)
n <- 300
dd <- data.frame(x = rnorm(n))
dd$y <- rnbinom(n, mu = exp(1 + 0.4 * dd$x), size = 2)
dd$cc <- ifelse(dd$y > 6, 1L, 0L)
dd$yc <- pmin(dd$y, 6)
probe("cens(): negbinomial vs negbinomial2 logLik", {
  f1 <- frm(yc | cens(cc) ~ x, family = negbinomial(), data = dd)
  f2 <- frm(yc | cens(cc) ~ x, family = nb2, data = dd)
  c(nb = logLik(f1), nb2 = logLik(f2),
    shape_vs_inv_sigma = c(exp(fixef_by_dpar(f1)$shape), 1 / sigma(f2)))
})
probe("trunc(): negbinomial vs negbinomial2 logLik and fitted", {
  d2 <- dd[dd$y >= 1, ]
  f1 <- frm(y | trunc(lb = 1) ~ x, family = negbinomial(), data = d2)
  f2 <- frm(y | trunc(lb = 1) ~ x, family = nb2, data = d2)
  c(nb = logLik(f1), nb2 = logLik(f2),
    fitted_ratio = max(abs(fitted(f2)[, 1] / fitted(f1)[, 1] - 1)))
})
probe("frm_compat cens()/trunc() rows", {
  rows <- function(fam) {
    t <- frm_compat(fam)
    t[t$feature_b %in% c("cens()", "trunc()"), c("feature_a", "feature_b",
                                                  "status")]
  }
  rbind(rows("negbinomial"), rows("negbinomial2"))
})
probe("residuals quantile/osa nb2 on the merge", {
  f2 <- frm(y ~ x, family = nb2, data = dd)
  head(dharma_residuals(f2, nsim = 20)$scaledResiduals)
})
