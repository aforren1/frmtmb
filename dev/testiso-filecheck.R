# Does one test file leave frmtmb's registries as it found them? Every
# extension is loaded first, so a registration the file's own loads would
# make is already in place and only the test's own change shows.
#
#   Rscript dev/testiso-filecheck.R <package> <test file>
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
a <- commandArgs(TRUE)
suppressMessages(library(a[1], character.only = TRUE))
for (p in c("frmtmb.coupling", "frmtmb.eam", "frmtmb.latent",
            "frmtmb.learn", "frmtmb.ode", "frmtmb.sample",
            "frmtmb.spline")) loadNamespace(p)
ns <- asNamespace("frmtmb")
reg <- function() {
  list(features = ns$frmtmb_compat_contrib$features,
       rules = ns$frmtmb_compat_contrib$rules,
       aterms = ns$frmtmb_aterm_registry$reg,
       checks = ns$frmtmb_frame_checks$fns,
       priors = ns$frmtmb_prior_defaults$providers,
       notices = sort(ls(ns$frmtmb_notice_state, all.names = TRUE)))
}
# a notice an earlier file gave, so a file that clears the state shows
if (nzchar(Sys.getenv("TESTISO_SEED_NOTICE"))) {
  assign("testiso_earlier_notice", TRUE, envir = ns$frmtmb_notice_state)
}
r0 <- reg()
r <- testthat::test_file(a[2], package = a[1],
                         env = testthat::test_env(a[1]), reporter = "silent")
d <- as.data.frame(r)
r1 <- reg()
ch <- names(r0)[!mapply(identical, r0, r1)]
cat(basename(a[2]), "pass", sum(d$passed), "fail", sum(d$failed), "err",
    sum(d$error), "| changed:", if (length(ch)) ch else "nothing", "\n")
for (k in ch) {
  cat("  ", k, ": ", length(r0[[k]]), " -> ", length(r1[[k]]), "\n",
      sep = "")
}
