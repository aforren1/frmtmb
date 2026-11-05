# Review of lane gpgap, item 1: the new branch on specs it was not written for.
# eq_gp() on an approximate gp (no nugget line to rewrite) must not pass
# through the edited route, and an edit that matched nothing must be loud.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-misapply.R misapply.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
requireNamespace("brms"); requireNamespace("rstan")
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
smp <- read.delim("dev/equiv-data/sample.tsv")
show <- function(tag, rec) {
  cat(sprintf("MIS %-34s verdict %s, level %s, l1_ok %s, l1_error %s, by decision: %s\n",
              tag, rec$verdict, paste(rec$level, collapse = "+"), rec$l1_ok %||% NA,
              substr(rec$l1_error %||% "-", 1, 70), rec$l1_by_decision %||% "-"))
}
key <- "gaussian x hsgp"
row <- smp[smp$id == "p082", ]
orig <- E$eq_specs[[key]]
show("p082 as it stands", E$eq_measure(key, row))
E$eq_specs[[key]] <- function() {
  s <- orig()
  s$l4 <- FALSE
  s$l1_code <- E$eq_gp_frm_nugget
  s
}
show("p082 with the rewrite forced", E$eq_measure(key, row))
# an edit that changes a line and so reaches brms's model everywhere: a wrong
# rewrite on the exact pair must fail, not pass
key2 <- "lognormal x gp"
row2 <- smp[smp$id == "p094", ]
orig2 <- E$eq_specs[[key2]]
E$eq_specs[[key2]] <- function() {
  s <- orig2()
  s$l1_code <- function(code) {
    gsub("cov[n, n] += 1e-12;", "cov[n, n] += 1e-5 * square(sdgp);", code, fixed = TRUE)
  }
  s
}
show("p094 with a wrong rewrite (1e-5)", E$eq_measure(key2, row2))
cat("specs carrying l1_decision:",
    sum(vapply(names(E$eq_specs), function(k) {
      s <- tryCatch(E$eq_specs[[k]](), error = function(e) list())
      !is.null(s$l1_decision)
    }, NA)), "\n")
cat("specs carrying l1_code:", paste(names(Filter(function(k) TRUE, list())), collapse = ","),
    paste(names(E$eq_specs)[vapply(names(E$eq_specs), function(k) {
      s <- tryCatch(E$eq_specs[[k]](), error = function(e) list())
      is.function(s$l1_code)
    }, NA)], collapse = "; "), "\n")
