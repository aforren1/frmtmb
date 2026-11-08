# Release review 0.80.0: five sampled pairs that a lane's change can
# reach (residuals, fitted, predict), their l2 measure on both records.
#   Rscript dev/review-rel-survey2.R
rd <- function(f) utils::read.delim(f, colClasses = "character", na.strings = "NA")
Sys.setenv(GIT_OPTIONAL_LOCKS = "0")
old <- tempfile(fileext = ".tsv")
writeLines(system2("git", c("-c", "core.autocrlf=false", "show",
                            "3cd2a0c8:dev/equiv-data/results.tsv"), stdout = TRUE), old)
a <- rd(old); b <- rd("dev/equiv-data/results.tsv")
for (i in c("p074", "p101", "p093", "p079", "p018", "p096", "p092", "p028")) {
  x <- a[a$id == i, ]; y <- b[b$id == i, ]
  cat(sprintf("%s %s x %s | l2_what: %s\n  old l1_const %s l1_grad %s l2 %s l4_z %s warn %s\n  new l1_const %s l1_grad %s l2 %s l4_z %s warn %s\n",
              i, x$feature_a, x$feature_b, substr(y$l2_what, 1, 90),
              x$l1_const, x$l1_grad, x$l2_metric, x$l4_z, x$frm_warn,
              y$l1_const, y$l1_grad, y$l2_metric, y$l4_z, y$frm_warn))
}
