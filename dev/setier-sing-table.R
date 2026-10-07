# Lane setier: dev/setier-singular.R's tables in 80 columns, one row per
# design and arm. msg = fits with the boundary message (in brackets:
# of them lme4-singular); warn = the SE warning; conv = another warning
# (the convergence one); se100 = a finite theta SE above 100.
#   Rscript dev/setier-sing-table.R <arm>=<tsv> ...
cat(sprintf("%-10s %-7s %4s %5s %5s %9s %5s %5s %6s\n", "arm", "design",
            "n", "lme4s", "lost", "msg", "warn", "conv", "se100"))
for (a in commandArgs(TRUE)) {
  nm <- sub("=.*", "", a)
  r <- read.delim(sub("^[^=]*=", "", a), stringsAsFactors = FALSE)
  for (d in unique(r$design)) {
    x <- r[r$design == d, ]
    cat(sprintf("%-10s %-7s %4d %5d %5d %4d (%2d) %5d %5d %6d\n", nm, d,
                nrow(x), sum(x$lme4_singular), sum(x$n_lost > 0),
                sum(x$boundary_msg > 0),
                sum(x$boundary_msg > 0 & x$lme4_singular),
                sum(x$se_warn > 0), sum(x$other_warn > 0),
                sum(is.finite(x$se_theta_max) & x$se_theta_max > 100)))
  }
}
