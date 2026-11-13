# Lane ordstart, short check item 1: what each stricter "proper" test
# would do on the 40-fit table (dev/ordstart-strict.R output). Floors: an
# exact-Hessian eigenvalue floor just above the six flat points' largest
# smallest eigenvalue (5.07e-8, dpar_formula seed 4; the review's
# q4pts.log and moves.log), and a category-cover floor of 1e-8.
#   Rscript dev/ordstart-strict-sum.R strict40.tsv
X <- utils::read.delim(commandArgs(trailingOnly = TRUE)[1])
six <- c(cs_sr_acat_seed4 = -3.94e-31, cs_sr_acat_seed1 = 6.88e-9,
         dpar_formula_seed4 = 5.07e-8, cum_sratio_14 = 3.92e-8,
         cum_sratio_18 = 3.7e-8, acat_cumcll_K7_seed6 = 9.82e-9)
six_cover <- c(cum_sratio_14 = 6.57e-7, cum_sratio_18 = 9.07e-9,
               acat_cumcll_K7_seed6 = 1.18e-8, dpar_formula_seed4 = 3.2e-8)
eig_floor <- 5.1e-8
cov_floor <- 1e-8
cat("six flat points: rejected by the eigenvalue floor", eig_floor, ":",
    sum(six < eig_floor), "of", length(six), "\n")
cat("six flat points with a cover measured: rejected by cover <", cov_floor,
    ":", sum(six_cover < cov_floor), "of", length(six_cover), "\n")
passed <- !X$degen0
moved <- abs(X$ll - X$ll0) > 1e-6
op <- options(width = 200)
show <- function(sel, what) {
  cat("\n", what, ": ", sum(sel), "\n", sep = "")
  if (any(sel)) {
    print(X[sel, c("family", "seed", "ll0", "eig0", "cover0", "ll", "eig",
                   "cover")], digits = 4, row.names = FALSE)
  }
}
show(passed & (X$eig0 < eig_floor) %in% TRUE,
     "passed defaults the eigenvalue floor would call flat (they would stay)")
show(passed & moved & (X$eig0 < eig_floor) %in% TRUE,
     "of those, fits whose upward move would be lost")
show(passed & X$cover0 < cov_floor,
     "passed defaults the cover floor would call flat")
show(passed & moved & X$cover0 < cov_floor,
     "of those, fits whose upward move would be lost")
show(!passed & moved & ((X$eig < eig_floor) %in% TRUE),
     "step-to-proper moves whose new point the eigenvalue floor rejects")
show(passed & moved & ((X$eig < eig_floor) %in% TRUE),
     "upward moves to a point the eigenvalue floor would not call proper")
