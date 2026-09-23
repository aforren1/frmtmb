# Punch round 2 blocker check: every job that is NOT a declaring dpar
# must be bitwise identical to base, and the poisson scale sweep is
# checked against glm() on both builds.
d <- "C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/"
a <- readRDS(paste0(d, "noint-base.rds"))
b <- readRDS(paste0(d, "noint-fixed.rds"))
stopifnot(identical(names(a), names(b)))
fam_jobs <- setdiff(names(a), grep("^glm_", names(a), value = TRUE))
# The blocker is about families that did NOT ask for the fix. Every
# skew_normal job may move, because that family's own initializers
# changed (the residual start), which is the lane's purpose; what may
# NOT move is any other family. dev/skewinit-noint-why.R shows the
# starts: sigma ~ 0 + g and mu ~ 0 + g are all-zero on BOTH builds, so
# a non-declaring dpar is still skipped, and only alpha ~ 0 + g is
# placed.
declaring <- grep("^sn_", names(a), value = TRUE)

cat(sprintf("%-18s %-9s %14s %s\n", "job", "bitwise", "d logLik",
            "conv/warn base -> lane"))
same <- 0L
moved <- character(0)
for (nm in fam_jobs) {
  x <- a[[nm]]; y <- b[[nm]]
  if (!is.null(x$error) || !is.null(y$error)) {
    cat(sprintf("%-18s ERROR\n", nm)); next
  }
  bit <- identical(x$ll, y$ll) && identical(x$par, y$par)
  if (bit) same <- same + 1L else moved <- c(moved, nm)
  cat(sprintf("%-18s %-9s %14s %d/%d -> %d/%d\n", nm,
              if (bit) "identical" else "MOVED",
              format(y$ll - x$ll, digits = 6), x$conv, x$warns,
              y$conv, y$warns))
}
cat("\nbitwise identical:", same, "of", length(fam_jobs), "\n")
cat("moved:", if (length(moved)) paste(moved, collapse = " ") else "none",
    "\n")
unexpected <- setdiff(moved, declaring)
cat("NON-skew_normal jobs that moved:", length(unexpected),
    if (length(unexpected)) paste(unexpected, collapse = " "), "\n")
if (!length(unexpected)) {
  cat("PASS: the blast radius is skew_normal only\n")
} else {
  cat("FAIL: another family still moves\n")
}
for (nm in intersect(declaring, moved)) {
  cat(sprintf("  %s: base %.9f  lane %.9f  d %+.6g\n", nm,
              a[[nm]]$ll, b[[nm]]$ll, b[[nm]]$ll - a[[nm]]$ll))
}

cat("\n-- the poisson scale sweep against glm()\n")
cat(sprintf("%-8s %18s %18s %18s %12s %12s\n", "scale", "glm", "base",
            "lane", "base-glm", "lane-glm"))
for (s in c("s0", "s2", "s4", "s6")) {
  g <- a[[paste0("glm_", s)]]$ll
  cat(sprintf("%-8s %18.9f %18.9f %18.9f %12.6f %12.6f\n", s, g,
              a[[paste0("pois_", s)]]$ll, b[[paste0("pois_", s)]]$ll,
              a[[paste0("pois_", s)]]$ll - g,
              b[[paste0("pois_", s)]]$ll - g))
}
cat("\nlane equals base on every poisson scale:",
    all(vapply(c("s0", "s2", "s4", "s6"), function(s) {
      identical(a[[paste0("pois_", s)]]$ll, b[[paste0("pois_", s)]]$ll)
    }, TRUE)), "\n")
