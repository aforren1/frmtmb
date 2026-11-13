# Lane gpstop: summary of dev/gpstop-ng.R on the 150 rows. Logs: the lane
# with the 15 reference starts (ng-lane-a.txt, ng-lane-b.txt) and the base
# default fit (ng-base.txt). Reference "all": the best of both defaults and
# every start; "clean": the same over the defaults and the starts that
# ended at code 0 with no warning. Short: more than 1e-3 below; silent:
# code 0 and no warning. Moves: lane against base by more than 1e-3.
o <- "C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/gpstop-out"
rd <- function(f) {
  l <- grep("^ANG ", readLines(file.path(o, f)), value = TRUE)
  l <- l[!grepl("ERROR", l)]
  v <- strsplit(l, " +")
  num <- function(key) vapply(v, function(x) {
    i <- match(key, x)
    if (is.na(i)) NA_real_ else as.numeric(x[i + 1L])
  }, 0)
  data.frame(key = vapply(v, function(x) paste(x[3:6], collapse = " "), ""),
             case = vapply(v, `[`, "", 3L), curve = vapply(v, `[`, "", 4L),
             ll = num("ll"), best = num("best"), bestc = num("bestc"),
             conv = num("conv"), w = num("w"), wall = num("wall"),
             stringsAsFactors = FALSE)
}
l <- rbind(rd("ng-lane-a.txt"), rd("ng-lane-b.txt"))
b <- rd("ng-base.txt")
b <- b[match(l$key, b$key), ]
stopifnot(identical(b$key, l$key))
ref_all <- pmax(l$ll, b$ll, l$best, na.rm = TRUE)
ref_cln <- pmax(l$ll, b$ll, l$bestc, na.rm = TRUE)
tol <- 1e-3
silent <- function(x) x$conv == 0 & x$w == 0
for (grp in list(c("all", "."), c("student", "^student"),
                 c("step curve", " step "))) {
  i <- grepl(grp[2], paste0(l$key, " "))
  for (rf in c("all", "clean")) {
    ref <- if (rf == "all") ref_all else ref_cln
    sh <- function(x) i & x$ll < ref - tol
    cat(sprintf(paste0("NGS %-10s ref %-5s rows %3d | short base %2d (silent %2d)",
                       " worst %8.3f | lane %2d (silent %2d) worst %8.3f\n"),
                grp[1], rf, sum(i), sum(sh(b)), sum(sh(b) & silent(b)),
                min(c(0, (b$ll - ref)[i])), sum(sh(l)), sum(sh(l) & silent(l)),
                min(c(0, (l$ll - ref)[i]))))
  }
}
dw <- l$ll < b$ll - tol
up <- l$ll > b$ll + tol
cat(sprintf("NGS moves: lane below base %d, above %d; equal at 1e-6 %d of %d\n",
            sum(dw), sum(up), sum(abs(l$ll - b$ll) <= 1e-6), nrow(l)))
for (k in which(dw | up)) {
  cat(sprintf("NGSM %s base %.6f lane %.6f\n", l$key[k], b$ll[k], l$ll[k]))
}
for (k in which(l$ll < ref_all - tol)) {
  cat(sprintf(paste0("NGSS %s lane %.4f conv %d w %d | best %.4f (%.3f) | ",
                     "clean %.4f (%.3f)\n"),
              l$key[k], l$ll[k], l$conv[k], l$w[k], ref_all[k],
              l$ll[k] - ref_all[k], ref_cln[k], l$ll[k] - ref_cln[k]))
}
cat(sprintf("NGS wall of the default fit, separate processes: base %.1f s, lane %.1f s\n",
            sum(b$wall), sum(l$wall)))
