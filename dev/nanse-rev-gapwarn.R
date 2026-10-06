# Reviewer: on the mo() seeds whose fit is below the profile maximum
# (lane's dev/nanse-log/mo-profile-base.tsv), does the lane build warn?
# Reads the lane's mo-lane-final.tsv and the trial merge's mo-merge.tsv.
#   Rscript dev/nanse-rev-gapwarn.R
w <- "C:/Users/adf44/source/r/frmtmb-wt-nanse/dev/"
p <- read.delim(paste0(w, "nanse-log/mo-profile-base.tsv"))
for (f in c("nanse-log/mo-lane-final.tsv", "nanse-rev-log/mo-merge.tsv")) {
  x <- read.delim(paste0(w, f))
  x <- x[x$form == "int", ]
  m <- merge(x, p[, c("seed", "gap")], by = "seed")
  cat("==", f, "\n")
  for (thr in c(0.01, 0.1, 0.5)) {
    big <- m$gap > thr
    cat(sprintf(paste0("gap > %.2f: %d seeds; SE warning %d; any fit warning",
                       " %d; silent %d; every SE finite %d\n"),
                thr, sum(big), sum(m$se_warn[big]),
                sum(m$n_warn_fit[big] > 0), sum(m$n_warn_fit[big] == 0),
                sum(m$se_finite[big] == m$n_par[big])))
  }
  print(m[m$seed %in% c(7, 12, 54, 113, 175),
          c("seed", "se_finite", "n_par", "se_warn", "n_warn_fit", "gap")],
        row.names = FALSE)
}
# the base build on the same seeds: were the silent, finite lane fits
# all-NaN on base (the repair made them confident) or finite already?
b <- read.delim(paste0(w, "nanse-log/mo-base.tsv"))
b <- b[b$form == "int", c("seed", "se_finite", "n_par")]
names(b)[2] <- "base_se_finite"
x <- read.delim(paste0(w, "nanse-log/mo-lane-final.tsv"))
x <- x[x$form == "int", ]
m <- merge(merge(x, p[, c("seed", "gap")], by = "seed"), b[, 1:2],
           by = "seed")
for (thr in c(0.01, 0.1, 0.5)) {
  sel <- m$gap > thr & m$n_warn_fit == 0 & m$se_finite == m$n_par
  cat(sprintf(paste0("gap > %.2f, silent, every SE finite on the lane: %d;",
                     " of these all-NaN on base %d, some-NaN %d, finite %d\n"),
              thr, sum(sel), sum(m$base_se_finite[sel] == 0),
              sum(m$base_se_finite[sel] > 0 &
                    m$base_se_finite[sel] < m$n_par[sel]),
              sum(m$base_se_finite[sel] == m$n_par[sel])))
}
