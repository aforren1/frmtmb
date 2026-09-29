# Dry run of the dependency parser in the "Install the seven extensions"
# step of .github/workflows/pkgdown.yaml. Same code, minus the installs,
# so the names it would pass to install.packages() can be read.
ext <- c("frmtmb.eam", "frmtmb.latent", "frmtmb.ode",
         "frmtmb.spline", "frmtmb.coupling", "frmtmb.sample",
         "frmtmb.learn")
for (nm in ext) {
  dir <- file.path("extensions", nm)
  d <- read.dcf(file.path(dir, "DESCRIPTION"))
  f <- intersect(c("Imports", "Suggests"), colnames(d))
  p <- unlist(strsplit(paste(d[, f], collapse = ","), ","))
  p <- trimws(sub("[(].*", "", p))
  p <- p[nzchar(p)]
  p <- p[!startsWith(p, "frmtmb")]
  cat(sprintf("%-18s %2d  %s\n", nm, length(p), paste(p, collapse = " ")))
}
