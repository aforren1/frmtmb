# Compare the draws of dev/optima-p1-nomo.R, base against lane.
a <- readRDS("dev/optima-log/p1-nomo-base.rds")
b <- readRDS("dev/optima-log/p1-nomo-lane.rds")
for (nm in names(a)) {
  same_names <- identical(colnames(a[[nm]]), colnames(b[[nm]]))
  cat(sprintf("%-18s columns identical %s; draws identical %s; max |diff| %.3g\n",
              nm, same_names, identical(unname(a[[nm]]), unname(b[[nm]])),
              max(abs(a[[nm]] - b[[nm]]))))
}
