# Lane wt-mvprior: compare dev/mvprior-sample-draws.R's two arms.
root <- "C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-log"
a <- readRDS(file.path(root, "sample-draws-lane.rds")); b <- readRDS(file.path(root, "sample-draws-base.rds"))
for (nm in names(a)) {
  cat(sprintf("%-13s cols identical %s | draws identical %s | max abs diff %s | prior_summary identical %s (rows lane %d base %d)\n",
    nm, identical(a[[nm]]$cols, b[[nm]]$cols), identical(a[[nm]]$draws, b[[nm]]$draws),
    format(max(abs(a[[nm]]$draws - b[[nm]]$draws))), identical(a[[nm]]$ps, b[[nm]]$ps),
    NROW(a[[nm]]$ps), NROW(b[[nm]]$ps)))
}
