a <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-log/mo-lane.rds")
b <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-log/mo-base.rds")
for (nm in names(a)) cat(sprintf("%-22s identical par %s vcov %s fitted %s warn %s\n", nm,
  identical(a[[nm]]$par, b[[nm]]$par), identical(a[[nm]]$vcov, b[[nm]]$vcov),
  identical(a[[nm]]$fitted, b[[nm]]$fitted), identical(a[[nm]]$warn, b[[nm]]$warn)))
