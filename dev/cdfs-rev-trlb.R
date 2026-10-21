# Reviewer: the trunc(lb = ) normalizer when F(lb) rounds to one
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
set.seed(3)
d <- data.frame(y = 10 + abs(rnorm(50)))
for (fm in list(gaussian(), student(), Gamma(link = "log"))) {
  ob <- frm(bf(y | trunc(lb = 10) ~ 1), family = fm, data = d, dry_run = "objective")$obj
  p <- ob$par
  p[1] <- if (identical(fm$family, "Gamma")) log(1) else 0   # mean far below lb
  cat(sprintf("%-8s mean far below lb: fn %s  gradient %s\n", fm$family, format(ob$fn(p)),
              paste(format(ob$gr(p)), collapse = " ")))
}
