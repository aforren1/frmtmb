# summary() of the fuzz cumulative spec under the other modes.
#   Rscript dev/fuzzfix-cum-modes.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
d <- readRDS("dev/fuzzfix-out/data-cum.rds")
fo <- bf(y | weights(w) ~ 1 + x + (1 | ga/gb))
for (fam in c("cumulative", "sratio", "cratio", "acat")) {
  for (m in c("ml", "profile", "autoscale")) {
    ctl <- switch(m, ml = frmtmb_control(),
                  profile = frmtmb_control(profile = TRUE),
                  autoscale = frmtmb_control(autoscale = TRUE))
    f <- tryCatch(suppressMessages(suppressWarnings(
      frm(fo, family = get(fam)(), data = d, control = ctl))),
      error = function(e) e)
    s <- if (inherits(f, "error")) f else
      tryCatch(suppressMessages(suppressWarnings(summary(f))),
               error = function(e) e)
    cat(sprintf("MODE %-10s %-9s %s\n", fam, m,
                if (inherits(s, "error")) paste("ERROR", conditionMessage(s))
                else "summary ok"))
  }
}
