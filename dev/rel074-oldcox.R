# A cox() fit saved by frmtmb 0.73.0 (rellib-r11), read on 0.74.0
# (rellib-r12): the release review's finding 3, and the release fix
# old_core_family_home() (R/extension-home.R). Before the fix update()
# failed with "could not find function \"bhaz_spec\"".
#   Rscript dev/rel074-oldcox.R save <rds>
#   Rscript dev/rel074-oldcox.R read <rds> <attach frmtmb.survival: yes|no>
a <- commandArgs(trailingOnly = TRUE)
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
if (a[1] == "save") {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r11", user))
  suppressMessages(library(frmtmb))
  set.seed(20261023)
  d <- data.frame(x = stats::rnorm(150))
  d$time <- stats::rexp(150, exp(-0.5 + 0.7 * d$x))
  d$out <- stats::rbinom(150, 1, 0.3)
  fit <- frm(bf(time | cens(out) ~ x), family = cox(df = 5), data = d)
  cat("saved on frmtmb", format(packageVersion("frmtmb")),
      dirname(find.package("frmtmb")), "logLik",
      format(as.numeric(logLik(fit)), digits = 12), "\n")
  saveRDS(fit, a[2])
  quit(save = "no")
}
.libPaths(c("C:/Users/adf44/source/r/rellib-r12", user))
suppressMessages(library(frmtmb))
if (identical(a[3], "yes")) suppressMessages(library(frmtmb.survival))
cat("read on frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "| frmtmb.survival attached:",
    "package:frmtmb.survival" %in% search(), "\n")
fit <- suppressWarnings(readRDS(a[2]))
calls <- list(print = quote(utils::capture.output(print(fit))),
              summary = quote(summary(fit)), logLik = quote(logLik(fit)),
              vcov = quote(vcov(fit)), fitted = quote(fitted(fit)),
              update = quote(update(fit, . ~ 1)),
              simulate = quote(simulate(fit, nsim = 1)))
for (nm in names(calls)) {
  r <- tryCatch({ eval(calls[[nm]]); "ok" }, error = function(e) {
    paste(if (inherits(e, "frmtmb_extension_missing")) "REFUSED[ext]" else
      "ERR", substr(gsub("\\s+", " ", conditionMessage(e)), 1, 150))
  })
  cat(sprintf("  %-9s %s\n", nm, r))
}
