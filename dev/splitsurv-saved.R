# Saved cox() and royston_parmar() fits read back in a session where
# frmtmb.survival is not installed (review finding 2). The scratch
# library holds only the lane's frmtmb, ahead of rellib-r11, which has
# no frmtmb.survival.
#   Rscript dev/splitsurv-saved.R save <dir>   (lane lib, writes fits)
#   Rscript dev/splitsurv-saved.R absent <dir> <scratch lib>
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
O <- a[2]
first <- if (arm == "absent") a[3] else
  "C:/Users/adf44/source/r/wt-splitsurv-lib"
.libPaths(c(first, "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "; survival installed:",
    nzchar(system.file(package = "frmtmb.survival")), "\n")
if (arm == "save") {
  suppressMessages(library(frmtmb.survival))
  e <- new.env()
  utils::data("bc", package = "flexsurv", envir = e)
  bc <- e$bc
  bc$censored <- 1 - bc$censrec
  fr <- frm(bf(recyrs | cens(censored) ~ group),
            family = royston_parmar(df = 2), data = bc)
  set.seed(1)
  dd <- data.frame(x = rnorm(150), z = rnorm(150))
  dd$time <- rexp(150, exp(-0.5 + 0.7 * dd$x))
  dd$out <- rbinom(150, 1, 0.3)
  fc <- frm(bf(time | cens(out) ~ x), family = cox(df = 6), data = dd)
  saveRDS(list(fr = fr, fc = fc, bc = bc, dd = dd),
          file.path(O, "saved-fits.rds"))
  cat("saved\n")
  quit(save = "no")
}
s <- withCallingHandlers(readRDS(file.path(O, "saved-fits.rds")),
                         warning = function(w) {
                           cat("readRDS WARNING:", conditionMessage(w), "\n")
                           invokeRestart("muffleWarning")
                         })
show <- function(lab, expr) {
  r <- tryCatch({ force(expr); "OK" }, error = function(e) {
    paste0("ERROR [", class(e)[1], "]: ", conditionMessage(e))
  })
  cat(sprintf("%-24s %s\n", lab, substr(r, 1, 170)))
}
g <- data.frame(recyrs = seq(0.5, 6, length.out = 5),
                group = factor("Good", levels = levels(s$bc$group)))
for (nm in c("fc", "fr")) {
  f <- s[[nm]]
  cat("==", nm, "\n")
  show("print", utils::capture.output(print(f)))
  show("summary", utils::capture.output(print(summary(f))))
  show("logLik", logLik(f))
  show("vcov", vcov(f))
  show("fixef", fixef(f))
  show("predict(link)", predict(f, type = "link"))
  show("frm_linpred", frm_linpred(f, type = "link"))
  show("fitted", fitted(f))
  show("residuals", residuals(f))
  show("simulate", simulate(f, nsim = 1))
  show("update", update(f, . ~ 1))
  show("confint", confint(f))
  show("family", family(f))
  if (nm == "fr") {
    show("frm_curve", frm_curve(f, newdata = g, dpar = "mu",
                                simultaneous = FALSE))
  }
}
