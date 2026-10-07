# Lane optima, item 4a: does frmtmb_control(restarts = ) take
# fit_loss2 to the best code-0 optimum of dev/optima-loss2.R?
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
loss <- utils::read.csv("dev/optima-log/ClarkTriangle.csv")
nlform2 <- bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
              ult ~ 1 + (1 | ID1 | AY), omega ~ 1 + (1 | ID1 | AY),
              theta ~ 1 + (1 | ID1 | AY), nl = TRUE)
for (r in c(1, 3, 10)) {
  f <- suppressWarnings(frm(nlform2, data = loss, family = gaussian(),
                            start = list(beta = c(5000, 1, 45)),
                            control = frmtmb_control(restarts = r)))
  cat(sprintf("restarts %2d: logLik %.6f code %d evals %d\n", r,
              as.numeric(logLik(f)), f$opt$convergence, f$opt$evals))
}
for (o in c("optim")) {
  f <- suppressWarnings(frm(nlform2, data = loss, family = gaussian(),
                            start = list(beta = c(5000, 1, 45)),
                            control = frmtmb_control(optimizer = o)))
  cat(sprintf("optimizer %s: logLik %.6f code %d\n", o,
              as.numeric(logLik(f)), f$opt$convergence))
}
