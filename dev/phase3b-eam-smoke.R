# Development smoke test for items 3.4 and 3.5. Not evidence.
.libPaths(c("C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all("extensions/frmtmb.eam", quiet = TRUE))
library(frmtmb)
set.seed(2)
d <- ddm_simulate(600, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
dl <- 1.2
d$cens <- ifelse(d$rt > dl, "right", "none")
d$rtc <- pmin(d$rt, dl)
cat("censored:", sum(d$cens == "right"), "\n")
fit <- frm(bf(rtc | dec(upper) + cens(cens) ~ 1, bias = 0.5),
           family = wiener(), data = d)
print(fixef(fit)); print(logLik(fit))
fc <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
          family = wiener(contaminant = TRUE), data = d)
print(fixef(fc)); print(logLik(fc))
fp <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5), family = wiener(), data = d)
print(logLik(fp))
v <- try(frm(bf(rtc | dec(upper) + cens(cens) ~ 1, bias = 0.5),
             family = wiener(variability = "sv"), data = d))
