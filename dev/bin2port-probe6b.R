# Probe 6b: on the same wiener fit as probe 6 (seed 420), fitted() with
# newdata reads the boundary from newdata and predict() does not.
# Log: dev/bin2port-out/probe6b.txt.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
set.seed(420)
dat <- ddm_simulate(200, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
fit <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5), family = wiener(),
           data = dat)
nd <- dat[1:2, , drop = FALSE]
r <- tryCatch(fitted(fit, newdata = nd), error = conditionMessage)
cat("fitted(fit, newdata = 2 rows):\n")
print(r)
r <- tryCatch(predict(fit, newdata = nd), error = conditionMessage)
cat("predict(fit, newdata = 2 rows):", if (is.character(r)) r else "ran",
    "\n")
