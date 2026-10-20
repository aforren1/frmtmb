# Lane ndaterm, fix round 2: the two refusals' texts on the lane build.
# Seed 7. Log dev/ndaterm-out/r2-msg.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ndaterm-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
set.seed(7)
d <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.3)
d$side <- factor(ifelse(d$upper == 1, "up", "down"), c("up", "down"))
ff <- frm(bf(rt | dec(side) ~ 1, bias = 0.3), family = wiener(), data = d)
fw <- frm(bf(rt | dec(upper) ~ 1, bias = 0.3), family = wiener(), data = d)
cat(tryCatch(fitted(ff, newdata = data.frame(side = "left")),
             error = conditionMessage), "\n")
cat(tryCatch(fitted(fw, newdata = data.frame(upper = 0.5)),
             error = conditionMessage), "\n")
