wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
b <- readRDS(file.path(wt, "dev/ordinal-rev-bitwise-base.rds"))
l <- readRDS(file.path(wt, "dev/ordinal-rev-bitwise-lane.rds"))
x <- b$mv_ord_gauss$coef; y <- l$mv_ord_gauss$coef
str(attributes(x)); str(attributes(y)); print(all.equal(x, y))
x <- b$cum_cloglog_gr$fixef; y <- l$cum_cloglog_gr$fixef
print(abs(x - y) / abs(x))
