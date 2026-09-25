# Reviewer 2: what does dry_run = "objective" return (no optimization).
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
set.seed(1)
d <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
d$code <- as.integer(d$rt > 1.2); d$rt[d$code == 1] <- 1.2
o <- frm(bf(rt | dec(upper) + cens(code) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
         family = wiener(), data = d, dry_run = "objective")
print(class(o)); print(names(o))
str(o$obj$par)
print(o$obj$fn(o$obj$par))
fam <- stats::family(o)
print(names(fam)); print(fam$ndt_bound$ub)
