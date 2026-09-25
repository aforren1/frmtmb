# Reviewer 2, item 8: gddm() with cens() and trunc(), with a valid
# condition index so that the family's own refusal is what fires.
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
msg <- function(expr) tryCatch({ expr; "ACCEPTED" }, error = function(e)
  paste0("REFUSED: ", substr(conditionMessage(e), 1, 160)))
set.seed(3)
dg <- ddm_simulate(200, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
dg$cond <- 1L; dg$code <- 0L; dg$code[1:5] <- 1L
cat("gddm plain:  ", msg(frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5), family = gddm(),
                            data = dg, dry_run = "objective")), "\n")
cat("gddm cens:   ", msg(frm(bf(rt | vint(upper, cond) + cens(code) ~ 1, bias = 0.5),
                            family = gddm(), data = dg, dry_run = "objective")), "\n")
cat("gddm trunc:  ", msg(frm(bf(rt | vint(upper, cond) + trunc(ub = 5) ~ 1, bias = 0.5),
                            family = gddm(), data = dg, dry_run = "objective")), "\n")
