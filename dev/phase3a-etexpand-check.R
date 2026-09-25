# Which of dev/phase3a-ode-rxode2.R's random schedules does
# rxode2::etExpand() fail on, and do they all carry an evid 3 record?
source("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-lib.R")
suppressMessages(library(rxode2))
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-ode-rxode2.R")
a <- grep("^rand_sched <- function", src)
b <- a - 1L + which(src[a:length(src)] == "}")[1L]
eval(parse(text = src[a:b]))
res <- t(vapply(1:200, function(k) {
  e <- rand_sched(20260923L + k)
  fails <- inherits(try(rxode2::etExpand(e), silent = TRUE), "try-error")
  c(fails = fails, reset = any(as.data.frame(e)$evid == 3))
}, logical(2)))
print(table(etExpand_fails = res[, "fails"], has_evid3 = res[, "reset"]))
