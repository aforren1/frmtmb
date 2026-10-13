# Reviewer: mo_search()'s cost on a model with three simplexes from a
# factor (mo(x) * f), interleaved arms, minimum of 8 rounds, with a
# control arm that repeats the mo_search = FALSE fit.
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
n <- 600
d <- data.frame(x = sample(0:3, n, TRUE), f = factor(sample(c("a", "b", "c"), n, TRUE)),
                g = factor(sample(1:30, n, TRUE)))
d$y <- 1 + c(a = 0.4, b = -0.6, c = 0.02)[as.character(d$f)] *
  c(0, 0, 1, 1)[d$x + 1] + rnorm(30, 0, 0.4)[d$g] + rnorm(n, 0, 0.5)
fitit <- function(ms) suppressWarnings(frm(bf(y ~ mo(x) * f + (1 | g)),
  family = gaussian(), data = d, control = frmtmb_control(mo_search = ms)))
tm <- function(ms) { t0 <- proc.time()[[3]]; for (i in 1:5) f <- fitit(ms); list(t = proc.time()[[3]] - t0, f = f) }
res <- list(on = Inf, off = Inf, ctl = Inf)
for (r in 1:8) {
  a <- tm(TRUE); b <- tm(FALSE); c <- tm(FALSE)
  res$on <- min(res$on, a$t); res$off <- min(res$off, b$t); res$ctl <- min(res$ctl, c$t)
}
cat(sprintf("5 fits each, min of 8 rounds: mo_search on %.2f s, off %.2f s, control (off again) %.2f s | ratio on/off %.2f, control/off %.2f\n",
            res$on, res$off, res$ctl, res$on / res$off, res$ctl / res$off))
cat("mo_search slot:", a$f$opt$mo_search, "| logLik on", format(as.numeric(logLik(a$f)), digits = 12),
    "off", format(as.numeric(logLik(b$f)), digits = 12), "\n")
cat("weights:", vapply(frmtmb::mo_frame_terms(a$f), function(t) paste(round(frmtmb::mo_simplex(a$f$estimates[[t$zeta]]), 4), collapse = ","), ""), "\n")
