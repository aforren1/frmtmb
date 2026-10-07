.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib", "C:/Users/adf44/source/r/rellib-r6", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
for (s in 1:6) {
set.seed(s); n <- 300; x <- rnorm(n)
tt <- stats::rexp(n, exp(-0.5 + 0.7 * x)); ct <- stats::rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
w <- character()
fit <- withCallingHandlers(frm(bf(time | cens(cens) ~ x), family = cox(), data = dd), warning = function(z) {w <<- c(w, conditionMessage(z)); invokeRestart("muffleWarning")})
sdr <- frmtmb:::sdr_of(fit)
cat("seed", s, "code", fit$opt$convergence, "minw", format(min(cox_baseline(fit)), digits=3), "lost", paste(names(sdr$se_lost), sdr$se_lost), "tier", format(sdr$se_tier), "warn", substr(paste(w, collapse=" || "), 1, 200), "\n")
}
