# Lane sparsese, fix round 1: the cumulative fit of dev/sparsese-fix1.R at
# seed 2, whose standard-error check loses a parameter as "nonfinite" and
# whose fitted() Est.Error reached 5.3 on both routes. Prints the warnings
# and messages, the lost parameters and the largest Est.Error.
#   Rscript sparsese-cumlost.R <lib: lane|base> <seed> <limit>
a <- commandArgs(trailingOnly = TRUE)
lib <- if (a[1] == "lane") "C:/Users/adf44/source/r/wt-sparsese-lib" else {
  "C:/Users/adf44/source/r/rellib-r22"
}
seed <- as.integer(a[2])
options(frmtmb.joint_root_max = as.numeric(a[3]))
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "seed", seed, "limit", a[3], "\n")
say <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    cat("W:", substr(conditionMessage(w), 1, 300), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) {
    cat("M:", substr(conditionMessage(m), 1, 300), "\n")
    invokeRestart("muffleMessage")
  })
}
set.seed(seed)
ng <- 30L
n <- 300L
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)),
                f = factor(rep(c("a", "b"), length.out = n)))
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
d$u <- stats::rnorm(ng, 0, 0.5)[d$g]
d$e <- stats::rnorm(n, 0, 0.3)
d$z <- stats::rlogis(n)
eta <- sin(3 * d$x) + d$u + 0.3 * (d$f == "b")
d$yc <- factor(cut(0.8 * eta + d$z, c(-Inf, 0, 0.8, Inf),
                   labels = c("lo", "mid", "hi")), ordered = TRUE)
cat("-- fit\n")
fit <- say(frm(yc ~ x1 + x2 + x3 + x4 + x5 + (1 | g), data = d,
               family = cumulative(), REML = TRUE))
cat("par:", format(fit$opt$par, digits = 5), "\n")
cat("code:", fit$opt$convergence, "\n")
cat("-- sdr\n")
s <- say(frmtmb:::sdr_of(fit))
cat("lost:", paste(names(s$se_lost), s$se_lost), "\n")
cat("-- fitted\n")
fi <- say(fitted(fit))
cat("max Est.Error:", format(max(fi[, "Est.Error", ]), digits = 5),
    " rows above 0.5:", sum(fi[, "Est.Error", ] > 0.5), "of",
    length(fi[, "Est.Error", ]), "\n")
cat("-- summary\n")
print(say(summary(fit))$fixed)
