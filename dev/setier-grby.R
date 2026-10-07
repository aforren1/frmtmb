# Lane setier, punch round 1: the gr(g, by = f) fixture of
# test-se-check.R (seed 11): its thetas, what moving each toward its
# end does to the objective, and the verdict.
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(11)
dd <- data.frame(x = stats::rnorm(160), g = factor(rep(1:16, 10)))
dd$f <- factor(ifelse(as.integer(dd$g) <= 8, "a", "b"))
u <- cbind(stats::rnorm(16, 0, 0.7), stats::rnorm(16, 0, 0.4))
dd$y <- stats::rnorm(160, 1 + 0.5 * dd$x + u[dd$g, 1] +
                       u[dd$g, 2] * dd$x, 1)
w <- character()
fit <- withCallingHandlers(frm(bf(y ~ x + (1 + x | gr(g, by = f))),
                               family = gaussian(), data = dd),
  warning = function(c) {w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")},
  message = function(c) {w <<- c(w, paste("M:", conditionMessage(c)))
    invokeRestart("muffleMessage")})
cat(substr(w, 1, 200), sep = "\n")
nm <- ns$outer_par_names(fit)
p <- fit$opt$par
cat("par", paste0(nm, "=", signif(p, 4)), "\n")
print(ns$sdr_of(fit)$se_lost)
f0 <- fit$obj$fn(p)
for (k in grep("theta", nm)) for (s in c(-2, -0.5, 0.5, 2)) {
  q <- p
  q[k] <- q[k] + s
  cat(nm[k], s, signif(fit$obj$fn(q) - f0, 3), "\n")
}
