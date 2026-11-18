# rel084-kt copy of dev/reviews/sparsese-review-cum4500.R: arm lane = rellib-r23 then user lib (lane lib and rellib-r22 removed). No other change.
# Review of lane sparsese: a cumulative REML fit on a degree-5 raw
# polynomial with (1 | g), 4500 groups of 2 rows (4507 joint parameters,
# above joint_root_max), fitted() of 6 new rows.
#   Rscript sparsese-review-cum4500.R <base|lane> <limit|default>
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
lim <- a[2]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "lib", dirname(find.package("frmtmb")), "limit", lim, "\n")
if (lim != "default") options(frmtmb.joint_root_max = as.numeric(lim))
set.seed(1)
ng <- 4500L
n <- 2L * ng
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)))
for (k in 1:5) d[[paste0("x", k)]] <- d$x^k
eta <- sin(3 * d$x) + stats::rnorm(ng, 0, 0.5)[d$g]
z <- 0.8 * eta + stats::rlogis(n)
d$y <- factor(cut(z, c(-Inf, 0, 0.8, Inf), labels = c("lo", "mid", "hi")),
              ordered = TRUE)
w <- character()
fit <- withCallingHandlers(suppressMessages(frm(
  y ~ x1 + x2 + x3 + x4 + x5 + (1 | g), data = d, family = cumulative(),
  REML = TRUE)), warning = function(x) {
    w <<- c(w, substr(conditionMessage(x), 1, 160))
    invokeRestart("muffleWarning")
  })
sdr <- frmtmb:::sdr_of(fit)
cat("whitened", !is.null(fit$beta_whiten), "joint", nrow(sdr$jointPrecision),
    "lost [", sdr$se_lost, "] rows attr",
    length(attr(sdr$joint_root, "rows")), "\n")
nd <- d[1:6, ]
t0 <- proc.time()[["elapsed"]]
fp <- withCallingHandlers(fitted(fit, newdata = nd), warning = function(x) {
  w <<- c(w, substr(conditionMessage(x), 1, 160))
  invokeRestart("muffleWarning")
})
cat(sprintf("fitted(newdata) %.2f s\n", proc.time()[["elapsed"]] - t0))
print(round(fp[, "Est.Error", ], 6))
lp <- suppressWarnings(frm_linpred(fit, newdata = nd, se.fit = TRUE))
cat("linpred se.fit:", format(lp$se.fit, digits = 7), "\n")
if (length(w)) cat("W:", paste(unique(w), collapse = " | "), "\n")
