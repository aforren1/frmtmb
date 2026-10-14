# Is the degenerate mixture fit in dev/famuni-out/wiener-links.txt
# (logLik 7713.9) the unbounded ndt link's doing, or the documented
# mixture recipe's on data with no contaminant? Same data, seed 20261013,
# plus a data set that has 5 percent fast guesses.
# Output: dev/famuni-out/wiener-mix.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all(
  "C:/Users/adf44/source/r/frmtmb-wt-famuni/extensions/frmtmb.eam",
  quiet = TRUE))
seed <- 20261013L
set.seed(seed)
n <- 400
d <- ddm_simulate(n, mu = 0.8, bs = 1.3, ndt = 0.3, bias = 0.45)
d$x <- rnorm(n)
g <- d
k <- sample(n, 20)
g$rt[k] <- runif(20, 0.1, 0.3)
run <- function(lab, fam, data) {
  w <- character(0)
  f <- withCallingHandlers(
    tryCatch(frm(bf(rt | dec(upper) ~ 1), family = fam, data = data),
             error = function(e) e),
    warning = function(cw) {
      w <<- c(w, substr(conditionMessage(cw), 1, 60))
      invokeRestart("muffleWarning")
    })
  if (inherits(f, "error")) {
    cat(sprintf("%-28s ERROR %s\n", lab, conditionMessage(f))); return()
  }
  ndt <- frm_linpred(f, dpar = "ndt1", type = "response")[1]
  th <- frm_linpred(f, dpar = "theta1", type = "response")[1]
  cat(sprintf("%-28s conv %d logLik %10.3f ndt1 %.4g theta1 %.4g warn %d\n",
              lab, f$opt$convergence, as.numeric(logLik(f)), ndt, th,
              length(w)))
  if (length(w)) cat("   ", unique(w), sep = "\n    ")
}
for (dd in list(list("no guesses", d), list("20 fast guesses", g))) {
  cat("\n== data:", dd[[1]], " min rt", min(dd[[2]]$rt), "\n")
  run("max_ndt = 0.4 (documented)",
      frmtmb::mixture(wiener(max_ndt = 0.4, allow_unreachable = TRUE),
                      lognormal()), dd[[2]])
  run("link_ndt = log",
      frmtmb::mixture(wiener(link_ndt = "log", allow_unreachable = TRUE),
                      lognormal()), dd[[2]])
}
