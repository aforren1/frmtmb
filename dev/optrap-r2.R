# Lane optrap, punch round 1: the review's high-R^2 design
# multipliers, with one of the lane's three parts taken out
# (OPTRAP_ABLATE as in dev/optrap-scale.R: inner, theta, units), and the
# first warning of each fit.
#   Rscript dev/optrap-rev-r2b.R <lib or none> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
lab <- a[2]
ab <- Sys.getenv("OPTRAP_ABLATE")
if (grepl("inner", ab)) {
  utils::assignInNamespace("predictor_inner_control",
                           function(units) list(maxit = 1000), "frmtmb")
}
if (grepl("theta", ab)) {
  utils::assignInNamespace("unit_theta_start",
                           function(frame, tpl, units = NULL) tpl, "frmtmb")
}
if (grepl("units", ab)) {
  utils::assignInNamespace("predictor_par_units", function(...) NULL,
                           "frmtmb")
}
if (nzchar(Sys.getenv("OPTRAP_ROUND0"))) source("dev/optrap-round0-unit.R")
cat(lab, "ablate", ab, "lib", dirname(find.package("frmtmb")), "\n")
for (k in c(300, 1e3, 3e3, 1e4)) {
  nw <- 0L; nc <- 0L; sh <- 0L; msgs <- character()
  for (s in 1:40) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
    d$y <- 1 + k * d$x + rnorm(20, 0, 0.5)[d$g] + rnorm(100)
    w <- character()
    f <- withCallingHandlers(suppressMessages(frm(y ~ x + (1 | g),
                                                  data = d)),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      })
    l4 <- suppressMessages(lmer(y ~ x + (1 | g), data = d, REML = FALSE,
      control = lmerControl(check.conv.singular = "ignore")))
    nw <- nw + (length(w) > 0)
    nc <- nc + (f$opt$convergence != 0)
    sh <- sh + (as.numeric(logLik(f)) - as.numeric(logLik(l4)) < -1e-4)
    if (length(w)) msgs <- c(msgs, substr(w[1], 1, 70))
    if (length(w)) cat("  warned: k", k, "seed", s, "
")
  }
  cat(sprintf("%s k %-6g | warned %d | code != 0 %d | short %d\n", lab, k,
              nw, nc, sh))
  if (length(msgs)) print(table(msgs))
}
