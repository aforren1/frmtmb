# Reviewer of lane optrap: a copy of dev/optrap-scale.R writing to dev/optrap-rev-log/, with OPTRAP_REV_RESID=1 sourcing dev/optrap-rev-resid.R.
# Lane optrap, item 1: does a gaussian mixed fit reach lme4's maximum at
# every scale of the response, and is the fit of y * c the fit of y
# scaled? Designs: "ri20" y ~ x + (1 | g), 20 groups of 5; "slope"
# y ~ t + (1 + t | s), 18 subjects of 10 days (sleepstudy-like, t in
# 0..9, intercept sd 25, slope sd 6, sigma 25, then divided by 25).
# Scales 1e-6 to 1e6. Per fit: code, logLik - lme4's at the same scale,
# and the gap to the same seed's fit at scale 1 after unscaling.
#   Rscript dev/optrap-scale.R <lib or ""> <tag> [seeds]
args <- commandArgs(TRUE)
.libPaths(c(if (nzchar(args[1])) args[1],
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
tag <- args[2]
seeds <- if (length(args) > 2) eval(parse(text = args[3])) else 1:40
suppressMessages({library(frmtmb); library(lme4)})
# ablations of the lane's three parts, on its own build
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
  utils::assignInNamespace("predictor_par_units",
                           function(...) NULL, "frmtmb")
}
if (nzchar(Sys.getenv("OPTRAP_REV_RESID"))) source("dev/optrap-rev-resid.R")
cat("ablated:", ab, "\n")
cat("lib", find.package("frmtmb"), "\n")
q <- function(expr) suppressMessages(suppressWarnings(expr))
scales <- 10^c(-6, -3, 0, 3, 6)
mk <- function(design, s) {
  set.seed(s)
  if (design == "ri20") {
    # the setier construction: no group effect in the truth
    d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
    d$y <- 1 + 0.5 * d$x + rnorm(100)
    list(d = d, f = y ~ x + (1 | g))
  } else {
    ns <- 18
    d <- data.frame(s = factor(rep(seq_len(ns), each = 10)),
                    t = rep(0:9, ns))
    b0 <- rnorm(ns, 0, 25)
    b1 <- rnorm(ns, 0, 6)
    d$y <- (250 + b0[d$s] + (10 + b1[d$s]) * d$t + rnorm(nrow(d), 0, 25)) /
      25
    list(d = d, f = y ~ t + (1 + t | s))
  }
}
X <- list()
for (design in c("ri20", "slope")) for (s in seeds) {
  m <- mk(design, s)
  ref <- NULL
  for (sc in scales) {
    d <- m$d
    d$y <- d$y * sc
    l4 <- q(lmer(m$f, data = d, REML = FALSE,
                 control = lmerControl(check.conv.singular = "ignore")))
    f <- q(frm(m$f, data = d))
    ll <- as.numeric(logLik(f))
    th <- f$estimates$theta
    sg <- f$estimates$betad[1]
    n <- nrow(d)
    row <- data.frame(design = design, seed = s, scale = sc,
                      code = f$opt$convergence,
                      dll4 = ll - as.numeric(logLik(l4)),
                      ll_unscaled = ll + n * log(sc),
                      th1 = th[1] - log(sc), sig = sg - log(sc),
                      b1 = f$estimates$beta[2] / sc)
    X[[length(X) + 1L]] <- row
  }
}
X <- do.call(rbind, X)
# gap to the same seed's scale-1 fit, after unscaling
key <- paste(X$design, X$seed)
one <- X[X$scale == 1, ]
i1 <- match(key, paste(one$design, one$seed))
X$dll1 <- X$ll_unscaled - one$ll_unscaled[i1]
X$dth1 <- X$th1 - one$th1[i1]
dir.create("dev/optrap-rev-log", showWarnings = FALSE)
write.table(X, sprintf("dev/optrap-rev-log/scalecopy-%s.tsv", tag), sep = "\t",
            quote = FALSE, row.names = FALSE)
for (design in unique(X$design)) for (sc in scales) {
  x <- X[X$design == design & X$scale == sc, ]
  cat(sprintf(paste0("%-5s scale %-6g | below lme4 by > 1e-4: %2d ",
                     "(min %.3g) | code != 0: %2d | max |logLik - ",
                     "scale-1 fit| %.3g | max |dlog sd| %.3g\n"),
              design, sc, sum(x$dll4 < -1e-4), min(x$dll4),
              sum(x$code != 0), max(abs(x$dll1)), max(abs(x$dth1))))
}
