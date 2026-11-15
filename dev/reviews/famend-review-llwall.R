# Review of lane famend, attack 1: the fit of test-linpred-layers.R
# ("several location predictors stack, and cs() adds its thresholds") on
# base and lane: code, logLik, warnings, the narrowest own category, and
# whether other starts reach a higher point off the wall.
#   Rscript dev/reviews/famend-review-llwall.R <base|lane>
a <- commandArgs(trailingOnly = TRUE)
lanelib <- "C:/Users/adf44/source/r/wt-famend-lib"
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] == "lane") lanelib, base, ulib))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "\n")
cap <- function(expr) {
  w <- character(); m <- character()
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
    }), error = function(e) structure(conditionMessage(e), class = "err"))
  list(v = v, w = w, m = m)
}
ll_data <- function() {
  set.seed(13)
  n <- 240
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$y <- cut(d$x - 0.5 * d$z + stats::rlogis(n),
             c(-Inf, -0.8, 0.2, 1.2, Inf), labels = FALSE)
  d
}
d <- ll_data()
fo <- bf(y ~ 1, mu1 ~ x + cs(z), mu2 ~ x)
mixfam <- mixture(cumulative(), cumulative())
r <- cap(frm(fo, family = mixfam, data = d,
             control = frmtmb_control(check_se = "ignore")))
f <- r$v
cat(sprintf("test fit: code %d message '%s' logLik %.10f\n",
            f$opt$convergence, f$opt$message, as.numeric(logLik(f))))
for (x in r$w) cat("   W:", substr(gsub("\\s+", " ", x), 1, 400), "\n")
for (x in r$m) cat("   M:", substr(gsub("\\s+", " ", x), 1, 200), "\n")
dg <- ns$mixture_ord_degeneracy(f, f$spec$responses[[1]]$resp_name)
print(dg)
cat("max |grad|:", signif(max(abs(f$obj$gr(f$opt$par))), 4), "\n")
cat("par:", paste(names(f$opt$par), signif(f$opt$par, 5), collapse = "; "), "\n")
# single cumulative with cs(z), for scale
r1 <- cap(frm(bf(y ~ x + cs(z)), family = cumulative(), data = d))
cat(sprintf("single cumulative cs(z): code %d logLik %.6f\n",
            r1$v$opt$convergence, as.numeric(logLik(r1$v))))
# other starts: jitter the reported optimum and refit from each
set.seed(99)
best <- as.numeric(logLik(f))
for (k in 1:12) {
  st <- f$obj$env$parList(f$opt$par)
  st <- st[intersect(names(st), c("beta", "betad", "tau_raw1", "tau_raw2"))]
  st <- lapply(st, function(v) v + rnorm(length(v), 0, 0.7))
  rk <- cap(frm(fo, family = mixfam, data = d, start = st,
                control = frmtmb_control(check_se = "ignore")))
  if (inherits(rk$v, "err")) {
    cat("start", k, "ERR", substr(rk$v, 1, 100), "\n")
    next
  }
  dk <- ns$mixture_ord_degeneracy(rk$v, rk$v$spec$responses[[1]]$resp_name)
  cat(sprintf("start %2d: code %d logLik %.6f closed %s warns %d\n", k,
              rk$v$opt$convergence, as.numeric(logLik(rk$v)),
              if (!is.null(dk$closed)) paste(signif(dk$closed, 3), collapse = ",")
              else "-", length(rk$w)))
}
cat("\ndone\n")
