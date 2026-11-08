# Release review 0.80.0: two cross-lane cases.
# A. fitedge x bigx: an exact gp() in a body's parameter beside a body
#    term linear in a shifted column, `a + b * (z + m)` with a ~ gp(x). The
#    shift is absorbed by a's intercept, so every m is one model.
# B. osa x latperf: residuals(type = "osa") on a latent residual
#    (poisson, ar(cov = TRUE)), on this build and on the osa lane's build;
#    latperf changed no objective, so the residuals should agree.
#
#   Rscript dev/review-rel-cross.R <lib> [A|B]
a <- commandArgs(trailingOnly = TRUE)
lib <- a[1]
part <- if (length(a) > 1L) a[2] else "AB"
.libPaths(c(strsplit(lib, ";")[[1]], "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
cnd <- function(expr) {
  cs <- character()
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(w) {
      cs <<- c(cs, paste0("W:", class(w)[1L]))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cs <<- c(cs, paste0("M:", class(m)[1L]))
      invokeRestart("muffleMessage")
    }), error = function(e) {
      cs <<- c(cs, paste0("E:", substr(conditionMessage(e), 1, 100)))
      NULL
    })
  list(v = v, cs = unique(cs))
}
if (grepl("A", part)) for (seed in 1:3) {
  set.seed(9100 + seed)
  n <- 120
  d <- data.frame(x = sort(runif(n, 0, 10)), z = rnorm(n))
  d$y <- 2 * sin(d$x) + 0.8 * d$z + rnorm(n, 0, 0.4)
  base <- NULL
  for (m in c(0, 1e4, 1e6)) {
    d$w <- d$z + m
    r <- cnd(frm(bf(y ~ a + b * w, a ~ gp(x), b ~ 1, nl = TRUE), data = d))
    ll <- if (is.null(r$v)) NA else as.numeric(logLik(r$v))
    if (m == 0) base <- ll
    b <- if (is.null(r$v)) NA else r$v$estimates$beta[["b_(Intercept)"]]
    cat(sprintf("A seed %d m %.0e ll %.4f (centered %.4f, short %.4f) b %.5f | %s\n",
                seed, m, ll, base, base - ll, b,
                if (length(r$cs)) paste(r$cs, collapse = " ") else "-"))
  }
}
if (grepl("B", part)) {
  set.seed(5521)
  n <- 300
  d <- data.frame(t = 1:n, x = rnorm(n))
  e <- as.numeric(stats::arima.sim(list(ar = 0.6), n, sd = 0.4))
  d$y <- rpois(n, exp(0.5 + 0.3 * d$x + e))
  f <- cnd(frm(y ~ x + ar(time = t, cov = TRUE), data = d, family = poisson()))
  t0 <- proc.time()[["elapsed"]]
  r <- cnd(residuals(f$v, type = "osa"))
  cat(sprintf("B fit %s | osa %d rows, finite %d, %.1f s | %s\n",
              paste(f$cs, collapse = " "), length(r$v), sum(is.finite(r$v)),
              proc.time()[["elapsed"]] - t0, paste(r$cs, collapse = " ")))
  out <- file.path("dev/review-rel-out",
                   paste0("osa-", basename(strsplit(lib, ";")[[1]][1]), ".rds"))
  saveRDS(r$v, out)
  cat("B saved", out, "\n")
}
cat("CROSS DONE\n")
