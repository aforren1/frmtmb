# Reviewer of lane optstart: a partial `start` that names mu's slope
# keeps sigma at the least-squares residual spread, which no longer
# matches the residuals at the user's slope. Does the fit still reach the
# maximum? 20 seeds per case; the default fit is the reference.
#   Rscript dev/optstart-rev-partial.R <lib> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat(a[2], "lib", dirname(find.package("frmtmb")), "\n")
one <- function(fo, d, fam, start = NULL) {
  w <- 0L
  f <- tryCatch(withCallingHandlers(suppressMessages(
    frm(fo, data = d, family = fam, start = start)),
    warning = function(x) {
      w <<- w + 1L
      invokeRestart("muffleWarning")
    }), error = function(e) NULL)
  c(ll = if (is.null(f)) NA else as.numeric(logLik(f)),
    code = if (is.null(f)) NA else f$opt$convergence, w = w)
}
rows <- list()
for (s in 1:20) {
  set.seed(s)
  n <- 400
  d <- data.frame(x = rnorm(n), g = factor(rep(1:20, 20)))
  d$y <- 5 + 1e3 * d$x + rnorm(20, 0, 2)[d$g] + rnorm(n)
  d$yt <- 5 + 1e3 * d$x + rnorm(20, 0, 2)[d$g] + rt(n, 4)
  for (cs in list(list("gaussian", y ~ x + (1 | g), gaussian()),
                  list("student", yt ~ x + (1 | g), student()))) {
    ref <- one(cs[[2]], d, cs[[3]])
    for (sl in c(0, 500, 2000)) {
      v <- one(cs[[2]], d, cs[[3]], list(beta = c(x = sl)))
      rows[[length(rows) + 1L]] <- data.frame(seed = s, fam = cs[[1]],
        slope_start = sl, gap = v[["ll"]] - ref[["ll"]], code = v[["code"]],
        warned = v[["w"]] > 0, ref_code = ref[["code"]])
    }
  }
}
X <- do.call(rbind, rows)
print(aggregate(cbind(error = is.na(gap), short = is.na(gap) | gap < -1e-4,
                      silent_short = !is.na(gap) & gap < -1e-4 &
                        code == 0 & !warned,
                      warned = warned) ~ fam + slope_start, data = X,
                FUN = sum, na.action = na.pass))
print(summary(X$gap))
