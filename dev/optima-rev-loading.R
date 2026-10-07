# Reviewer of lane optima, claim 7: a non-mo() construction of a removed
# Hessian direction (flat or concave) that loads a little on a parameter
# which keeps its standard error. se_tier3() is traced; for every fit
# the largest such loading is logged. Families whose simplex is still a
# softmax (cox()'s baseline) and mixtures' mixing weights are scanned.
#   Rscript dev/optima-rev-loading.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
rec <- new.env()
suppressMessages(trace("se_tier3", where = ns, print = FALSE,
  exit = quote({
    r <- get("rec", envir = globalenv())
    r$lost <- why
    r$small <- 0
    r$kind <- ""
    if (exists("rem", inherits = FALSE) && length(rem)) {
      for (k in rem) {
        v <- abs(e$vectors[, k])
        kept <- mark == ""
        if (any(kept) && any(!kept)) {
          s <- max(v[kept])
          if (s > r$small) {
            r$small <- s
            r$kind <- if (ev[k] < 0) "neg" else "pos"
            r$who <- nm[free][kept][which.max(v[kept])]
          }
        }
      }
    }
  })))
assign("rec", rec, envir = globalenv())
scan1 <- function(label, f) {
  rec$small <- 0
  rec$kind <- ""
  rec$who <- ""
  fit <- tryCatch(suppressWarnings(f()), error = function(e) NULL)
  if (is.null(fit)) return(cat(label, "ERROR\n"))
  sl <- tryCatch(frmtmb:::sdr_of(fit)$se_lost, error = function(e) NULL)
  if (!length(sl)) return(invisible())
  cat(sprintf("%-22s code %d lost %s | max kept loading %.3g (%s, %s)\n",
              label, fit$opt$convergence,
              paste(names(sl), sl, sep = "=", collapse = ","),
              rec$small, rec$kind, rec$who))
}
for (s in 1:30) {
  set.seed(s)
  n <- 300
  x <- rnorm(n)
  tt <- stats::rexp(n, exp(-0.5 + 0.7 * x))
  ct <- stats::rexp(n, 0.2)
  dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
  scan1(paste("cox seed", s), function()
    frm(bf(time | cens(cens) ~ x), family = cox(), data = dd))
  set.seed(s)
  dm <- data.frame(x = rnorm(200))
  dm$y <- 1 + 0.5 * dm$x + rnorm(200)
  scan1(paste("mix2 seed", s), function()
    frm(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = dm))
}
