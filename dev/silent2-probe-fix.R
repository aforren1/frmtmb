# Lane silent2 fix round: F3 (cs_sr_acat seed 1's warnings) and 5214 (1)
# (a random-effects ordinal mixture that ends degenerate).
#   Rscript dev/silent2-probe-fix.R
.libPaths(c("C:/Users/adf44/source/r/wt-silent2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
gen <- function(seed, n = 300, K = 4) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  e <- rlogis(n)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + e
  cuts <- stats::quantile(lat, seq_len(K - 1) / K)
  d$y <- 1L + rowSums(outer(lat, cuts, ">"))
  d
}
ws <- function(e) {
  w <- character()
  withCallingHandlers(e, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  w
}
w <- ws(f <- frm(bf(y ~ cs(x) + z), family = mixture(sratio(), acat()),
                 data = gen(1)))
cat("cs_sr_acat seed 1 logLik", as.numeric(logLik(f)), "\n")
for (x in w) cat(" W:", substr(gsub("\\s+", " ", x), 1, 400), "\n")
mk <- function(seed, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
for (s in c(9, 12, 16, 5)) {
  for (fam in c("cum_cum", "cum_sratio")) {
    d <- mk(s)
    d$g <- factor(rep_len(1:20, 300))
    fm <- if (fam == "cum_cum") mixture(cumulative(), cumulative()) else
      mixture(cumulative(), sratio())
    t0 <- proc.time()[["elapsed"]]
    w <- tryCatch(ws(fr <- frm(bf(y ~ x + z + (1 | g)), family = fm,
                               data = d)), error = function(e) e)
    if (inherits(w, "error")) {
      cat("RE", fam, "seed", s, "ERROR", substr(conditionMessage(w), 1, 80),
          "\n")
      next
    }
    cat(sprintf("RE %s seed %d %.1fs code %d np %d an %s : %s\n", fam, s,
                proc.time()[["elapsed"]] - t0, fr$opt$convergence,
                length(fr$opt$par), !is.null(fr$cache$se_fit_an),
                paste(substr(gsub("\\s+", " ", w), 1, 200), collapse = " || ")))
  }
}
