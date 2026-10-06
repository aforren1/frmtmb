# Reviewer, item 2: how often a cs() fit on the cumulative families
# dies in the optimizer, against the same model on sratio() (whose
# thresholds cannot cross). 20 seeds each (1..20), n = 300, cs() effect
# null or small. Also the (cs(1) | g) spelling.
#   Rscript dev/relrev-csfail.R > dev/relrev-log/csfail.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
mk <- function(seed, n = 300, mix = FALSE) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  lat <- if (mix) {
    cls <- rbinom(n, 1, 0.4)
    ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  } else 0.8 * d$x + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
one <- function(expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }), error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(r, "err")) return(if (grepl("NA/NaN gradient", r)) "nan_grad" else "error")
  ow <- w[!grepl("Category specific effects for this family", w, fixed = TRUE)]
  if (r$opt$convergence != 0) return("nonconv")
  if (length(ow)) return("warn")
  "ok"
}
designs <- list(
  cum = function(d) frm(y ~ x + cs(z), family = cumulative(), data = d),
  cum_probit = function(d) frm(y ~ x + cs(z), family = cumulative("probit"), data = d),
  sratio = function(d) frm(y ~ x + cs(z), family = sratio(), data = d),
  mix_cum_sratio = function(d) frm(bf(y ~ x + cs(z)), family = mixture(cumulative(), sratio()), data = d),
  mix_cum_cum = function(d) frm(bf(y ~ x + cs(z)), family = mixture(cumulative(), cumulative()), data = d),
  mix_sratio_sratio = function(d) frm(bf(y ~ x + cs(z)), family = mixture(sratio(), sratio()), data = d),
  mix_cum_sratio_nocs = function(d) frm(bf(y ~ x + z), family = mixture(cumulative(), sratio()), data = d)
)
for (nm in names(designs)) {
  mix <- startsWith(nm, "mix")
  res <- vapply(1:20, function(s) one(designs[[nm]](mk(s, mix = mix))), "")
  cat(sprintf("%-22s %s\n", nm, paste(names(table(res)), table(res), sep = "=", collapse = " ")))
}
d <- mk(1)
d$g <- factor(sample(letters[1:6], nrow(d), TRUE))
for (fam in c("sratio", "cumulative", "poisson")) {
  r <- tryCatch({frm(y ~ x + (cs(1) | g), family = get(fam)(), data = d); "FIT"},
                error = function(e) conditionMessage(e))
  cat(sprintf("(cs(1) | g) on %-10s %s\n", fam, substr(r, 1, 120)))
  r <- tryCatch({frm(y ~ x + (cs(x) | g), family = get(fam)(), data = d); "FIT"},
                error = function(e) conditionMessage(e))
  cat(sprintf("(cs(x) | g) on %-10s %s\n", fam, substr(r, 1, 120)))
}
