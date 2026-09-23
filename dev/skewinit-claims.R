# Two claims in dev/skewinit-findings.md that were not measured when
# they were written.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")

# 1. The residual-skew start value at seed 34 of the dead-draw stream,
#    which the wt-drmtmb reviewer recorded as 2.5524.
dd <- make_data(34, TRUE)
m3 <- skew(residuals(stats::lm(y ~ xs, data = dd)))
cat("seed 34 dead draw: residual skew", format(m3, digits = 8),
    " start rule", format(2 * sign(m3) + 0.5 * m3, digits = 8), "\n")

# 2. Whether diagnose() already says something about a skew_normal fit
#    whose alpha is weakly identified.
make_sym <- function(seed, n = 50) {
  set.seed(1000 + seed)
  xs <- -abs(rnorm(n)) * 3
  data.frame(y = xs + rnorm(n, 0, 1.5), xs = xs)
}
for (s in c(1, 3, 6)) {
  f <- suppressWarnings(frm_sn(make_sym(s)))
  d <- diagnose(f, quiet = TRUE)
  se <- suppressWarnings(sqrt(diag(vcov(f)))[["alpha_Intercept"]])
  cat(sprintf("\nsym seed %d: alpha %12.4g  se %12.4g\n", s,
              f$estimates$betad[["alpha_(Intercept)"]], se))
  for (nm in names(d)) {
    v <- d[[nm]]
    if (is.null(v) || (is.character(v) && !length(v))) next
    cat("  ", nm, ": ", paste(utils::head(format(v), 3), collapse = " | "),
        "\n", sep = "")
  }
}
