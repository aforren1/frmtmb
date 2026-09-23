# The two remaining punch minors that need a measurement: how tight the
# quadrature rescue really is, and whether mixture() now carries a
# component's stationary declaration.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")
op <- options(digits = 12)
gap <- function(f, ref) as.numeric(logLik(f)) - ref

cat("\n=== quadrature rescue, acceptance seed 1 plus (1 | g)\n")
dd <- make_data(1, FALSE)
set.seed(99)
dd$g <- factor(rep(1:20, length.out = nrow(dd)))
m3 <- skew(dd$y)
bad <- list(betad = c(log(stats::sd(dd$y)), 2 * sign(m3) + 0.5 * m3))
fo <- bf(y ~ xs + (1 | g), sigma ~ 1, alpha ~ 1)
for (q in c(FALSE, TRUE)) {
  d0 <- suppressWarnings(frm(fo, family = skew_normal(), data = dd,
                             quadrature = q))
  d1 <- suppressWarnings(frm(fo, family = skew_normal(), data = dd,
                             quadrature = q, start = bad))
  e <- d1$opt[["stationary_escape"]]
  cat(sprintf("  quadrature=%-5s default %16.9f  rescued %16.9f  d %11.3e  esc %s\n",
              q, as.numeric(logLik(d0)), as.numeric(logLik(d1)),
              as.numeric(logLik(d1)) - as.numeric(logLik(d0)),
              if (is.null(e)) "-" else format(e[["gain"]], digits = 5)))
  # is the difference an escape shortfall, or does the DEFAULT fit also
  # sit that far from the best of several starts?
  alt <- vapply(c(1, 3, 6), function(a) {
    as.numeric(logLik(suppressWarnings(
      frm(fo, family = skew_normal(), data = dd, quadrature = q,
          start = list(betad = c(log(stats::sd(dd$y)), a))))))
  }, 0)
  cat(sprintf("           best of starts 1,3,6 %16.9f  spread %11.3e\n",
              max(alt), max(alt) - min(alt)))
}

cat("\n=== mixture() carries a component's stationary declaration\n")
fam <- mixture(skew_normal(), skew_normal())
st <- fam[["post"]][["stationary"]]
cat("  declared dpars:", if (is.null(st)) "NONE" else
  paste(names(st), collapse = " "), "\n")
set.seed(4)
n <- 400
z <- rbinom(n, 1, 0.5)
ym <- ifelse(z == 1, 6, 0) +
  (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
dm <- data.frame(y = ym)
f <- try(suppressWarnings(frm(bf(y ~ 1), family = fam, data = dm)),
         silent = TRUE)
if (inherits(f, "try-error")) {
  cat("  mixture fit ERROR:", sub("\n.*", "", attr(f, "condition")$message),
      "\n")
} else {
  a <- f$estimates$betad[grep("^alpha", names(f$estimates$betad))]
  cat(sprintf("  mixture fit ll %14.6f  alphas %s  esc %s\n",
              as.numeric(logLik(f)), paste(format(a, digits = 4),
                                           collapse = ", "),
              if (is.null(f$opt[["stationary_escape"]])) "-" else "fired"))
}
options(op)
