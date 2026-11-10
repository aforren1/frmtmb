# Reviewer scratch: what brms 2.23.0 does with each addition term on the
# ordinal families, at stancode (no compile). Seed 20261110.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")
for (f in c("cumulative", "sratio", "cratio", "acat", "hurdle_cumulative")) {
  cat(f, "ad:", paste(brms:::family_info(f, "ad"), collapse = ","), "\n")
}
set.seed(20261110)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n), nt = 4L, cc = 0L, s = 1,
                k = 1L, rr = 0.5, tt = 2, w = 1)
d$y <- as.integer(cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf)))
d$y2 <- d$y
d$y2[c(3, 9)] <- NA
spell <- c(trials = "trials(nt)", cens = "cens(cc)", trunc = "trunc(lb = 0)",
           se = "se(s)", mi = "mi()", vint = "vint(k)", vreal = "vreal(rr)",
           rate = "rate(tt)", weights = "weights(w)")
try1 <- function(fml, fam, ...) {
  r <- tryCatch({
    suppressWarnings(suppressMessages(stancode(fml, data = d, family = fam,
                                               ...)))
    "BUILDS"
  }, error = function(e) paste("ERR:", conditionMessage(e)))
  substr(gsub("\n", " ", r), 1, 160)
}
for (at in names(spell)) {
  yv <- if (at == "mi") "y2" else "y"
  lhs <- paste0(yv, " | ", spell[[at]])
  cat("\n==", at, "\n")
  cat(" plain cumulative:  ", try1(bf(as.formula(paste(lhs, "~ x"))),
                                    cumulative()), "\n")
  cat(" thres(3):          ", try1(bf(as.formula(paste0(yv, " | thres(3) + ",
                                    spell[[at]], " ~ x"))), cumulative()),
      "\n")
  cat(" disc ~ 0 + z:      ", try1(bf(as.formula(paste(lhs, "~ x")),
                                    disc ~ 0 + z), cumulative()), "\n")
  cat(" equidistant:       ", try1(bf(as.formula(paste(lhs, "~ x"))),
                                    cumulative(threshold = "equidistant")),
      "\n")
  cat(" sum_to_zero:       ", try1(bf(as.formula(paste(lhs, "~ x"))),
                                    cumulative(threshold = "sum_to_zero")),
      "\n")
  cat(" sratio cs(z):      ", try1(bf(as.formula(paste(lhs, "~ x + cs(z)"))),
                                    sratio()), "\n")
}
# the absent case in brms: the same terms beside cs() on a gaussian response
cat("\n== gaussian with cs(z) and se():",
    try1(bf(y | se(s) ~ x + cs(z)), gaussian()), "\n")
# mi() as a predictor beside an ordinal response: a different feature
# (mi_pred()), and brms builds it
d$xm <- d$x; d$xm[c(5, 7)] <- NA
cat("== mvbf ordinal y ~ mi(xm), xm | mi():",
    try1(bf(y ~ mi(xm) + cs(z), family = sratio()) +
           bf(xm | mi() ~ z, family = gaussian()) + set_rescor(FALSE),
         NULL), "\n")
