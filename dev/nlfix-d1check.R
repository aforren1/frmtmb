# D1 and D7 quick check on the lane's build: nl vs linear maxima under each
# ordinal family, seed 49 (dev/equiv-specs.R eq_dat), and every core autocor
# term with an nl body.
lib <- Sys.getenv("NLFIX_LIB",
                  "C:/Users/adf44/source/r/wt-nlfix-lib2;C:/Users/adf44/source/r/rellib-r16")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("dev/equiv-specs.R", envir = E)
d <- E$eq_dat("ordinal", seed = 49)
set.seed(49)
d$y0 <- ordered(ifelse(runif(nrow(d)) < 0.2, 0L, as.integer(d$y)),
                levels = 0:4)
ll <- function(expr) {
  r <- tryCatch(withCallingHandlers(expr,
    warning = function(w) {
      cat("  warning:", substr(conditionMessage(w), 1, 110), "\n")
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage")),
    error = function(e) e)
  if (inherits(r, "error")) {
    cat("  ERROR:", substr(conditionMessage(r), 1, 200), "\n")
    return(NA_real_)
  }
  as.numeric(logLik(r))
}
for (fam in c("cumulative", "sratio", "cratio", "acat", "hurdle_cumulative")) {
  Y <- if (fam == "hurdle_cumulative") "y0" else "y"
  F <- get(fam, envir = asNamespace("frmtmb"))
  ff <- function(s) stats::as.formula(gsub("Y", Y, s, fixed = TRUE))
  a <- ll(frm(bf(ff("Y ~ x")), family = F(), data = d))
  b <- ll(frm(bf(ff("Y ~ b1 * x"), b1 ~ 1, nl = TRUE), family = F(), data = d))
  c0 <- ll(frm(bf(ff("Y ~ 1")), family = F(), data = d))
  e <- ll(frm(bf(ff("Y ~ b0 + b1 * x"), b0 + b1 ~ 1, nl = TRUE),
              family = F(), data = d))
  f <- ll(frm(bf(ff("Y ~ b1 * exp(b2 * x)"), b1 ~ 1, b2 ~ 1 + z, nl = TRUE),
              family = F(), data = d,
              start = list(beta = c("b1_(Intercept)" = 0.5,
                                    "b2_(Intercept)" = 0.5))))
  g <- ll(frm(bf(ff("Y ~ b1 * x"), b1 ~ 1, disc ~ 0 + z, nl = TRUE),
              family = F(), data = d))
  h <- ll(frm(bf(ff("Y ~ x"), disc ~ 0 + z), family = F(), data = d))
  cat(sprintf(paste0("%-18s linear %.4f nl %.4f (diff %.2e) thresholds-only ",
                     "%.4f; b0+b1x %.4f; nonlinear %.4f; disc nl %.4f, ",
                     "linear %.4f\n"), fam, a, b, b - a, c0, e, f, g, h))
}
for (fam in c("cumulative", "sratio", "cratio", "acat")) {
  F <- get(fam, envir = asNamespace("frmtmb"))
  a <- ll(frm(bf(y | thres(3) ~ x), family = F(), data = d))
  b <- ll(frm(bf(y | thres(3) ~ b1 * x, b1 ~ 1, nl = TRUE), family = F(),
              data = d))
  cat(sprintf("%-18s thres(3): linear %.4f nl %.4f\n", fam, a, b))
}
cat("\ncs:\n")
ll(frm(bf(y ~ b1 * x + cs(z), b1 ~ 1, nl = TRUE), family = sratio(), data = d))
ll(frm(bf(y ~ b1 * x, b1 ~ 1 + cs(z), nl = TRUE), family = sratio(), data = d))
fit <- frm(bf(y ~ b1 * x, b1 ~ 1 + z, nl = TRUE), family = cumulative(),
           data = d)
print(fixef(fit))
print(head(fitted(fit), 2))
print(dim(predict(fit, newdata = d[1:3, ])))
