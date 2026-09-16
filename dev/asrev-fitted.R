## REVIEW 2.5e, risk: fitted() is now ONE call to predict(). Did the
## answer change for any existing call?
## Usage: Rscript dev/asrev-fitted.R <lib> <outfile>
args <- commandArgs(trailingOnly = TRUE)
LIB <- args[1]
OUT <- args[2]
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", dirname(system.file(package = "frmtmb")), "\n")

set.seed(11)
n <- 90
dd <- data.frame(x = rnorm(n), g = factor(rep(1:9, each = 10)))
dd$lin <- 0.3 + 0.5 * dd$x + rnorm(9, 0, 0.6)[dd$g]
dd$yg <- dd$lin + rnorm(n, 0, 0.5)
dd$yp <- rpois(n, exp(dd$lin))
dd$yb <- rbinom(n, 1, plogis(dd$lin))
dd$trials <- rep(10L, n)
dd$ybt <- rbinom(n, 10, plogis(dd$lin))
dd$ylog <- exp(dd$lin + rnorm(n, 0, 0.3))
dd$yzi <- ifelse(runif(n) < 0.3, 0L, rpois(n, exp(dd$lin)))
dd$yord <- factor(cut(dd$lin + rnorm(n, 0, 0.5), 3), labels = 1:3,
                  ordered = TRUE)
dd$ycat <- factor(sample(c("a", "b", "c"), n, TRUE))
dd$yg2 <- dd$lin * 0.7 + rnorm(n, 0, 0.5)
ddna <- dd
ddna$x[c(3, 17, 40)] <- NA

cases <- list(
  gaussian_re = function() frm(bf(yg ~ x + (1 | g)) + gaussian(), data = dd),
  poisson_re  = function() frm(bf(yp ~ x + (1 | g)) + poisson(), data = dd),
  bernoulli   = function() frm(bf(yb ~ x + (1 | g)) + bernoulli(), data = dd),
  binom_tr    = function() frm(bf(ybt | trials(trials) ~ x) + binomial(),
                               data = dd),
  lognormal   = function() frm(bf(ylog ~ x) + lognormal(), data = dd),
  zip         = function() frm(bf(yzi ~ x, zi ~ 1) + zero_inflated_poisson(),
                               data = dd),
  ordinal     = function() frm(bf(yord ~ x) + cumulative(), data = dd),
  categorical = function() frm(bf(ycat ~ x) + categorical(), data = dd),
  sigma_dpar  = function() frm(bf(yg ~ x, sigma ~ x) + gaussian(), data = dd),
  na_exclude  = function() frm(bf(yg ~ x + (1 | g)) + gaussian(), data = ddna,
                               na.action = stats::na.exclude),
  na_omit     = function() frm(bf(yg ~ x + (1 | g)) + gaussian(), data = ddna,
                               na.action = stats::na.omit),
  mv          = function() frm(bf(yg ~ x) + bf(yg2 ~ x) + gaussian(),
                               data = dd)
)

res <- list()
for (nm in names(cases)) {
  fit <- tryCatch(cases[[nm]](), error = function(e) e)
  if (inherits(fit, "error")) {
    res[[nm]] <- list(err = conditionMessage(fit)); next
  }
  f <- tryCatch(stats::fitted(fit), error = function(e) conditionMessage(e))
  p <- tryCatch(stats::predict(fit, type = "response"),
                error = function(e) conditionMessage(e))
  res[[nm]] <- list(fitted = f, predict_resp = p,
                    coef = tryCatch(unlist(fixef(fit)),
                                    error = function(e) NA))
  cat(sprintf("%-12s fitted: %s  %s\n", nm,
              paste(class(f), collapse = "/"),
              if (is.numeric(f)) sprintf("len=%d first=%.10g", length(f), f[1])
              else if (is.character(f)) substr(f, 1, 70) else ""))
}
saveRDS(res, OUT)
cat("wrote", OUT, "\n")

cat("\n---- identity fitted() == predict(type='response') ----\n")
for (nm in names(res)) {
  r <- res[[nm]]
  if (!is.null(r$err)) { cat(sprintf("%-12s FIT ERROR\n", nm)); next }
  if (is.character(r$fitted) || is.character(r$predict_resp)) {
    cat(sprintf("%-12s fitted=%s predict=%s\n", nm,
                if (is.character(r$fitted)) "ERR" else "ok",
                if (is.character(r$predict_resp)) "ERR" else "ok"))
    next
  }
  cat(sprintf("%-12s identical=%-5s all.equal=%s\n", nm,
              identical(r$fitted, r$predict_resp),
              isTRUE(all.equal(r$fitted, r$predict_resp))))
}
