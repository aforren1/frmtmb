# Review of lane wcheck: common weighted models must fit, and fit as on
# rellib-r20. Prints one line per model: outcome, logLik, the fixed and
# dispersion coefficients' digest, nlminb code.
#   Rscript dev/reviews/wcheck-review-common.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(20261112)
ng <- 30; m <- 6; n <- ng * m
d <- data.frame(g = factor(rep(seq_len(ng), each = m)),
                h = factor(rep(seq_len(n / 2), each = 2)),
                x = rnorm(n), f = factor(sample(c("a", "b", "c"), n, TRUE)))
d$obs <- factor(seq_len(n))
ug <- rnorm(ng, 0, 0.5)[d$g]
eta <- 0.3 + 0.5 * d$x + ug
d$y <- eta + rnorm(n, 0, 0.7)
d$w <- runif(n, 0.3, 3)                     # survey weights, sum > n
d$wn <- d$w / mean(d$w)                     # normalized to mean 1
d$fw <- sample(1:4, n, TRUE)                # frequency weights
d$w0 <- ifelse(seq_len(n) %% 10 == 0, 0, d$w)
d$tr <- sample(5:15, n, TRUE)
d$k <- rbinom(n, d$tr, plogis(eta))
d$b <- rbinom(n, 1, plogis(eta))
d$cnt <- rpois(n, exp(eta + rnorm(n, 0, 0.4)))
d$yg <- rgamma(n, 3, 3 / exp(eta))
d$yb <- rbeta(n, plogis(eta) * 6, (1 - plogis(eta)) * 6)
d$yw <- rweibull(n, 2.5, exp(eta))
q8 <- quantile(d$y, 0.85)
d$cc <- ifelse(d$y > q8, "right", "none"); d$yc <- pmin(d$y, q8)
d$y2 <- d$y + rnorm(n, 0, 0.5)
# a few singleton groups: g2 equal to g except 6 rows in their own group
d$g2 <- as.character(d$g); d$g2[1:6] <- paste0("s", 1:6)
d$g2 <- factor(d$g2)

G <- gaussian()
cases <- list(
  list("gauss (1|g) survey w", bf(y | weights(w) ~ x + (1 | g)), G),
  list("gauss (1|g) normalized w", bf(y | weights(wn) ~ x + (1 | g)), G),
  list("gauss (1|g) freq w", bf(y | weights(fw) ~ x + (1 | g)), G),
  list("gauss (1+x|g) survey w", bf(y | weights(w) ~ x + (1 + x | g)), G),
  list("gauss (1|g/h) survey w", bf(y | weights(w) ~ x + (1 | g / h)), G),
  list("gauss (1|g2) singletons w", bf(y | weights(w) ~ x + (1 | g2)), G),
  list("gauss sigma ~ f (1|g) w", bf(y | weights(w) ~ x + (1 | g),
                                     sigma ~ f), G),
  list("gauss (1|g) zero weights", bf(y | weights(w0) ~ x + (1 | g)), G),
  list("gauss cens (1|g) w", bf(yc | weights(w) + cens(cc) ~ x + (1 | g)), G),
  list("gauss (1|obs) w normalized", bf(y | weights(wn) ~ x + (1 | g) +
                                          (1 | obs)), G),
  list("binomial (1|g) w", bf(k | trials(tr) + weights(w) ~ x + (1 | g)),
       binomial()),
  list("bernoulli (1|g) freq w", bf(b | weights(fw) ~ x + (1 | g)),
       bernoulli()),
  list("poisson (1|g) w", bf(cnt | weights(w) ~ x + (1 | g)), poisson()),
  list("poisson (1|obs) survey w", bf(cnt | weights(w) ~ x + (1 | obs)),
       poisson()),
  list("negbinomial (1|g) w", bf(cnt | weights(w) ~ x + (1 | g)),
       negbinomial()),
  list("Gamma (1|g) w", bf(yg | weights(w) ~ x + (1 | g)), Gamma("log")),
  list("Gamma shape ~ f (1|g) w", bf(yg | weights(w) ~ x + (1 | g),
                                     shape ~ f), Gamma("log")),
  list("lognormal (1|g) w", bf(yg | weights(w) ~ x + (1 | g)), lognormal()),
  list("Beta (1|g) w", bf(yb | weights(w) ~ x + (1 | g)), Beta()),
  list("student (1|g) w", bf(y | weights(w) ~ x + (1 | g)), student()),
  list("weibull (1|g) w", bf(yw | weights(w) ~ x + (1 | g)), weibull()),
  list("asym_laplace (1|g) w", bf(y | weights(w) ~ x + (1 | g)),
       asym_laplace()),
  list("mv gauss (1|p|g) w", mvbf(bf(y | weights(w) ~ x + (1 | p | g)),
                                  bf(y2 | weights(w) ~ x + (1 | p | g))), G),
  list("mv gauss (1|obs) + (1|g), w normalized",
       mvbf(bf(y | weights(wn) ~ x + (1 | g)),
            bf(cnt | weights(wn) ~ x + (1 | obs)) + poisson()), NULL))
for (cs in cases) {
  ws <- character()
  r <- tryCatch(withCallingHandlers(
    if (is.null(cs[[3]])) frm(cs[[2]], data = d) else
      frm(cs[[2]], family = cs[[3]], data = d),
    warning = function(w) {ws <<- c(ws, conditionMessage(w));
      invokeRestart("muffleWarning")},
    message = function(m) invokeRestart("muffleMessage")),
    error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("CASE %-38s ERROR %s\n", cs[[1]],
                substr(gsub("\\s+", " ", conditionMessage(r)), 1, 120)))
  } else {
    est <- unlist(r$estimates[c("beta", "betad", "theta")])
    cat(sprintf("CASE %-38s ok logLik %.10g | sum|est| %.10g | conv %s | warn %d\n",
                cs[[1]], as.numeric(logLik(r)), sum(abs(est)),
                r$opt$convergence, length(ws)))
  }
}
# lme4 comparison: weighted LMM, REML
if (requireNamespace("lme4", quietly = TRUE)) {
  lm4 <- lme4::lmer(y ~ x + (1 | g), data = d, weights = w)
  fr <- suppressWarnings(frm(y | weights(w) ~ x + (1 | g), data = d,
                             REML = TRUE))
  cat(sprintf("LME4 REML logLik lmer %.8f frm %.8f | fixef %s vs %s\n",
              as.numeric(logLik(lm4)), as.numeric(logLik(fr)),
              paste(format(lme4::fixef(lm4), digits = 8), collapse = " "),
              paste(format(fr$estimates$beta, digits = 8), collapse = " ")))
}
cat("DONE\n")
