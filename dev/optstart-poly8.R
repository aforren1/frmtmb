# Lane optstart: test-se-flags.R's degree-8 raw polynomial with (1 | g),
# seed 3, from the lane's start and from 0.70.0's (given as `start`):
# objective against lme4, code, warnings, the SE check's verdicts, and
# the coefficient SEs against the same likelihood fitted with poly(x, 8)
# (well conditioned), its covariance mapped to the raw basis.
#   Rscript dev/optstart-poly8.R <lib> [seed]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
s <- if (length(a) > 1) as.integer(a[2]) else 3L
cat("lib", dirname(find.package("frmtmb")), "seed", s, "\n")
set.seed(s)
d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
  stats::rnorm(200, 0, 0.2)
for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + (1 | g)
l4 <- frm(y ~ poly(x, 8) + (1 | g), data = d)
Xr <- model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, d)
Xp <- model.matrix(~ poly(x, 8), d)
Tm <- qr.coef(qr(Xp), Xr)
Ti <- solve(Tm)
Vp <- as.matrix(vcov(l4))[seq_len(9), seq_len(9)]
se4 <- sqrt(diag(Ti %*% Vp %*% t(Ti)))
names(se4) <- c("Intercept", paste0("x", 1:8))
old <- list(beta = c(mean(d$y), rep(0, 8)), betad = log(stats::sd(d$y)))
for (nm in c("lane start", "0.70.0 start")) {
  w <- character()
  f <- withCallingHandlers(suppressMessages({
    m <- frm(fo, data = d, start = if (nm == "0.70.0 start") old)
    vcov(m)
    m
  }), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  se <- suppressWarnings(fixef(f)[, "Est.Error"])
  lost <- frmtmb:::sdr_of(f)$se_lost
  cat(sprintf("%s: code %d | logLik - poly fit %.4g | lost %s | warnings %s\n",
              nm, f$opt$convergence,
              as.numeric(logLik(f)) - as.numeric(logLik(l4)),
              paste(names(lost), "=", lost, collapse = ","),
              substr(paste(w, collapse = " | "), 1, 200)))
  print(signif(cbind(frmtmb = se, poly = se4[names(se)],
                     ratio = se / se4[names(se)]), 4))
}
