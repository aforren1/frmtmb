# Lane optstart: test-v11.R's scalar-on-function regression (seed 602),
# s(Tmat, by = Lmat), from the lane's start and from 0.70.0's (given as
# `start`): objective against mgcv's ML, code, and the smoothing sd.
#   Rscript dev/optstart-sof.R <lib>
.libPaths(c(commandArgs(TRUE)[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
set.seed(602)
n <- 250; K <- 30
tt <- seq(0, 1, length.out = K)
Xc <- matrix(rnorm(n * K), n, K)
for (i in seq_len(n)) Xc[i, ] <- cumsum(rnorm(K, 0, 0.3))
y <- 1 + drop(Xc %*% sin(2 * pi * tt)) / K + rnorm(n, 0, 0.3)
dd <- data.frame(y = y)
dd$Tmat <- matrix(tt, n, K, byrow = TRUE)
dd$Lmat <- Xc / K
ref <- mgcv::gam(y ~ s(Tmat, by = Lmat), data = dd, method = "ML")
st <- par_template(bf(y ~ s(Tmat, by = Lmat)), data = dd, family = gaussian())
print(unlist(st[c("beta", "betad", "theta")]))
old <- list(beta = c(mean(dd$y), rep(0, length(st$beta) - 1)),
            betad = log(stats::sd(dd$y)))
for (nm in c("lane", "old")) {
  f <- suppressWarnings(frm(bf(y ~ s(Tmat, by = Lmat)), family = gaussian(),
                            data = dd, start = if (nm == "old") old,
                            control = frmtmb_control(verbose = TRUE)))
  cat(sprintf("%s: code %d | logLik - mgcv %.5g | theta %s | sigma %.4g\n", nm,
              f$opt$convergence, as.numeric(logLik(f)) + ref$gcv.ubre,
              paste(signif(f$opt$par[names(f$opt$par) == "theta"], 4),
                    collapse = " "),
              exp(f$opt$par[names(f$opt$par) == "betad"])))
}
