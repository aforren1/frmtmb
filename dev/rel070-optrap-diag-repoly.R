# Release 0.70.0: the probe found frmtmb short of the exact maximum on
# a raw polynomial with (1 | g) (degree 6, seed 6: 0.204 short at code
# 0, silent, where 0.69.1 reaches it). One seed is not a rate, so this
# repeats the probe's construction (x on [1, 2], 20 groups of 10) over
# 20 seeds per degree and counts fits short of an exact profile maximum
# (GLS by Householder QR of the whitened design, maximized over log
# sigma and log tau), and which of them said anything.
#   Rscript dev/rel070-optrap-diag-repoly.R <lib> <label> [degrees]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
degs <- if (length(a) > 2) eval(parse(text = a[3])) else 4:7
cat("lib", dirname(find.package("frmtmb")), "frmtmb",
    as.character(packageVersion("frmtmb")), "label", a[2], "\n")
ll_at <- function(y, X, Z, sig, tau) {
  V <- sig^2 * diag(length(y)) + tau^2 * tcrossprod(Z)
  R <- chol(V)
  yw <- backsolve(R, y, transpose = TRUE)
  Xw <- backsolve(R, X, transpose = TRUE)
  # the residual is the trailing part of Q'y, which a LAPACK QR gives
  r <- qr.qty(qr(Xw, LAPACK = TRUE), yw)[-seq_len(ncol(X))]
  -0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(R))) + sum(r^2))
}
X0 <- NULL
for (deg in degs) for (s in 1:20) {
  set.seed(1000 * deg + s)
  d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + "),
                         "+ (1 | g)"))
  X <- cbind(1, as.matrix(d[paste0("x", 1:deg)]))
  Z <- model.matrix(~ 0 + g, d)
  nw <- 0L; nm <- 0L; w1 <- ""
  f <- withCallingHandlers(frm(fo, data = d),
    warning = function(x) {
      nw <<- nw + 1L
      if (!nzchar(w1)) w1 <<- substr(conditionMessage(x), 1, 60)
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      nm <<- nm + 1L; invokeRestart("muffleMessage")
    })
  m <- suppressWarnings(suppressMessages(lmer(fo, data = d, REML = FALSE)))
  vl <- as.data.frame(VarCorr(m))
  st <- log(c(sigma(m), max(vl$sdcor[vl$grp == "g"], 1e-4)))
  prof <- function(p) -ll_at(d$y, X, Z, exp(p[1]), exp(p[2]))
  o <- optim(st, prof, method = "Nelder-Mead",
             control = list(reltol = 1e-14, maxit = 5000))
  o <- optim(o$par, prof, method = "BFGS", control = list(reltol = 1e-15))
  best <- -o$value
  X0 <- rbind(X0, data.frame(
    deg = deg, seed = s, gap = as.numeric(logLik(f)) - best,
    gap_lmer = as.numeric(logLik(m)) - best, code = f$opt$convergence,
    nwarn = nw, nmsg = nm, unit = (f$par_units %||% 1)[1], w1 = w1))
}
write.table(X0, sprintf("dev/rel070-log/studies/optrap/diag-repoly-%s.tsv",
                        a[2]), sep = "\t", quote = FALSE, row.names = FALSE)
for (dg in degs) {
  x <- X0[X0$deg == dg, ]
  sh <- x$gap < -1e-4
  cat(sprintf(paste0("deg %d | short %2d of 20 (min %9.3g) | silent short ",
                     "%2d | code!=0 %2d | warned %2d | lmer short %2d ",
                     "(min %9.3g) | units %s\n"),
              dg, sum(sh), min(x$gap),
              sum(sh & x$code == 0 & x$nwarn == 0 & x$nmsg == 0),
              sum(x$code != 0), sum(x$nwarn > 0), sum(x$gap_lmer < -1e-4),
              min(x$gap_lmer),
              paste(sort(unique(signif(x$unit, 3))), collapse = ",")))
}
w <- X0$w1[nzchar(X0$w1)]
if (length(w)) print(table(w))
