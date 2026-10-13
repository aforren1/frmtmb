# Lane optrap, consolidation fixes: the release's
# dev/rel070-optrap-diag-repoly.R (raw polynomials, x on [1, 2], 20
# groups of 10, 20 seeds per degree, against an exact profile maximum:
# GLS by Householder QR of the whitened design over log sigma and log
# tau, or lm() without the group term), with and without (1 | g), with
# the lane's parts taken out one at a time (OPTRAP_ABLATE as in
# dev/optrap-scale.R, plus "all" for predictor_units() empty).
#   Rscript dev/optrap-c3-repoly.R <lib> <label> [degrees] [re: 1|0]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
degs <- if (length(a) > 2) eval(parse(text = a[3])) else 5:8
re <- if (length(a) > 3) a[4] == "1" else TRUE
ab <- Sys.getenv("OPTRAP_ABLATE")
if (grepl("inner", ab)) {
  utils::assignInNamespace("predictor_inner_control",
                           function(units) list(maxit = 1000), "frmtmb")
}
if (grepl("theta", ab)) {
  utils::assignInNamespace("unit_theta_start",
                           function(frame, tpl, units = NULL) tpl, "frmtmb")
}
if (grepl("units", ab)) {
  utils::assignInNamespace("predictor_par_units",
                           function(...) NULL, "frmtmb")
}
if (grepl("all", ab)) {
  utils::assignInNamespace("predictor_units",
                           function(frame) numeric(0), "frmtmb")
}
cat("lib", dirname(find.package("frmtmb")), "frmtmb",
    as.character(packageVersion("frmtmb")), "label", a[2], "re", re,
    "ablated", ab, "\n")
ll_at <- function(y, X, Z, sig, tau) {
  V <- sig^2 * diag(length(y)) + tau^2 * tcrossprod(Z)
  R <- chol(V)
  yw <- backsolve(R, y, transpose = TRUE)
  Xw <- backsolve(R, X, transpose = TRUE)
  r <- qr.qty(qr(Xw, LAPACK = TRUE), yw)[-seq_len(ncol(X))]
  -0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(R))) + sum(r^2))
}
ll_lm <- function(y, X) {
  r <- qr.qty(qr(X, LAPACK = TRUE), y)[-seq_len(ncol(X))]
  n <- length(y)
  -0.5 * n * (log(2 * pi) + log(sum(r^2) / n) + 1)
}
X0 <- NULL
for (deg in degs) for (s in 1:20) {
  set.seed(1000 * deg + s)
  d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + "),
                         if (re) "+ (1 | g)"))
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
  if (re) {
    # lmer itself fails at degree 8 ("Downdated VtV is not positive
    # definite"); the profile then starts from the generating values
    m <- tryCatch(suppressWarnings(suppressMessages(lmer(fo, data = d,
      REML = FALSE))), error = function(e) NULL)
    st <- if (is.null(m)) log(c(0.2, 0.3)) else {
      vl <- as.data.frame(VarCorr(m))
      log(c(sigma(m), max(vl$sdcor[vl$grp == "g"], 1e-4)))
    }
    prof <- function(p) -ll_at(d$y, X, Z, exp(p[1]), exp(p[2]))
    o <- optim(st, prof, method = "Nelder-Mead",
               control = list(reltol = 1e-14, maxit = 5000))
    o <- optim(o$par, prof, method = "BFGS", control = list(reltol = 1e-15))
    best <- -o$value
    gl <- if (is.null(m)) NA else as.numeric(logLik(m)) - best
  } else {
    best <- ll_lm(d$y, X)
    gl <- NA
  }
  X0 <- rbind(X0, data.frame(
    deg = deg, seed = s, gap = as.numeric(logLik(f)) - best,
    gap_ref = gl, code = f$opt$convergence, nwarn = nw, nmsg = nm,
    unit = (f$par_units %||% 1)[1],
    steps = (f$opt$polish %||% c(steps = 0))[["steps"]], w1 = w1))
}
dir.create("dev/optrap-log/c3", showWarnings = FALSE)
write.table(X0, sprintf("dev/optrap-log/c3/repoly-%s.tsv", a[2]),
            sep = "\t", quote = FALSE, row.names = FALSE)
for (dg in degs) {
  x <- X0[X0$deg == dg, ]
  sh <- x$gap < -1e-4
  cat(sprintf(paste0("deg %d re %d | short %2d of 20 (min %9.3g) | silent ",
                     "short %2d | code!=0 %2d | warned %2d | polished %2d | ",
                     "lmer short %s | units %s\n"),
              dg, re, sum(sh), min(x$gap),
              sum(sh & x$code == 0 & x$nwarn == 0 & x$nmsg == 0),
              sum(x$code != 0), sum(x$nwarn > 0), sum(x$steps > 0),
              if (re) sum(x$gap_ref < -1e-4, na.rm = TRUE) else "-",
              paste(sort(unique(signif(x$unit, 3))), collapse = ",")))
}
w <- X0$w1[nzchar(X0$w1)]
if (length(w)) print(table(w))
