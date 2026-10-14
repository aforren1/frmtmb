# Lane optstart: the two backlog cases (dev/test-backlog.md, "Filed at
# the 0.70.0 release": the wide high-R^2 design and the group sd 1000
# times the residual sd), fitted from candidate starts given through
# `start =`, before any change to the start code.
#   Rscript dev/optstart-cases.R <lib or ""> <case> [seeds] [arms]
# (dev/optstart-cand.R with the arms to run as a fourth argument,
# comma-separated)
# case: wide (n 101000, 50 columns, y = 1 + 1e4 z1 + b_g + e, sd(b) 0.5)
#       re1000 (20 groups of 5, y = 1 + 0.5 x + b_g + e, sd(b) 1000)
# arms: default; ls_mad (mu at least squares, sigma at the MAD of its
# residuals); ls_within (mu at least squares, sigma at the MAD of the
# residuals with the group means taken out)
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
case <- a[2]
seeds <- if (length(a) > 2) eval(parse(text = a[3])) else 1:20
cat("lib", dirname(find.package("frmtmb")), "case", case, "\n")
mk <- function(s) {
  set.seed(s)
  if (case == "wide") {
    n <- 101000
    p <- 49
    Z <- matrix(rnorm(n * p), n, p)
    colnames(Z) <- paste0("z", seq_len(p))
    d <- data.frame(Z, g = factor(rep(seq_len(1000), length.out = n)))
    d$y <- 1 + 1e4 * d$z1 + rnorm(1000, 0, 0.5)[d$g] + rnorm(n)
    fo <- as.formula(paste("y ~", paste(colnames(Z), collapse = " + "),
                           "+ (1 | g)"))
  } else {
    S <- as.numeric(sub("re", "", case))
    d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
    d$y <- 1 + 0.5 * d$x + rnorm(20, 0, S)[d$g] + rnorm(100)
    fo <- y ~ x + (1 | g)
  }
  list(d = d, fo = fo)
}
for (s in seeds) {
  m <- mk(s)
  d <- m$d
  fo <- m$fo
  l4 <- as.numeric(logLik(lmer(fo, data = d, REML = FALSE,
    control = lmerControl(check.conv.singular = "ignore"))))
  X <- model.matrix(lme4::nobars(fo), d)
  q <- qr(X)
  r <- qr.resid(q, d$y)
  rw <- r - ave(r, d$g)
  arms <- list(default = NULL,
               ls_mad = list(beta = qr.coef(q, d$y), betad = log(mad(r))),
               ls_within = list(beta = qr.coef(q, d$y),
                                betad = log(mad(rw))))
  if (length(a) > 3) arms <- arms[strsplit(a[4], ",")[[1]]]
  for (nm in names(arms)) {
    w <- character()
    t0 <- proc.time()[[3]]
    f <- withCallingHandlers(suppressMessages(frm(fo, data = d,
                                                  start = arms[[nm]])),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      })
    cat(sprintf(paste0("seed %2d %-10s code %d | logLik - lmer %10.4g | ",
                       "iter %3d | warnings %d %s | secs %.1f\n"),
                s, nm, f$opt$convergence, as.numeric(logLik(f)) - l4,
                f$opt$iterations %||% NA, length(w),
                substr(paste(unique(substr(w, 1, 40)), collapse = " | "),
                       1, 80),
                proc.time()[[3]] - t0))
  }
}
