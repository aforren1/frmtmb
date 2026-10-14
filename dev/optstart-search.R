# Lane optstart: a design small enough for a test where 0.70.0's start
# stops short of lme4 in silence: y = 1 + k z1 + ... + b_g + e, sd(b)
# 0.5, residual sd 1, over n, the number of columns p, k and seeds.
# Per setting: fits short of lmer by more than 1e-4 log-likelihood
# units, with their codes and warnings.
#   Rscript dev/optstart-search.R <lib or ""> [seeds]
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:20
cat("lib", dirname(find.package("frmtmb")), "\n")
for (n in c(500, 2000, 10000)) for (p in c(2, 10)) for (k in c(1e3, 1e5)) {
  short <- character()
  for (s in seeds) {
    set.seed(s)
    Z <- matrix(rnorm(n * p), n, p)
    colnames(Z) <- paste0("z", seq_len(p))
    G <- max(10, n / 100)
    d <- data.frame(Z, g = factor(rep(seq_len(G), length.out = n)))
    d$y <- 1 + k * d$z1 + rnorm(G, 0, 0.5)[d$g] + rnorm(n)
    fo <- as.formula(paste("y ~", paste(colnames(Z), collapse = " + "),
                           "+ (1 | g)"))
    w <- character()
    f <- withCallingHandlers(suppressMessages(frm(fo, data = d)),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      })
    l4 <- as.numeric(logLik(lmer(fo, data = d, REML = FALSE,
      control = lmerControl(check.conv.singular = "ignore"))))
    dl <- as.numeric(logLik(f)) - l4
    if (dl < -1e-4) {
      short <- c(short, sprintf("seed %d %.4g code %d warn %d", s, dl,
                                f$opt$convergence, length(w)))
    }
  }
  cat(sprintf("n %5d p %2d k %g: short %d of %d %s\n", n, p, k,
              length(short), length(seeds), paste(short, collapse = "; ")))
}
