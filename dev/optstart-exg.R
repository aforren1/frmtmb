# Lane optstart: the battery's exgaussian with a group sd 30 times sigma
# (dev/optstart-battery.R exgaussian rebig), where fits stop at "NA/NaN
# gradient". The start, the objective and gradient there, and where the
# optimizer was when the gradient went NaN.
#   Rscript dev/optstart-exg.R <lib or ""> [seeds]
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 4
cat("lib", dirname(find.package("frmtmb")), "\n")
for (s in seeds) {
  set.seed(s)
  n <- 400
  g <- factor(rep(1:20, length.out = n))
  x <- rnorm(n)
  bg <- rnorm(20, 0, 30)
  mu <- 0.3 + 0.5 * x + bg[g]
  y <- frmtmb:::family_registry$exgaussian()$sim(
    list(mu = mu, sigma = rep(1, n), beta = rep(1.5, n)), list(), n)
  d <- data.frame(y = y, x = x, g = g)
  st <- par_template(bf(y ~ x + (1 | g)), data = d, family = exgaussian())
  cat("seed", s, "start:", format(unlist(st[c("beta", "betad", "theta")]),
                                 digits = 4), "\n")
  r <- tryCatch(frm(y ~ x + (1 | g), data = d, family = exgaussian(),
                    control = frmtmb_control(verbose = TRUE)),
                error = function(e) conditionMessage(e))
  if (is.character(r)) cat("  ERROR", substr(r, 1, 200), "\n") else
    cat("  code", r$opt$convergence, "logLik", as.numeric(logLik(r)), "\n")
}
