# Lane famend, 3003: exgaussian() with a group sd 30 times sigma
# (dev/optstart-exg.R's design), from the default start and from the
# true values, on one build.
#   Rscript dev/famend-probe-exg.R <lib or "base"> [seeds]
a <- commandArgs(trailingOnly = TRUE)
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] != "base") a[1], base, ulib))
suppressPackageStartupMessages(library(frmtmb))
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:40
cat("lib", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
one <- function(expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }), error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(r, "err")) return(c(cls = "ERROR", ll = NA, code = NA))
  c(cls = if (length(w)) "warn" else "ok",
    ll = as.numeric(logLik(r)), code = r$opt$convergence)
}
tab <- NULL
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
  a1 <- one(frm(y ~ x + (1 | g), data = d, family = exgaussian()))
  a2 <- one(frm(y ~ x + (1 | g), data = d, family = exgaussian(),
                start = list(beta = c(0.3, 0.5), betad = log(c(1, 1.5)),
                             theta = log(30))))
  tab <- rbind(tab, data.frame(seed = s, default = a1[["cls"]],
                               code = a1[["code"]], ll = a1[["ll"]],
                               truth = a2[["cls"]], ll_truth = a2[["ll"]]))
}
print(tab, row.names = FALSE)
cat(sprintf(paste0("SUMMARY seeds %d: default errors %d, code != 0 %d, ",
                   "warned %d; from the truth errors %d; max |ll - ",
                   "ll_truth| %.3g\n"), nrow(tab),
            sum(tab$default == "ERROR"), sum(tab$code != "0", na.rm = TRUE),
            sum(tab$default == "warn"), sum(tab$truth == "ERROR"),
            max(abs(as.numeric(tab$ll) - as.numeric(tab$ll_truth)),
                na.rm = TRUE)))
