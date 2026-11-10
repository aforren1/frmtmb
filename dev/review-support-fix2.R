# Reviewer probe, short check of the fix round (lane support):
# shifted_lognormal with a bound between 0 and ndt, von Mises with
# bounds past both ends, and the option's restore after an error.
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
fitq <- function(fo, fam, d) {
  tryCatch(suppressWarnings(frm(fo, family = fam, data = d)),
           error = function(e) conditionMessage(e))
}
cmp <- function(tag, f0, f1) {
  if (is.character(f0) || is.character(f1)) {
    cat(sprintf("%-40s ERR %s\n", tag,
                substr(if (is.character(f1)) f1 else f0, 1, 120)))
  } else {
    cat(sprintf("%-40s dlogLik %.3g  max|dpar| %.3g\n", tag,
                as.numeric(logLik(f1)) - as.numeric(logLik(f0)),
                max(abs(f1$opt$par - f0$opt$par))))
  }
}
set.seed(402)
n <- 300
d <- data.frame(y = 2 * rexp(n) + 0.05)
f0 <- fitq(y ~ 1, shifted_lognormal(), d)
cat("shifted_lognormal ndt_hat",
    if (!is.character(f0)) format(frmtmb:::eval_dpars(f0)[[1]]$ndt[1]), "\n")
for (lb in c(0.001, 0.01, 0.03)) {
  d$lb <- lb
  cmp(paste0("shifted_lognormal trunc(lb=", lb, ")"), f0,
      fitq(y | trunc(lb = lb) ~ 1, shifted_lognormal(), d))
}
set.seed(406)
v <- data.frame(y = atan2(sin(rnorm(n, 0.5, 0.8)), cos(rnorm(n, 0.5, 0.8))))
v0 <- fitq(y ~ 1, von_mises(), v)
v$lb <- -4
v$ub <- 4
cmp("von_mises trunc(lb=-4)", v0, fitq(y | trunc(lb = lb) ~ 1, von_mises(), v))
cmp("von_mises trunc(ub=4)", v0, fitq(y | trunc(ub = ub) ~ 1, von_mises(), v))
cmp("von_mises trunc(lb=-4, ub=4)", v0,
    fitq(y | trunc(lb = lb, ub = ub) ~ 1, von_mises(), v))
v$cc <- as.numeric(seq_len(n) %% 5 == 0)
vc <- fitq(y | cens(cc) ~ 1, von_mises(), v)
cmp("von_mises cens + trunc(lb=-4, ub=4)", vc,
    fitq(y | cens(cc) + trunc(lb = lb, ub = ub) ~ 1, von_mises(), v))
# the option after an error inside an exact reader
bad <- frmtmb::frmtmb_family("bad", dpars = "mu", links = list(mu = "identity"),
                             lpdf = function(y, dpars, aterms) stop("boom"))
r <- try(frmtmb:::row_lpdf(bad, 1, 1, list(mu = 0), list(), NULL,
                           exact = TRUE), silent = TRUE)
cat("after error: option", format(getOption("frmtmb.exact_support")), "\n")
