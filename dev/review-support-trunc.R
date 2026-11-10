# Reviewer probe (lane support): trunc(lb <= 0) and cens() on positive
# families. FRMTMB_PORT_LIB picks the build. Within one process, a fit
# with a non-binding bound against the same fit without it.
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "\n")
set.seed(402)
n <- 300
e <- rexp(n)
u <- runif(n)
ys <- list(pos = 2 * e + 0.05, unit = qbeta(u, 2, 3))
fams <- list(
  exponential = list(exponential(), "pos"),
  gamma = list(Gamma(link = "log"), "pos"),
  weibull = list(weibull(), "pos"),
  lognormal = list(lognormal(), "pos"),
  frechet = list(frechet(), "pos"),
  inverse.gaussian = list(inverse.gaussian(), "pos"),
  shifted_lognormal = list(shifted_lognormal(), "pos"),
  beta = list(Beta(), "unit"),
  xbeta = list(xbeta(), "unit"),
  gev = list(gen_extreme_value(), "pos")
)
show <- function(nm, tag, f0, f1) {
  if (is.character(f0) || is.character(f1)) {
    cat(sprintf("%-18s %-26s ERR %s\n", nm, tag,
                substr(if (is.character(f1)) f1 else f0, 1, 110)))
    return(invisible())
  }
  l0 <- as.numeric(logLik(f0))
  l1 <- as.numeric(logLik(f1))
  b0 <- fixef(f0)[, 1]
  b1 <- fixef(f1)[, 1]
  cat(sprintf("%-18s %-26s dlogLik %.3g  max|dbeta| %.3g  (%.8g)\n", nm, tag,
              l1 - l0, max(abs(b1 - b0)), l0))
}
fitq <- function(fo, fam, d) {
  tryCatch(suppressWarnings(frm(fo, family = fam, data = d)),
           error = function(e) conditionMessage(e))
}
for (nm in names(fams)) {
  fam <- fams[[nm]][[1]]
  d <- data.frame(y = ys[[fams[[nm]][[2]]]])
  f0 <- fitq(y ~ 1, fam, d)
  for (lb in c(-1, 0)) {
    d$lb <- lb
    show(nm, paste0("trunc(lb=", lb, ")"), f0,
         fitq(y | trunc(lb = lb) ~ 1, fam, d))
  }
  # a bound that is not binding on both sides
  d$lb <- -1
  d$ub <- if (fams[[nm]][[2]] == "unit") 2 else 1e6
  show(nm, "trunc(lb=-1, ub=big)", f0,
       fitq(y | trunc(lb = lb, ub = ub) ~ 1, fam, d))
  # right censoring of 20% rows with and without a non-binding lb
  d$cc <- as.numeric(seq_len(n) %% 5 == 0)
  fc <- fitq(y | cens(cc) ~ 1, fam, d)
  show(nm, "cens + trunc(lb=-1)", fc,
       fitq(y | cens(cc) + trunc(lb = lb) ~ 1, fam, d))
  d$lb <- 0
  show(nm, "cens + trunc(lb=0)", fc,
       fitq(y | cens(cc) + trunc(lb = lb) ~ 1, fam, d))
}
