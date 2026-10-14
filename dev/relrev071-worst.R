# Release reviewer of 0.71.0: the merged worst case of lanes optstart
# and illse. A gaussian raw polynomial whose design condition number
# (unit-length columns) is inside optstart's least-squares start window
# (<= 1e5) and whose Hessian is ill-conditioned enough for illse's route
# (se_reparam_kappa 1e8), with (1 | g) on data with no group effect, so
# the sd runs to zero. Per fit: design kappa, whether the LS start
# applied, code, warnings, logLik against lmer (ML), max |SE / lmer - 1|
# over the coefficients, route, lost set.
#   Rscript dev/relrev071-worst.R <lib> <seeds> [sd_g]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "| R", R.home(), "\n")
seeds <- eval(parse(text = a[2]))
sdg <- if (length(a) > 2) as.numeric(a[3]) else 0
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e)
    structure(conditionMessage(e), class = "err")),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(x) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
wtag <- function(w) {
  if (!length(w)) return("-")
  paste(substr(gsub("\n", " ", w), 1, 60), collapse = " || ")
}
# x on [1, 1 + width]: narrower is worse conditioned
cases <- list(c(deg = 3, wid = 1), c(deg = 4, wid = 1), c(deg = 2, wid = 0.1),
              c(deg = 3, wid = 0.3), c(deg = 2, wid = 0.03))
for (cs in cases) for (seed in seeds) {
  deg <- cs[["deg"]]; wid <- cs[["wid"]]
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::runif(n, 1, 1 + wid),
                  g = factor(rep(1:20, length.out = n)))
  d$y <- sin(3 * d$x / wid) + stats::rnorm(20, 0, sdg)[d$g] +
    stats::rnorm(n, 0, 0.3)
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  X <- cbind(1, as.matrix(d[xn]))
  kap <- kappa(sweep(X, 2, sqrt(colSums(X^2)), "/"), exact = TRUE)
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                                "+ (1 | g)"))
  r <- cap(frm(fo, data = d))
  f <- r$v
  if (inherits(f, "err")) {
    cat(sprintf("deg %d wid %g s%2d | ERROR %s\n", deg, wid, seed, f)); next
  }
  m <- suppressWarnings(suppressMessages(lme4::lmer(fo, data = d,
                                                    REML = FALSE)))
  ls <- exists("start_ls_plan", ns) && length(ns$start_ls_plan(f$frame)) > 0
  cn <- c("(Intercept)", xn)
  V <- suppressWarnings(vcov(f))
  se <- sqrt(diag(as.matrix(V)))[seq_along(cn)]
  sr <- sqrt(diag(as.matrix(vcov(m))))
  l <- ns$sdr_of(f)$se_lost
  cat(sprintf(paste0("deg %d wid %-4g s%2d | kappa %.1e | ls %s | code %d |",
                     " sd %.1e | dLL %+.2e | maxrelSE %.2e | route %s |",
                     " lost %s | warn %s\n"),
              deg, wid, seed, kap, ls, f$opt$convergence,
              exp(f$opt$par[names(f$opt$par) == "theta"][1]),
              as.numeric(logLik(f)) - as.numeric(logLik(m)),
              max(abs(se / sr - 1)),
              !is.null(ns$sdr_of(f)$cov_fixed_root),
              if (length(l)) paste(names(l), l, sep = ":", collapse = ",")
              else "-", wtag(r$w)))
}
