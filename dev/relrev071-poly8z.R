# Release reviewer of 0.71.0: lane ci701's "illcond" composition case
# (dev/ci701-rev-compose.R: degree-8 raw polynomial with (1 + z | g),
# the correlation running toward -1) on a build. The same model in
# poly()'s orthogonal basis is the reference: logLik, and the
# population fitted values' SEs on a grid (re_formula = NA).
#   Rscript dev/relrev071-poly8z.R <lib> [seeds]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "|", R.home(), "\n")
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 3
conds <- function(expr) {
  out <- character(0)
  v <- withCallingHandlers(expr, warning = function(w) {
    out <<- c(out, paste("W:", conditionMessage(w)))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    out <<- c(out, paste("M:", conditionMessage(m)))
    invokeRestart("muffleMessage")
  })
  list(value = v, conds = out)
}
for (s in seeds) {
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)),
                  z = stats::rnorm(200))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + z + (1 + z | g)
  r <- conds(frm(fo, data = d))
  f <- r$value
  ro <- conds(frm(y ~ poly(x, 8) + z + (1 + z | g), data = d))
  fo2 <- ro$value
  nd <- data.frame(x = seq(1.05, 1.95, length.out = 10), z = 0, g = d$g[1])
  for (k in 1:8) nd[[paste0("x", k)]] <- nd$x^k
  p1 <- tryCatch(frm_linpred(f, newdata = nd, re_formula = NA,
                             se.fit = TRUE), error = function(e) NULL)
  p2 <- frm_linpred(fo2, newdata = nd, re_formula = NA, se.fit = TRUE)
  g <- max(abs(drop(f$obj$gr(f$opt$par)) * (f$par_units %||% 1)))
  l <- ns$sdr_of(f)$se_lost
  se <- sqrt(diag(as.matrix(suppressWarnings(vcov(f)))))
  th <- f$estimates$theta
  cat(sprintf(paste0("seed %d raw: code %d t %.4g max|g u| %.3g | logLik raw",
                     " - orth %+.2e | orth code %d | lost %s | coef SE",
                     " finite %d of 10 | route %s\n"),
              s, f$opt$convergence, th[length(th)], g,
              as.numeric(logLik(f)) - as.numeric(logLik(fo2)),
              fo2$opt$convergence,
              if (length(l)) paste(names(l), l, sep = ":", collapse = ",")
              else "none", sum(is.finite(se[1:10])),
              !is.null(ns$sdr_of(f)$cov_fixed_root)))
  if (!is.null(p1)) {
    cat(sprintf("   fitted (re_formula NA): max |eta raw - orth| %.2e, se raw / orth in [%.4g, %.4g]\n",
                max(abs(p1$fit - p2$fit)), min(p1$se.fit / p2$se.fit),
                max(p1$se.fit / p2$se.fit)))
  }
  for (x in r$conds) cat("   raw", substr(gsub("\\s+", " ", x), 1, 230), "\n")
  for (x in ro$conds) cat("   orth", substr(gsub("\\s+", " ", x), 1, 160), "\n")
}
