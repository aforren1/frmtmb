# Log-density agreement of frmtmb.eam's wiener against brms 2.23.0's
# Stan density (its own generated functions block, exposed through
# rstan), RWiener, WienR, rtdists, and a 1200-bit Rmpfr series reference.
# Grid: decision time x boundary x bias x drift x boundary hit.
# Output: dev/famuni-out/wiener-density.txt and wiener-density.rds
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all(
  "C:/Users/adf44/source/r/frmtmb-wt-famuni/extensions/frmtmb.eam",
  quiet = TRUE))
suppressMessages({library(brms); library(rstan)})
out_dir <- "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-out"
cat("brms", format(packageVersion("brms")), "rstan",
    format(packageVersion("rstan")), "stan", rstan::stan_version(), "\n")
cat("RWiener", format(packageVersion("RWiener")), "WienR",
    format(packageVersion("WienR")), "rtdists",
    format(packageVersion("rtdists")), "Rmpfr",
    format(packageVersion("Rmpfr")), "\n")

## brms's own functions block, cut out of its generated program
set.seed(1)
d0 <- data.frame(rt = c(0.5, 0.7), d = c(0L, 1L))
sc <- strsplit(as.character(stancode(bf(rt | dec(d) ~ 1),
                                     family = wiener(), data = d0)),
               "\n")[[1]]
i0 <- grep("^functions [{]", sc); i1 <- grep("^data [{]", sc) - 1L
fun_code <- paste(sc[i0:i1], collapse = "\n")
cat("\n-- brms functions block\n", fun_code, "\n")
t_comp <- system.time(
  expose_stan_functions(stanc(model_code = fun_code), env = globalenv()))
cat("compiled in", t_comp[["elapsed"]], "s\n")

## grid
tt <- c(1e-4, 1e-3, 0.01, 0.05, 0.1, 0.3, 1, 3, 10)
aa <- c(0.5, 1, 2, 4)
ww <- c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
vv <- c(-20, -10, -5, -1, 0, 1, 5, 10, 20)
g <- expand.grid(t = tt, a = aa, w = ww, v = vv, up = c(0L, 1L))
tau <- 0.1
g$y <- g$t + tau
# the decision time each implementation sees, formed the way Stan forms
# it, so the comparison is not of two roundings of t
g$td <- g$y - tau
cat("grid points", nrow(g), "\n")

## 1200-bit reference. Lower-boundary density; the upper boundary is the
## lower one at v -> -v, w -> 1 - w. The small-time series for u < 1, the
## large-time series otherwise; both are exact series, truncated at 60
## terms where the next term is below 1e-1000 of the sum.
prec <- 1200
lref1 <- function(t, v, a, w) {
  t <- Rmpfr::mpfr(t, prec); v <- Rmpfr::mpfr(v, prec)
  a <- Rmpfr::mpfr(a, prec); w <- Rmpfr::mpfr(w, prec)
  u <- t / (a * a)
  pre <- -v * a * w - v * v * t / 2 - 2 * log(a)
  if (as.numeric(u) < 1) {
    k <- Rmpfr::mpfr(-60:60, prec)
    s <- sum((w + 2 * k) * exp(-(w + 2 * k)^2 / (2 * u)))
    as.numeric(pre + log(s) - 0.5 * log(2 * Rmpfr::Const("pi", prec)) -
                 1.5 * log(u))
  } else {
    k <- Rmpfr::mpfr(1:60, prec)
    p <- Rmpfr::Const("pi", prec)
    s <- sum(k * exp(-k * k * p * p * u / 2) * sin(k * p * w))
    as.numeric(pre + log(p) + log(s))
  }
}
t_ref <- system.time(g$ref <- vapply(seq_len(nrow(g)), function(i) {
  up <- g$up[i]
  lref1(g$td[i], if (up == 1) -g$v[i] else g$v[i], g$a[i],
        if (up == 1) 1 - g$w[i] else g$w[i])
}, 0))
cat("reference in", t_ref[["elapsed"]], "s\n")

## eam (the plain family's density is ddm_lpdf_both(y - ndt, ...))
g$eam <- ddm_lpdf_both(g$td, g$v, g$a, g$w, g$up)
## brms Stan
stan_err <- character(0)
g$stan <- vapply(seq_len(nrow(g)), function(i) {
  tryCatch(wiener_diffusion_lpdf(g$y[i], g$up[i], g$a[i], tau, g$w[i],
                                 g$v[i]),
           error = function(e) {
             stan_err <<- c(stan_err, conditionMessage(e)); NA_real_
           })
}, 0)
## brms's R-side dwiener, which log_lik() uses (RWiener backend)
g$brms_r <- brms::dwiener(g$y, alpha = g$a, tau = tau, beta = g$w,
                          delta = g$v, resp = g$up, log = TRUE)
g$rwiener <- vapply(seq_len(nrow(g)), function(i) {
  RWiener::dwiener(g$y[i], g$a[i], tau, g$w[i], g$v[i],
                   resp = if (g$up[i] == 1) "upper" else "lower",
                   give_log = TRUE)
}, 0)
g$wienr <- WienR::dWDM(g$td, ifelse(g$up == 1, "upper", "lower"), g$a, g$v,
                       g$w, precision = 12)$logvalue
g$rtdists <- log(rtdists::ddiffusion(g$y, ifelse(g$up == 1, "upper",
                                                 "lower"),
                                     a = g$a, v = g$v, t0 = tau,
                                     z = g$w * g$a, precision = 6))
saveRDS(g, file.path(out_dir, "wiener-density.rds"))

cat("\nStan errors:", length(stan_err), "\n")
if (length(stan_err)) print(head(unique(stan_err), 5))

## summaries
g$u <- g$td / g$a^2
fin <- is.finite(g$ref)
cat("reference finite at", sum(fin), "of", nrow(g), "\n")
cat("reference range", format(range(g$ref[fin])), "\n")
summ <- function(col, sel = rep(TRUE, nrow(g)), lab = "all") {
  x <- g[[col]][sel & fin]; r <- g$ref[sel & fin]
  e <- abs(x - r)
  bad <- !is.finite(x)
  sprintf("%-8s %-22s n=%4d nonfinite=%4d maxabs=%.3e median=%.3e maxrel=%.3e",
          col, lab, length(x), sum(bad), max(e[!bad]), median(e[!bad]),
          max((e / pmax(1, abs(r)))[!bad]))
}
cols <- c("eam", "stan", "brms_r", "rwiener", "wienr", "rtdists")
regions <- list(
  all = rep(TRUE, nrow(g)),
  `u < 1e-3` = g$u < 1e-3,
  `1e-3 <= u <= 10` = g$u >= 1e-3 & g$u <= 10,
  `u > 10` = g$u > 10,
  `w in {0.001, 0.999}` = g$w %in% c(0.001, 0.999),
  `|v| = 20` = abs(g$v) == 20,
  `ref > -700` = g$ref > -700)
for (cl in cols) {
  cat("\n")
  for (nm in names(regions)) cat(summ(cl, regions[[nm]], nm), "\n")
}
cat("\n-- eam against Stan directly, where both finite\n")
ok <- is.finite(g$eam) & is.finite(g$stan)
cat(sprintf("n=%d maxabs=%.3e\n", sum(ok), max(abs(g$eam - g$stan)[ok])))
cat("\n-- worst 12 rows for stan\n")
e <- abs(g$stan - g$ref); e[!is.finite(e)] <- Inf
print(head(g[order(-e), c("t", "a", "w", "v", "up", "u", "ref", "eam",
                          "stan", "wienr")], 12), digits = 10)
cat("\n-- worst 12 rows for eam\n")
e <- abs(g$eam - g$ref); e[!is.finite(e)] <- Inf
print(head(g[order(-e), c("t", "a", "w", "v", "up", "u", "ref", "eam",
                          "stan", "wienr")], 12), digits = 10)
cat("\n-- worst 6 rows for rwiener\n")
e <- abs(g$rwiener - g$ref); e[!is.finite(e)] <- Inf
print(head(g[order(-e), c("t", "a", "w", "v", "up", "u", "ref", "rwiener",
                          "brms_r")], 6), digits = 10)
cat("\n-- the edge: decision time at and below zero\n")
for (td in c(0, -1e-3)) {
  cat("td =", td, " eam:", ddm_lpdf_both(td, 1, 1.4, 0.5, 1),
      " stan:", tryCatch(wiener_diffusion_lpdf(tau + td, 1L, 1.4, tau, 0.5,
                                                1),
                         error = function(e) paste("ERR",
                                                   conditionMessage(e))),
      "\n")
}
