# False alarms and cost of the stationary-point escape, on the FIXED
# build. "Fires and gains nothing" is the false alarm; it costs two
# extra fits and no correctness.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")

# alpha truly 0, with the same covariate shape as the stall design so
# the marginal response is still skewed by the covariate.
make_sym <- function(seed, n = 200) {
  set.seed(1000 + seed)
  xs <- -abs(rnorm(n)) * 3
  data.frame(y = xs + rnorm(n, 0, 1.5), xs = xs)
}
make_mild <- function(seed, alpha = 1, n = 200) {
  set.seed(2000 + seed)
  xs <- -abs(rnorm(n)) * 3
  dl <- alpha / sqrt(1 + alpha^2)
  e <- 1.5 * (dl * abs(rnorm(n)) + sqrt(1 - dl^2) * rnorm(n))
  data.frame(y = xs + e, xs = xs)
}
make_strong <- function(seed, n = 200) {
  set.seed(3000 + seed)
  xs <- -abs(rnorm(n)) * 3
  data.frame(y = xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5, xs = xs)
}

one <- function(dd) {
  nw <- 0L
  t0 <- proc.time()[["elapsed"]]
  f <- withCallingHandlers(
    frm_sn(dd),
    warning = function(w) { nw <<- nw + 1L; invokeRestart("muffleWarning") })
  el <- proc.time()[["elapsed"]] - t0
  e <- f$opt[["stationary_escape"]]
  if (is.null(e)) e <- c(starts = 0, gain = 0)
  s <- sn_ll(dd)
  c(ll = as.numeric(logLik(f)), ll_sn = s[["ll"]],
    gap = as.numeric(logLik(f)) - s[["ll"]], alpha = f$opt$par[[4]],
    starts = e[["starts"]], gain = e[["gain"]], secs = el,
    conv = f$opt$convergence, warns = nw)
}

arms <- list(
  sym_n200 = function(s) make_sym(s, 200),
  sym_n50  = function(s) make_sym(s, 50),
  mild_a1  = function(s) make_mild(s, 1, 200),
  mild_a05 = function(s) make_mild(s, 0.5, 200),
  strong   = function(s) make_strong(s)
)
out <- list()
for (nm in names(arms)) {
  m <- as.data.frame(t(vapply(1:40, function(s) one(arms[[nm]](s)),
                              numeric(9))))
  m$arm <- nm
  m$seed <- 1:40
  out[[nm]] <- m
}
d <- do.call(rbind, out)
utils::write.csv(
  d,
  paste0("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/",
         "falsealarm-", Sys.getenv("SKEWINIT_TAG", "fixed"), ".csv"),
  row.names = FALSE)

cat("\n=== escape firing and what it bought, 40 seeds per arm\n")
cat(sprintf("%-9s %5s %6s %9s %9s %11s %9s %9s\n", "arm", "fires",
            "gain>0.01", "gain>1", "maxgain", "falsealarm", "med s fire",
            "med s no"))
for (nm in unique(d$arm)) {
  s <- d[d$arm == nm, ]
  fi <- s$starts > 0
  cat(sprintf("%-9s %5d %6d %9d %9.4f %11d %9.3f %9.3f\n", nm, sum(fi),
              sum(s$gain > 0.01), sum(s$gain > 1), max(s$gain),
              sum(fi & s$gain <= 0.01),
              stats::median(s$secs[fi]), stats::median(s$secs[!fi])))
}
cat("\ntotal fits", nrow(d), " total escape restarts", sum(d$starts), "\n")
cat("extra fits per fit, overall:", sum(d$starts) / nrow(d), "\n")
cat("\n=== against sn::selm on every arm\n")
cat("rows short of sn by more than 1e-6:", sum(d$gap < -1e-6), "\n")
cat("worst shortfall:", format(-min(d$gap), digits = 6), "\n")
cat("rows beating sn by more than 1e-6:", sum(d$gap > 1e-6),
    " best excess:", format(max(d$gap), digits = 6), "\n")
print(d[d$gap < -1e-6 | d$gap > 1e-6,
        c("arm", "seed", "alpha", "ll", "ll_sn", "gap", "starts", "gain")],
      row.names = FALSE)
