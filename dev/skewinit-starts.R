# Which START rule removes the stall, measured on the BASE build so the
# escape refit cannot mask the answer. Four arms over both 40-seed
# streams, 80 fits each. Every arm passes BOTH betad entries, so the
# sigma start each number came from is explicit.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/rellib-r3")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")

rule <- function(m3) {
  if (!is.finite(m3) || m3 == 0) m3 <- 1
  2 * sign(m3) + 0.5 * m3
}

arms <- list(
  raw_sdy = function(dd, r) c(log(stats::sd(dd$y)), rule(skew(dd$y))),
  res_sdy = function(dd, r) c(log(stats::sd(dd$y)), rule(skew(r))),
  res_sdr = function(dd, r) c(log(stats::sd(r)), rule(skew(r))),
  raw_sdr = function(dd, r) c(log(stats::sd(r)), rule(skew(dd$y)))
)

one <- function(dd) {
  r <- residuals(stats::lm(y ~ xs, data = dd))
  ll_ref <- sn_ll(dd)[["ll"]]
  out <- c(ll_sn = ll_ref)
  for (nm in names(arms)) {
    s <- arms[[nm]](dd, r)
    f <- frm_sn(dd, start = list(betad = s))
    out[[nm]] <- as.numeric(logLik(f)) - ll_ref
    out[[paste0(nm, "_a")]] <- f$opt$par[[4]]
  }
  # the default path, for the record: same as raw_sdy on the base build
  out[["default"]] <- as.numeric(logLik(frm_sn(dd))) - ll_ref
  out
}

res <- list()
for (dead in c(FALSE, TRUE)) {
  m <- t(vapply(1:40, function(s) one(make_data(s, dead)), one(make_data(1, FALSE))))
  d <- as.data.frame(m)
  d$seed <- 1:40
  d$stream <- if (dead) "dead_draw" else "stated"
  res[[length(res) + 1L]] <- d
}
d <- do.call(rbind, res)
utils::write.csv(
  d, "C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/starts.csv",
  row.names = FALSE)

cat("\nstalls (logLik below sn::selm by more than 1e-6), 40 seeds each\n")
cat(sprintf("%-10s %-10s %8s %12s %12s\n", "stream", "arm", "stalls",
            "worst gap", "max |a| kept"))
for (st in unique(d$stream)) {
  s <- d[d$stream == st, ]
  for (nm in c("default", names(arms))) {
    bad <- s[[nm]] < -1e-6
    cat(sprintf("%-10s %-10s %8d %12.6f %12.4f\n", st, nm, sum(bad),
                if (any(bad)) -min(s[[nm]][bad]) else 0,
                if (nm %in% names(arms)) max(abs(s[[paste0(nm, "_a")]])) else
                  NA_real_))
  }
}
cat("\nseeds that stall, by arm\n")
for (st in unique(d$stream)) {
  s <- d[d$stream == st, ]
  for (nm in c("default", names(arms))) {
    bad <- s[[nm]] < -1e-6
    if (any(bad)) {
      cat(sprintf("%-10s %-10s %s\n", st, nm,
                  paste(s$seed[bad], collapse = " ")))
    }
  }
}
cat("\nsigma start against the fitted sigma, arm res_sdy, stream stated\n")
s <- d[d$stream == "stated", ]
dd1 <- make_data(1, FALSE)
cat("seed 1: log(sd(y)) =", format(log(sd(dd1$y)), digits = 6),
    " log(sd(resid)) =",
    format(log(sd(residuals(lm(y ~ xs, data = dd1)))), digits = 6),
    " fitted sigma_Intercept =",
    format(frm_sn(dd1, start = list(betad = arms$res_sdy(
      dd1, residuals(lm(y ~ xs, data = dd1)))))$opt$par[[3]], digits = 6),
    "\n")
