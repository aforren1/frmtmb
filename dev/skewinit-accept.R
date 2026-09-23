# Acceptance: both 40-seed streams, 80 fits, every one against the
# optimum sn::selm reports. Also the cost: how many escape restarts ran
# and what they bought.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")

esc <- function(f) {
  e <- f$opt[["stationary_escape"]]
  if (is.null(e)) c(starts = 0, gain = 0) else e
}

one <- function(dd, seed, stream) {
  # conditions and the optimizer's own code are part of acceptance: a
  # fit that reaches the optimum and warns is still a defect
  nw <- 0L
  nm <- 0L
  t0 <- proc.time()[["elapsed"]]
  f <- withCallingHandlers(
    frm_sn(dd),
    warning = function(w) { nw <<- nw + 1L; invokeRestart("muffleWarning") },
    message = function(m) { nm <<- nm + 1L; invokeRestart("muffleMessage") })
  el <- proc.time()[["elapsed"]] - t0
  e <- esc(f)
  s <- sn_ll(dd)
  c(seed = seed, ll = as.numeric(logLik(f)), ll_sn = s[["ll"]],
    gap = as.numeric(logLik(f)) - s[["ll"]],
    alpha = f$opt$par[[4]], alpha_sn = s[["alpha"]],
    starts = e[["starts"]], gain = e[["gain"]], secs = el,
    conv = f$opt$convergence, warns = nw, msgs = nm,
    grad = max(abs(f$obj$gr(f$opt$par))))
}

res <- list()
for (dead in c(FALSE, TRUE)) {
  stream <- if (dead) "dead_draw" else "stated"
  m <- t(vapply(1:40, function(s) one(make_data(s, dead), s, stream),
                numeric(13)))
  d <- as.data.frame(m)
  d$stream <- stream
  res[[length(res) + 1L]] <- d
}
d <- do.call(rbind, res)
utils::write.csv(
  d, paste0("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/",
            "accept-", Sys.getenv("SKEWINIT_TAG", "fixed"), ".csv"),
  row.names = FALSE)

op <- options(digits = 12)
cat("\n=== ACCEPTANCE, 80 fits against sn::selm\n")
for (st in unique(d$stream)) {
  s <- d[d$stream == st, ]
  cat(sprintf("%-10s n=%2d  worst shortfall below sn %.3e  best excess %.3e\n",
              st, nrow(s), -min(s$gap), max(s$gap)))
  cat(sprintf("%-10s  |gap| max %.3e   fits reaching sn within 1e-6: %d\n",
              "", max(abs(s$gap)), sum(abs(s$gap) < 1e-6)))
  cat(sprintf("%-10s  escape ran on %d of %d, gained > 1e-6 on %d\n",
              "", sum(s$starts > 0), nrow(s), sum(s$gain > 1e-6)))
}
cat("\nworst |gap| over all 80:", max(abs(d$gap)), "\n")
cat("max relative |gap| / |ll_sn|:", max(abs(d$gap) / abs(d$ll_sn)), "\n")
cat("rows where frmtmb BEATS sn by more than 1e-6:", sum(d$gap > 1e-6), "\n")
cat("total escape restarts over 80 fits:", sum(d$starts), "\n")
options(op)
print(round(d[, c("seed", "gap", "alpha", "alpha_sn", "starts",
                  "gain", "conv", "warns", "msgs", "grad")], 6),
      row.names = FALSE)
