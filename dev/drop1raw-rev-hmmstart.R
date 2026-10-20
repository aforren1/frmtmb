# Reviewer: hmm_starts() on a fit whose call has `start = <local>`.
# rellib-r9's probe removed `start`; the lane's probe is update(), which
# keeps it. Arg: base | lane. Log: dev/drop1raw-rev-out/hmmstart-<arm>.txt
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.latent)})
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
small_data <- function(seed = 4501L, N = 10L, Tl = 20L) {
  set.seed(seed)
  G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(N), function(id) {
    s <- integer(Tl)
    s[1L] <- 1L
    for (t in seq_len(Tl)[-1L]) {
      s[t] <- sample.int(2L, 1L, prob = G[s[t - 1L], ])
    }
    data.frame(id = id, t = seq_len(Tl),
               y = stats::rnorm(Tl, c(0, 3)[s], 0.6))
  }))
}
dd <- small_data()
fam <- hmm(K = 2, gaussian(), time = t, group = id)
f0 <- q(frm(bf(y ~ 1), family = fam, data = dd))
fit <- local({
  st_gone <- f0$estimates
  q(frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t, group = id),
        data = dd, start = st_gone))
})
r <- tryCatch(q(hmm_starts(fit, n = 2, seed = 7L)), error = function(e) e)
cat("hmm_starts(), call has start = <gone local>:",
    if (inherits(r, "error")) paste("ERROR", substr(conditionMessage(r), 1,
                                                    200)) else
      paste("ok, status", paste(r$table$status, collapse = ",")), "\n")
