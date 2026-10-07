# Lane setier: backlog "fixef() and summary() report NaN standard errors
# without a warning" (lane fixes' dev/fixes-sx-conv3.R design: gamSim
# eg 6, n 200, seeds 1 to 20, seven smooth formulas, 140 fits). Per
# fit: the optimizer code, whether fixef() has a non-finite SE, and what
# the user was told at the fit and at fixef() and summary().
#   Rscript dev/setier-sx.R <lib or "base"> <out .tsv>
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
if (identical(args[1], "src")) {
  suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-setier",
                                     quiet = TRUE))
} else suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), " BLAS probe (s):",
    system.time({m <- matrix(1, 600, 600); m %*% m})[["elapsed"]], "\n")
conds <- function(expr) {
  w <- character(0); m <- character(0)
  val <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  }, message = function(c) {
    m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage")
  })
  list(value = val, warnings = w, messages = m)
}
forms <- list(y ~ s(x1) + s(x2), y ~ s(x1, bs = "cr", k = 6),
              y ~ s(x1, by = g) + g, y ~ s(x1, by = z), y ~ s(x1, x2),
              y ~ t2(x1, x2), y ~ s(x0) + s(x1) + s(x2) + s(x3))
rows <- list()
for (seed in 1:20) {
  set.seed(seed)
  d <- suppressMessages(mgcv::gamSim(eg = 6, n = 200, scale = 2,
                                     verbose = FALSE))
  d$z <- runif(200)
  d$g <- factor(sample(c("a", "b", "c"), 200, TRUE))
  for (k in seq_along(forms)) {
    r <- conds(frm(bf(forms[[k]]), data = d))
    fx <- conds(fixef(r$value))
    sm <- conds(summary(r$value))
    se <- fx$value[, "Est.Error"]
    lost <- frmtmb:::sdr_of(r$value)$se_lost
    rows[[length(rows) + 1L]] <- data.frame(
      seed = seed, form = k, conv = r$value$opt$convergence,
      fixef_se_ok = all(is.finite(se)),
      fit_warn = length(r$warnings), fit_msg = length(r$messages),
      use_warn = length(c(fx$warnings, sm$warnings)),
      lost = paste(names(lost), lost, sep = ":", collapse = ";"),
      first = substr(c(r$warnings, r$messages, fx$warnings,
                       sm$warnings, "")[1L], 1, 90))
  }
}
res <- do.call(rbind, rows)
write.table(res, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf(paste0("fits %d; conv != 0: %d; non-finite fixef SE: %d; ",
                   "of those told nothing (no warning, no message): %d\n"),
            nrow(res), sum(res$conv != 0), sum(!res$fixef_se_ok),
            sum(!res$fixef_se_ok & res$fit_warn + res$fit_msg +
                  res$use_warn == 0)))
print(res[!res$fixef_se_ok | res$lost != "", 1:8])
