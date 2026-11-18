# Lane drawdir: what firing the draw check does on a fit, against a known
# predictive truth. For each case id (dev/drawdir-cases.R) with newdata
# and truth (the wb designs G6, R*, S*): predict() at the case's newdata
# (new level), 4000 draws, seed 1, as it comes and with the check's
# verdict forced to fire; the exact probability the true predictive law
# puts inside the central 90% and 98% intervals, the interval widths, and
# what was said (classes) with the Pareto k.
# Usage: Rscript dev/drawdir-force.R <lib or "base"> <id> [<id> ...]
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
source("dev/drawdir-cases.R")
run <- function(fit, nd, force) {
  if (force) {
    chk <- fit$cache[["draw_check"]]
    chk$fire <- TRUE
    assign("draw_check", chk, envir = fit$cache)
    if (exists("draw_check_said", envir = fit$cache, inherits = FALSE)) {
      rm("draw_check_said", envir = fit$cache)
    }
  }
  said <- character(0)
  set.seed(1)
  p <- withCallingHandlers(
    predict(fit, newdata = nd, allow_new_levels = TRUE, ndraws = 4000,
            summary = FALSE),
    warning = function(w) {
      said <<- c(said, paste(class(w)[1],
                             sub(".*Pareto k is ([0-9.]+).*", "k=\\1",
                                 conditionMessage(w))))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      said <<- c(said, paste(class(m)[1],
                             sub(".*Pareto k ([0-9.]+).*", "k=\\1",
                                 conditionMessage(m))))
      invokeRestart("muffleMessage")
    })
  list(s = p[, 1], said = paste(unique(said), collapse = ";"))
}
for (id in args[-1]) {
  case <- case_fit(id)
  fit <- case$fit
  tr <- case$truth
  out <- character(0)
  for (force in c(FALSE, TRUE)) {
    rr <- tryCatch(run(fit, case$nd, force), error = function(e) NULL)
    s <- rr$s
    if (is.null(s)) {
      out <- c(out, "error")
      next
    }
    cov <- function(a) {
      q <- stats::quantile(s, c(a / 2, 1 - a / 2), names = FALSE)
      c(stats::pnorm(q[2], tr$mean, tr$sd) - stats::pnorm(q[1], tr$mean, tr$sd),
        q[2] - q[1])
    }
    c90 <- cov(0.1)
    c98 <- cov(0.02)
    chk <- fit$cache[["draw_check"]]
    out <- c(out, sprintf("%s cov90 %.3f w90 %.3f cov98 %.3f w98 %.3f [%s]",
                          if (force) "forced" else "as is ", c90[1], c90[2],
                          c98[1], c98[2], rr$said))
  }
  cat(sprintf("%-10s %s | %s\n", id, out[1], out[2]))
}
