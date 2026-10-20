# test-review-v29.R "a boundary correlation reports NA bounds, not a
# clamp" failed on the merged build (dev/ciharden-log/r072-ref-suite/).
# The test's three cases, with what each build says: the fit's code,
# the boundary messages, confint_varcorr()'s warnings, the sd rows.
#   Rscript dev/rel072-v29-diag.R <lib>
.libPaths(c(commandArgs(TRUE)[1], "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
set.seed(61)
ng <- 12
u <- stats::rnorm(ng, 0, 0.8)
per <- 5
dd <- data.frame(id = factor(rep(paste0("a", seq_len(ng)), each = per)),
                 x = stats::rnorm(ng * per))
dd$y1 <- 1 + 0.5 * dd$x + u[as.integer(dd$id)] + stats::rnorm(nrow(dd), 0, 0.4)
dd$y2 <- -1 + 0.2 * dd$x + 0.6 * u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.3)
dk <- dd
set.seed(8004)
for (v in c("x", "y1", "y2")) {
  dk[[v]] <- dk[[v]] * (1 + sample(0:4, nrow(dk), TRUE) * 2^-52)
}
cases <- list(list(d = dd, t0 = NULL), list(d = dd, t0 = 500),
              list(d = dk, t0 = 2000))
for (i in seq_along(cases)) {
  cs <- cases[[i]]
  msgs <- character(0); ws <- character(0)
  st <- if (!is.null(cs$t0)) list(theta = c(0, 0, cs$t0))
  fit <- withCallingHandlers(
    frm(mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
             bf(y2 ~ x + (1 | q | id)) + gaussian()), data = cs$d,
        start = st),
    warning = function(w) { ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning") },
    message = function(m) { msgs <<- c(msgs, class(m)[1], conditionMessage(m))
      invokeRestart("muffleMessage") })
  cat(sprintf("\n== case %d (t0 %s): code %d, theta %s, logLik %.8f\n", i,
              format(cs$t0), fit$opt$convergence,
              paste(format(fit$opt$par[grep("theta", names(fit$opt$par))],
                           digits = 6), collapse = " "),
              as.numeric(logLik(fit))))
  cat("fit warnings:", if (length(ws)) paste(substr(ws, 1, 300), collapse = " || ") else "none", "\n")
  cat("fit messages:", if (length(msgs)) paste(substr(msgs, 1, 300), collapse = " || ") else "none", "\n")
  cw <- character(0)
  ci <- withCallingHandlers(confint_varcorr(fit), warning = function(w) {
    cw <<- c(cw, conditionMessage(w)); invokeRestart("muffleWarning") })
  cat("confint_varcorr warnings:", if (length(cw)) paste(substr(cw, 1, 400), collapse = " || ") else "none", "\n")
  print(ci)
  lost <- frmtmb:::sdr_of(fit)$se_lost
  cat("se_lost:", if (length(lost)) paste(names(lost), lost) else "none", "\n")
}
