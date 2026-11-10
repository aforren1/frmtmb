# Lane code8 fix round: an R profile for the suite runs that logs every
# call of se_check_at_fit() whose answer the gradient budget decides
# differently at 2 np (0.80.0) and at 4 np (the lane), so the runs show
# which test fits now defer their fit-time SE check.
local({
  log <- "C:/Users/adf44/source/r/frmtmb-wt-code8/dev/code8-out/budget-defer.log"
  setHook(packageEvent("frmtmb", "onLoad"), function(...) {
    ns <- asNamespace("frmtmb")
    f <- get("se_check_at_fit", ns)
    g <- function(fit) {
      out <- f(fit)
      np <- length(fit$opt$par)
      ev <- fit$opt$evals
      if (np > get("se_check_np_free", ns) && is.numeric(ev) &&
            length(ev) == 1L && is.finite(ev)) {
        share <- get("se_check_share", ns)
        old <- 2 * np <= share * ev
        if (!identical(old, out)) {
          cat(sprintf("%s\tnp %d\tevals %d\told %s\tnew %s\n",
                      paste(commandArgs(TRUE), collapse = " "), np, as.integer(ev), old,
                      out), file = log, append = TRUE)
        }
      }
      out
    }
    environment(g) <- environment()
    unlockBinding("se_check_at_fit", ns)
    assign("se_check_at_fit", g, envir = ns)
    lockBinding("se_check_at_fit", ns)
  })
})
