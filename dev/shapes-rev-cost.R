# Reviewer, priority 5: what does the new fitted()/residuals() COST?
#
# The clock on this box ticks at 10 ms and swings with load
# (dev/lane-rules.md), so the primary instrument here is a COUNT, not a
# time: how many times `autoscale_sdreport()` runs per call. A count is
# load-independent and it answers the question the brief asks, which is
# whether a hot path now pays repeatedly. Times are reported beside it,
# interleaved, minimum of several rounds, with a control that must
# report 1.0.
#
#   Rscript dev/shapes-rev-cost.R lane
#   Rscript dev/shapes-rev-cost.R base

arg <- commandArgs(trailingOnly = TRUE)
which_lib <- if (identical(arg[1], "base")) "base" else "lane"
lib <- if (which_lib == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/shapes-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", lib, "\n")

CNT <- new.env(parent = emptyenv())
CNT$n <- 0L
ns <- asNamespace("frmtmb")
orig <- get("autoscale_sdreport", envir = ns)
unlockBinding("autoscale_sdreport", ns)
assign("autoscale_sdreport", function(...) {
  CNT$n <- CNT$n + 1L
  orig(...)
}, envir = ns)
lockBinding("autoscale_sdreport", ns)

reset <- function() CNT$n <<- 0L
count <- function(lab, e) {
  reset()
  t0 <- proc.time()
  v <- tryCatch(suppressWarnings(e), error = function(c)
    paste0("ERROR: ", substr(conditionMessage(c), 1, 60)))
  el <- (proc.time() - t0)[3]
  cat(sprintf("  %-46s sdreports=%3d  %6.2f s  %s\n", lab, CNT$n, el,
              if (is.character(v)) v else
                paste(class(v)[1], paste(dim(v) %||% length(v),
                                         collapse = "x"))))
  invisible(CNT$n)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

set.seed(4242)
N <- 4000L
dl <- data.frame(x = rnorm(N), z = rnorm(N),
                 g = factor(rep(seq_len(200), length.out = N)))
dl$y <- rnorm(N, 1 + 0.5 * dl$x + rnorm(200, 0, 0.7)[dl$g], 1)
dl$ord <- factor(cut(1 + 0.8 * dl$x + rnorm(N), c(-Inf, -0.3, 0.8, Inf),
                     labels = 1:3), ordered = TRUE)

cat("\n== large mixed fit, n =", N, ", 200 groups ==\n")
t0 <- proc.time()
fb <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dl)
cat("  fit took", round((proc.time() - t0)[3], 2), "s\n")
count("fitted(fit)  [first call]", fitted(fb))
count("fitted(fit)  [second call, cache warm]", fitted(fb))
count("residuals(fit)", residuals(fb))
count("frm_linpred(fit, type='response')", frm_linpred(fb, type = "response"))
count("fixef(fit)", fixef(fb))
count("fixef_by_dpar(fit)", if (which_lib == "lane")
  frmtmb::fixef_by_dpar(fb) else fixef(fb))
count("predict(fit, ndraws = 25)", predict(fb, ndraws = 25))
count("conditional_effects(fit)", conditional_effects(fb))
count("summary(fit)", summary(fb))
count("vcov(fit)", vcov(fb))

cat("\n== a FRESH fit each time (no warm cache) ==\n")
mk <- function() frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dl)
count("fitted(fresh fit)", fitted(mk()))
count("frm_linpred(fresh fit, 'response')",
      frm_linpred(mk(), type = "response"))
count("plot(fresh fit) [scatter path]", {
  f <- mk(); grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(f)
})

cat("\n== ordinal fit: the finite-difference delta method ==\n")
fo <- frm(bf(ord ~ x + z) + cumulative(), data = dl[seq_len(600), ])
cat("  outer parameters:", length(fo$opt$par), "\n")
count("fitted(ordinal)", fitted(fo))
count("fitted(ordinal) [second call]", fitted(fo))
count("residuals(ordinal)", residuals(fo))
count("predict(ordinal, ndraws = 25)", predict(fo, ndraws = 25))
count("conditional_effects(ordinal)", conditional_effects(fo))

cat("\n== the control: a call whose sdreport count must be 1 then 0 ==\n")
fc <- frm(bf(y ~ x) + gaussian(), data = dl[seq_len(400), ])
count("confint(fresh small fit)   [control, expect 1]", confint(fc))
count("confint(same fit again)    [control, expect 0]", confint(fc))
