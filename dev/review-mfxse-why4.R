# Review of lane mfxse: avg_slopes() at new levels in a fresh process, every
# warning with its class, and interop_rows_check() instrumented (in this
# process only) to print what it found. Data seed 7411.
.libPaths(c("C:/Users/adf44/source/r/wt-mfxse-lib", "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
set.seed(7411)
J <- 12L; m <- 20L; n <- J * m
d <- data.frame(x = rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)),
                g = factor(rep(seq_len(J), each = m)))
d$cnt <- rpois(n, exp(0.5 + 0.4 * d$x + rnorm(J, 0, 0.8)[d$g]))
fit <- suppressWarnings(frm(cnt ~ x + f + (1 | g), family = poisson(), data = d))
ndn <- data.frame(x = c(0, 1, 2), f = "a", g = factor(c("3", "new", "new2")))
orig <- ns$interop_rows_check
inst <- function(model, newdata, dots, asked) {
  armed <- isTRUE(model$cache[["mfx_rows_armed"]])
  r <- withCallingHandlers(orig(model, newdata, dots, asked), warning = function(w) {
    cat("  check raised:", class(w)[1], "\n")
  })
  if (armed) cat("  armed call: rows", nrow(newdata), "asked", sum(asked), "\n")
  r
}
environment(inst) <- environment()
utils::assignInNamespace("interop_rows_check", inst, ns = "frmtmb")
for (k in c("avg_slopes", "predictions")) {
  cat(k, ":\n")
  ws <- list()
  invisible(withCallingHandlers(
    if (k == "avg_slopes") avg_slopes(fit, newdata = ndn, variables = "x", allow_new_levels = TRUE)
    else predictions(fit, newdata = ndn, allow_new_levels = TRUE),
    warning = function(w) { ws[[length(ws) + 1]] <<- w; invokeRestart("muffleWarning") }))
  for (w in ws) cat("  reached the caller:", class(w)[1], "|", substr(conditionMessage(w), 1, 60), "\n")
}
