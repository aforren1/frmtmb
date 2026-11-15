# Lane parity2: brms 2.23.0 on a later nlf() or lf() on a parameter set
# earlier (no compile).
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
show <- function(f) {
  m <- character()
  r <- withCallingHandlers(tryCatch(f, error = function(e) e),
    message = function(c) { m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage") })
  if (inherits(r, "error")) cat("ERROR", conditionMessage(r), "\n") else
    cat("sigma:", deparse(r$pforms$sigma), " nl attr:", isTRUE(attr(r$pforms$sigma, "nl")), "\n")
  for (x in m) cat("  message:", trimws(x), "\n")
}
show(bf(y ~ x, sigma ~ z) + nlf(sigma ~ a * z))
show(bf(y ~ x) + nlf(sigma ~ a * z) + nlf(sigma ~ a + z))
show(bf(y ~ x) + nlf(sigma ~ a * z) + lf(sigma ~ z))
show(bf(y ~ x, sigma ~ z) + lf(sigma ~ w))
show(bf(y ~ x, sigma = 1) + lf(sigma ~ w))
