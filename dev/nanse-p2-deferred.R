# Lane, punch round 2: dev/nanse-rev2-deferred.R with the instrument
# moved to the new decision (se_check_at_fit() always FALSE). When the fit-time Hessian
# is over budget the warning is promised at the first SE use. Instrument:
# fd_hessian() is made to report "over budget" whenever it has a budget,
# in this process only (assignInNamespace), so every fit with random
# effects defers, as a large mixed model or a loaded machine does. Model: fixes'
# "c0 + exp(a)^k, a ~ 1 + (1 | g)", which loses 4 of 6 SEs. Per first
# accessor: SE warnings from it, then from vcov() and summary().
#   Rscript dev/nanse-rev2-deferred.R lane|merge
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
assignInNamespace("se_check_at_fit", function(fit) FALSE, "frmtmb")
set.seed(77)
G <- 8
dn <- data.frame(g = factor(rep(seq_len(G), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(G, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(v = v, w = sum(grepl("Standard errors are not available", w)))
}
first <- list(
  summary = function(f) summary(f),
  fixef = function(f) fixef(f),
  vcov = function(f) vcov(f),
  confint = function(f) confint(f),
  VarCorr = function(f) VarCorr(f)
)
for (nm in names(first)) {
  r <- cap(frm(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1,
                  nl = TRUE), data = dn))
  fit <- r$v
  deferred <- !is.null(fit$cache$se_deferred)
  a <- cap(first[[nm]](fit))
  b <- cap(vcov(fit))
  s <- cap(summary(fit))
  out <- capture.output(print(s$v))
  lost <- frmtmb:::sdr_of(fit)$se_lost
  cat(sprintf(paste0("first use %-10s: deferred %s | SE warnings frm() %d, ",
                     "%s() %d, then vcov() %d, summary() %d | summary ",
                     "lists them %s | lost %d\n"),
              nm, deferred, r$w, nm, a$w, b$w, s$w,
              any(grepl("without a standard error", out)), length(lost)))
}
