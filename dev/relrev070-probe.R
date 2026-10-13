# Reviewer of the 0.70.0 consolidation: details behind two rows of
# dev/relrev070-interplay.R. (a) `y ~ mo(x) + (mo(x) | g)` fits at a
# nonzero code: the optimizer's message and every condition, messages
# included. (b) `y ~ mo(x2)` with category 2 never observed (one exactly
# flat simplex direction): the simplex SEs and the conditions.
#   Rscript dev/relrev070-probe.R r8|lane
arm <- commandArgs(TRUE)[1]
.libPaths(c(switch(arm, r8 = "C:/Users/adf44/source/r/rellib-r8",
                   lane = "C:/Users/adf44/source/r/wt-moparity-lib"),
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, as.character(packageVersion("frmtmb")), "\n")
src <- readLines("dev/relrev070-interplay.R")
i <- grep("^mk <- function", src)
j <- grep("^designs <- list", src)
eval(parse(text = src[i:(j - 1L)]))
cond <- function(expr) {
  w <- character(); m <- character()
  v <- withCallingHandlers(expr,
    warning = function(x) { w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning") },
    message = function(x) { m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage") })
  list(v = v, w = w, m = m)
}
for (seed in c(1, 3, 8, 10)) {
  d <- mk(seed)
  r <- cond(frm(bf(y ~ mo(x) + (mo(x) | g)), family = gaussian(), data = d))
  s <- cond(summary(r$v))
  th <- r$v$opt$par[names(r$v$opt$par) == "theta"]
  cat(sprintf("mo_re seed %d: code %d '%s' | theta %s\n", seed,
              r$v$opt$convergence, r$v$opt$message,
              paste(format(th, digits = 3), collapse = " ")))
  for (x in unique(c(r$w, s$w))) cat("   W:", substr(x, 1, 150), "\n")
  for (x in unique(c(r$m, s$m))) cat("   M:", substr(x, 1, 150), "\n")
}
for (seed in c(10, 18, 1)) {
  d <- mk(seed)
  r <- cond(frm(bf(y ~ mo(x2)), family = gaussian(), data = d))
  s <- cond(summary(r$v))
  v <- suppressWarnings(tryCatch(vcov(r$v, full = TRUE), error = function(e) NULL))
  sdr <- suppressWarnings(frmtmb:::sdr_of(r$v))
  ss <- summary(sdr, "fixed")
  z <- ss[grepl("zeta", rownames(ss)), , drop = FALSE]
  cat(sprintf("dir_flat seed %d: code %d | zeta est %s | zeta SE %s\n", seed,
              r$v$opt$convergence,
              paste(format(z[, 1], digits = 4), collapse = " "),
              paste(format(z[, 2], digits = 4), collapse = " ")))
  for (x in unique(c(r$w, s$w))) cat("   W:", substr(x, 1, 150), "\n")
  for (x in unique(c(r$m, s$m))) cat("   M:", substr(x, 1, 150), "\n")
}
