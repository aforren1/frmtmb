# Is the Gamma + ar1 spec's failure the weighted likelihood being
# unbounded? ar1(tt + 0 | g) has one row per (g, tt) cell here, so every
# row owns a random effect, and then the weighted likelihood grows like
# shape^(sum(w - 1) / 2) as shape grows. Arms vary the weights only.
#   Rscript dev/fuzzfix-gam-weights.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
d <- readRDS("dev/fuzzfix-out/data-gam.rds")
cat("rows:", nrow(d), " (g, tt) cells:", nrow(unique(d[c("g", "tt")])),
    " sum(w - 1):", format(sum(d$w - 1), digits = 6), "\n")
d$obs <- factor(seq_len(nrow(d)))
G <- Gamma(link = "log")
run <- function(lbl, fo, dd, profile, fam = G) {
  w <- character()
  f <- withCallingHandlers(
    tryCatch(frm(fo, family = fam, data = dd,
                 control = frmtmb_control(profile = profile)),
             error = function(e) e),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  if (inherits(f, "error")) {
    cat(sprintf("ARM %-28s ERROR %s\n", lbl,
                substr(conditionMessage(f), 1, 90)))
    return(invisible())
  }
  e <- f$estimates
  cat(sprintf(paste0("ARM %-28s logLik %14.4f betad %s theta %s",
                     " | nwarn %d\n"), lbl, as.numeric(logLik(f)),
              paste(format(e$betad, digits = 4), collapse = ","),
              paste(format(e$theta, digits = 4), collapse = ","),
              length(w)))
}
fo_ar <- bf(y | weights(ww) ~ 1 + x + ar1(tt + 0 | g))
fo_ol <- bf(y | weights(ww) ~ 1 + x + (1 | obs))
for (sc in c("w", "w/mean", "0.9", "1", "1.1")) {
  dd <- d
  dd$ww <- switch(sc, w = d$w, "w/mean" = d$w / mean(d$w),
                  rep(as.numeric(sc), nrow(d)))
  for (pr in c(FALSE, TRUE)) {
    run(sprintf("ar1 w=%s profile=%s", sc, pr), fo_ar, dd, pr)
    run(sprintf("obs w=%s profile=%s", sc, pr), fo_ol, dd, pr)
  }
}
# the same unboundedness on a gaussian response: sigma -> 0
dd <- d; dd$ww <- rep(1.1, nrow(d)); dd$y <- log(d$y)
run("gaussian obs w=1.1", fo_ol, dd, FALSE, gaussian())
dd$ww <- rep(0.9, nrow(d))
run("gaussian obs w=0.9", fo_ol, dd, FALSE, gaussian())
cat("END\n")
