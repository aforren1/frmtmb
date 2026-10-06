# Reviewer, punch round 1: fixes' "c0 + exp(a)^k, a ~ 1 + (1 | g)" row
# showed 3 of 4 SEs NaN with no warning from frm() or fixef() on the
# trial merge. Repeats it, recording frm() warnings, whether the check
# was deferred, the warnings of the first SE use, and the SEs.
#   Rscript dev/nanse-rev2-c0k.R lane|merge [reps]
args <- commandArgs(TRUE)
arm <- args[1]
reps <- if (length(args) > 1) as.integer(args[2]) else 5L
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
set.seed(77)
G <- 8
dn <- data.frame(g = factor(rep(seq_len(G), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(G, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
for (r in seq_len(reps)) {
  w <- character()
  fit <- withCallingHandlers(
    frm(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1,
           nl = TRUE), data = dn),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  deferred <- !is.null(fit$cache$se_deferred)
  explained <- format(fit$cache$se_explained)
  wv <- character()
  se <- withCallingHandlers(fixef(fit)[, "Est.Error"], warning = function(x) {
    wv <<- c(wv, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  lost <- ns$sdr_of(fit)$se_lost
  cat(sprintf(paste0("rep %d: code %d, deferred %s, explained %s | frm() %s | ",
                     "fixef() %s | SE %s | lost %s\n"), r,
              fit$opt$convergence, deferred, explained,
              if (length(w)) substr(paste(w, collapse = " || "), 1, 70) else "none",
              if (length(wv)) substr(paste(wv, collapse = " || "), 1, 70) else "none",
              paste(signif(se, 3), collapse = " "),
              paste(sprintf("%s(%s)", names(lost), lost), collapse = " ")))
}
