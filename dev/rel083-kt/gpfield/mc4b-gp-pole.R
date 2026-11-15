# Merge check 4, pole paths: gp() end warnings on the pole-restart paths
# of nl_fit_probed() on the merged build.
#  (a) given start between the poles (frmtmb_nl_pole path): gpfield's
#      mm_V row with K started inside the band of -w.
#  (b) default start: nlstop2's binomial mm data (seeds 11, 44, 55; the
#      rows whose first fit ends between the poles), with gp(x) added to Vm.
#   Rscript dev/rel083-kt/gpfield/mc4b-gp-pole.R   (run from the tree root)
source("C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-kt/gpfield/gpfield-lib.R")
gf_libs("lane")
source("C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-kt/gpfield/gpfield-rows.R")
src <- parse("dev/nlstop-fresh.R")
keep <- vapply(src, function(e) {
  is.call(e) && identical(e[[1L]], as.name("<-")) &&
    as.character(e[[2L]])[1] %in% c("F", "pos", "B", "gen", "form")
}, NA)
for (e in src[keep]) eval(e)
run <- function(label, expr) {
  cond <- character(0)
  t0 <- proc.time()[[3]]
  f <- tryCatch(withCallingHandlers(expr,
    warning = function(w) {
      cond <<- c(cond, paste0("W:", class(w)[1]))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cond <<- c(cond, paste0("M:", class(m)[1]))
      invokeRestart("muffleMessage")
    }), error = function(e) e)
  wall <- proc.time()[[3]] - t0
  if (inherits(f, "error")) {
    cat(sprintf("MC4B %s ERROR %s\n", label, substr(conditionMessage(f), 1, 200)))
    return(invisible())
  }
  gp <- grep("^W:frmtmb_gp", cond, value = TRUE)
  b <- f$estimates$beta
  cat(sprintf("MC4B %s ll %.6f conv %d wall %.1f gp_warnings %d duplicated %s said %s\n",
              label, as.numeric(logLik(f)), f$opt$convergence, wall,
              length(gp), any(duplicated(gp)), isTRUE(f$cache$gp_end_said)))
  cat("   beta:", paste(names(b), sprintf("%.5g", b), collapse = " "), "\n")
  cat("   conditions in order:", paste(cond, collapse = " "), "\n")
  flush.console()
}
# (a) given start inside the pole band
r <- gf_row("sweep:mm_V:bernoulli:40:1")
pt <- par_template(r$fml, data = r$data, family = r$family)
bn <- names(pt$beta)
cat("beta names:", bn, "| range w:", range(r$data$w), "\n")
kn <- grep("^K", bn, value = TRUE)[1]
for (K in c(-2.5, -1)) {
  st <- pt$beta
  st[] <- 0
  st[grep("^V", bn)] <- 1
  st[kn] <- K
  run(sprintf("mm_V given start %s=%g", kn, K),
      frm(r$fml, data = r$data, family = r$family, start = list(beta = st)))
}
# (b) default start, binomial mm with a gp() on Vm
for (seed in c(11, 44, 55)) {
  d <- gen("binomial", "mm", seed * 1000 + match("mm", names(B)))
  fm <- bf(y | trials(nt) ~ Vm * xp / (K + xp), Vm ~ 1 + gp(x), K ~ 1,
           nl = TRUE)
  run(sprintf("binomial mm %d + gp(x) on Vm, default start", seed),
      frm(fm, family = binomial(), data = d))
}
cat("DONE\n")
