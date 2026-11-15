# Merge check 4 (gpfield x nlstop2) on the merged build: the gp() end
# warnings of frm() reach the caller on the sign-probe and pole-restart
# paths of nl_fit_probed(), once each. Every condition is kept in order.
#   Rscript dev/rel083-kt/gpfield/mc4-gp-nl.R
source("C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-kt/gpfield/gpfield-lib.R")
gf_libs("lane")
source("C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-kt/gpfield/gpfield-rows.R")
run <- function(tag, probe = TRUE, start = NULL) {
  old <- options(frmtmb.nl_sign_probe = probe)
  on.exit(options(old))
  r <- gf_row(tag)
  cond <- character(0)
  t0 <- proc.time()[[3]]
  f <- tryCatch(withCallingHandlers(
    if (is.null(start)) frm(r$fml, data = r$data, family = r$family) else
      frm(r$fml, data = r$data, family = r$family, start = start),
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
    cat(sprintf("MC4 %s probe=%s ERROR %s\n", tag, probe,
                substr(conditionMessage(f), 1, 200)))
    return(invisible())
  }
  gp <- grep("^W:frmtmb_gp", cond, value = TRUE)
  dup <- any(duplicated(gp))
  path <- if (any(cond == "M:frmtmb_nl_multistart")) "sign-probe" else
    if (any(cond %in% c("M:frmtmb_nl_pole_start", "W:frmtmb_nl_pole")))
      "pole-restart" else "first-fit"
  cat(sprintf("MC4 %s probe=%s path %s ll %.6f conv %d wall %.1f gp_warnings %d duplicated %s said %s\n",
              tag, probe, path, as.numeric(logLik(f)), f$opt$convergence,
              wall, length(gp), dup, isTRUE(f$cache$gp_end_said)))
  cat("   conditions in order:", paste(cond, collapse = " "), "\n")
  flush.console()
}
run("sweep:mm_V:bernoulli:40:1")
run("sweep:mm_V:bernoulli:40:1", probe = FALSE)
for (tg in c("sweep:mm_K:bernoulli:120:1", "sweep:mm_K:bernoulli:120:2",
             "sweep:mm_K:bernoulli:40:2", "sweep:mm_K:poisson:40:2")) {
  run(tg)
}
cat("DONE\n")
