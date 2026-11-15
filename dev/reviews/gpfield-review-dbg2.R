# Review of lane gpfield: why sweep:mm_V:bernoulli:40:1 (tk 60, importance
# sampling 22 below the Laplace value) gives no frmtmb_gp_laplace.
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-rows.R")
ns <- asNamespace("frmtmb")
tags <- commandArgs(TRUE)
for (tg in tags) {
  r <- gf_row(tg)
  calls <- 0L
  orig <- ns$gp_end_messages
  gf_put("gp_end_messages", function(fit) {
    calls <<- calls + 1L
    out <- orig(fit)
    cat("  gp_end_messages call", calls, "returned", length(out), "\n")
    e1 <- tryCatch(ns$gp_inner_pd(fit$obj, fit$opt$par, 1e-2),
                   error = function(e) paste("ERROR", conditionMessage(e)))
    e2 <- tryCatch(ns$gp_laplace_tk(fit),
                   error = function(e) paste("ERROR", conditionMessage(e)))
    cat("  at fit end: inner_pd", format(e1), "| tk", format(e2),
        "| par names", paste(unique(names(fit$opt$par)), collapse = ","),
        "| beta", paste(format(fit$estimates$beta, digits = 4),
                        collapse = ","), "\n")
    out
  })
  q <- gf_fit(r$fml, data = r$data, family = r$family)
  gf_put("gp_end_messages", orig)
  f <- q$fit
  cat("ROW", tg, "conv", f$opt$convergence, "nwarn", length(q$warn),
      "calls", calls, "| autoscale:", !is.null(f$autoscale),
      paste(names(f)[grepl("autos|scale", names(f))], collapse = ","), "\n")
  m <- ns$gp_end_messages(f)
  cat("  direct gp_end_messages:", length(m),
      vapply(m, `[[`, "", "class"), "\n")
  cat("  tk", ns$gp_laplace_tk(f), "\n")
  print(q$warn)
}
cat("final beta", format(f$estimates$beta, digits = 6), "\n")
str(f$autoscaled, max.level = 1)
