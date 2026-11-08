# Lane arstep: traces on gp_second_descent() and gp_marginal_alts(), sourced by dev/arstep-ng3.R and dev/arstep-r2ng3.R after rev_libs().
ns <- asNamespace("frmtmb"); sd_log <- character(0); sd_t0 <- 0
logit_ <- function(s) assign("sd_log", c(get("sd_log", globalenv()), s), envir = globalenv())
trace("gp_second_descent", tracer = quote(assign("sd_t0", proc.time()[[3]], envir = globalenv())),
      exit = quote({ r <- returnValue(); if (length(alt)) get("logit_", globalenv())(sprintf("ran first %.6f second %.6f kept %.6f secs %.2f", -opt$objective, if (exists("o") && !is.null(o)) -o$objective else NA, -r$objective, proc.time()[[3]] - get("sd_t0", globalenv()))) }),
      print = FALSE, where = ns)
# the working marginal's best value and each alternative's, from inside
# gp_marginal_alts(): how much lower the second maximum is
trace("gp_marginal_alts", exit = quote({ r <- returnValue(); if (length(r)) { av <- vapply(r, function(th) { i <- which(vapply(found, function(f) identical(f$theta, th), NA))[1]; found[[i]]$value }, 0); get("logit_", globalenv())(sprintf("marginal best %.4f alts %s", out$value, paste(sprintf("%.4f", av), collapse = ","))) } }),
      print = FALSE, where = ns)
