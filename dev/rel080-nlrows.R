# The three fresh rows of the nlstart review where 0.80.0 reads lower than
# the nlstart lane's build (dev/rel080-log.md): each fit default, with
# bigx's optimizer basis off (options(frmtmb.opt_basis = FALSE)), and with
# the sign probe off, to tell which change moved it.
#
#   REV_ARM=<lib> Rscript dev/rel080-nlrows.R
src <- readLines("dev/review-nlstart-fresh2.R")
stop_at <- grep("^cat\\(\"ROWHEAD", src)
eval(parse(text = src[1:(stop_at - 1L)]))
rows <- list(c("gaussian", "gomp", 22), c("poisson", "two_sat_rev", 33),
             c("bernoulli", "mm", 11))
one <- function(fm, fam, d, opts) {
  ms <- character()
  r <- withr::with_options(opts, tryCatch(withCallingHandlers(
    frm(fm, family = F[[fam]], data = d),
    warning = function(w) {
      ms <<- c(ms, paste0("W:", class(w)[1L]))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      ms <<- c(ms, paste0("M:", class(m)[1L]))
      invokeRestart("muffleMessage")
    }), error = function(e) NULL))
  if (is.null(r)) return("error")
  sprintf("ll %.6f conv %d | %s", as.numeric(logLik(r)), r$opt$convergence,
          paste(unique(ms), collapse = " "))
}
for (rw in rows) {
  fam <- rw[1]; body <- rw[2]; seed <- as.numeric(rw[3])
  d <- gen(fam, body, seed * 1000 + match(body, names(B)))
  fm <- form(fam, body)
  cat(fam, body, seed, "\n")
  cat("  default      ", one(fm, fam, d, list()), "\n")
  cat("  no basis     ", one(fm, fam, d, list(frmtmb.opt_basis = FALSE)), "\n")
  cat("  probe off    ", one(fm, fam, d, list(frmtmb.nl_sign_probe = FALSE)), "\n")
}
cat("NLROWS DONE\n")
