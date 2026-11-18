# Lane nlband fix round: which existing check tells a separated hole fit
# (bernoulli band data, y = 1 for x > 2, review item 1) from a hole
# maximum (band data)? For the fit the current build returns: the limit
# check nl_limit_message() as is, the same without its gate on one
# conditional standard error, the fit's standard errors, and the SE
# check's analysis (tier, lost).
# NLBAND_ARM=after Rscript dev/nlband-dbg3.R
source("dev/nlband-common.R")
nogate <- ns$nl_limit_message
b <- body(nogate)
txt <- deparse(nogate)
txt <- sub("if (!all(abs(D[big]) * sqrt(abs(diag(H))[big]) >= 1)) return(NULL)",
           "", txt, fixed = TRUE)
nogate <- eval(parse(text = txt))
environment(nogate) <- ns
probe <- function(tag, f) {
  if (iserr(f)) return(cat(tag, "ERROR", f$msg, "\n"))
  v <- tryCatch(suppressWarnings(vcov(f)), error = function(e) NULL)
  an <- tryCatch(ns$se_fit_analysis(f, force = TRUE), error = function(e) NULL)
  m1 <- tryCatch(ns$nl_limit_message(f), error = function(e) "ERR")
  m2 <- tryCatch(nogate(f), error = function(e) paste("ERR", conditionMessage(e)))
  cat(sprintf("%s | ll %.6g code %d beta %s | se %s | rownames %s | tier %s lost %s | limit %s | nogate %s | msg %s warn %s\n",
              tag, llv(f), f$opt$convergence,
              paste(signif(f$estimates$beta, 6), collapse = " "),
              paste(signif(sqrt(diag(v)), 4), collapse = " "),
              paste(rownames(v), collapse = ","),
              if (is.null(an)) "NULL" else an$tier,
              if (is.null(an)) "-" else paste(an$lost, collapse = ","),
              if (is.null(m1)) "NULL" else substr(m1, 1, 60),
              if (is.null(m2)) "NULL" else substr(m2, 1, 60), msum(f), wsum(f)))
}
for (s in 2:4) {
  d <- band_data(s)
  d$y <- as.integer(d$x > 2)
  probe(paste("SEP", s), tt(frm(band_form, data = d, family = bernoulli())))
}
for (s in 1:2) probe(paste("BAND", s), tt(frm(band_form, data = band_data(s))))
for (s in 1:2) {
  d <- band_data(s)
  d$y1 <- rbinom(nrow(d), 1, plogis(2 / (d$x - 2)))
  probe(paste("BERN", s), tt(frm(bf(y1 ~ a / (x + b), a + b ~ 1, nl = TRUE),
                                  family = bernoulli(), data = d)))
}
cat("DONE\n")
