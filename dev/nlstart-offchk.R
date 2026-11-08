# Why do 7 ordinal logi rows differ by up to 6e-6 between the probe-off fit
# and 0.79.1? Check whether nl_par_units() fires on them, and refit twice.
source("dev/nlstart-common.R")
fam <- "cumulative"
src <- readLines("dev/nlstart-scan.R")
eval(parse(text = src[grep("^llv <-", src):(grep("^cat[(]\"ROWHEAD", src) - 1L)]))
options(frmtmb.nl_sign_probe = FALSE)
for (r in list(c(2, -0.6), c(4, -0.6), c(5, 0.6))) {
  d <- gen(r[1] * 100 + 1, r[2])
  fm <- bf(y ~ b1 / (1 + exp(b2 * x)), b1 + b2 ~ 1, nl = TRUE)
  f <- suppressWarnings(suppressMessages(frm(fm, family = F[[fam]], data = d)))
  f2 <- suppressWarnings(suppressMessages(frm(fm, family = F[[fam]], data = d)))
  cat(sprintf("seed %d slope %.1f ll %.9f again %.9f par_units %s beta %s\n",
              r[1], r[2], as.numeric(logLik(f)), as.numeric(logLik(f2)),
              paste(format(f$par_units %||% 1), collapse = ","),
              paste(format(f$estimates$beta, digits = 6), collapse = ",")))
}
