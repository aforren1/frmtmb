# Lane wcheck2: every row right censored (backlog 4998's actual model),
# asym_laplace and gaussian with (1 | g4), every warning in full. Data of
# dev/wcheck-cases.R part cens (seed 20261021), k = 80.
#   Rscript dev/wcheck2-allcens.R <base|lane|lib>
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
set.seed(20261021)
n <- 80
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)),
                g4 = factor(rep(1:20, each = 4)))
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$ya <- 1 + 0.5 * d$x + rnorm(n)
d$ca <- "right"
d$yac <- rep(min(d$ya), n)
for (fam in list(asym_laplace(), gaussian())) {
  ws <- character()
  f <- withCallingHandlers(frm(yac | cens(ca) ~ x + (1 | g4), family = fam,
                               data = d),
                           warning = function(w) {
                             ws <<- c(ws, paste(class(w)[1],
                                                conditionMessage(w)))
                             invokeRestart("muffleWarning")
                           }, message = function(m) {
                             ws <<- c(ws, paste("MESSAGE",
                                                conditionMessage(m)))
                             invokeRestart("muffleMessage")
                           })
  cat(sprintf("%s: logLik %.4g, conv %d, betad %s\n", fam$family,
              as.numeric(logLik(f)), f$opt$convergence,
              paste(format(f$estimates$betad, digits = 4), collapse = " ")))
  for (w in ws) cat("  ", substr(gsub("\\s+", " ", w), 1, 300), "\n")
}
cat("DONE\n")
