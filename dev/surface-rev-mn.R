# Reviewer: method = "predict" on multinomial (claim 6 message scope)
source("dev/surface-rev-env.R")
arm <- rev_env(commandArgs(TRUE)[1])
suppressPackageStartupMessages(library(frmtmb))
grDevices::pdf(NULL)
set.seed(3); n <- 200
d <- data.frame(x = rnorm(n))
pr <- cbind(1, exp(0.5 + d$x), exp(-0.3 - d$x)); pr <- pr / rowSums(pr)
d$Y <- t(sapply(seq_len(n), function(i) rmultinom(1, 10, pr[i, ])))
colnames(d$Y) <- c("a", "b", "c"); d$size <- 10
f <- frm(bf(Y | trials(size) ~ x), family = multinomial(K = 3), data = d)
cat("dpars:", names(f$spec$responses[[1]]$dpars), "type:", f$spec$responses[[1]]$family$type, "\n")
rev_show("multinomial CE predict", conditional_effects(f, "x", method = "predict"))
rev_show("multinomial CE epred", conditional_effects(f, "x"))
