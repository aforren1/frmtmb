# Part C: the Hessian each se_tier3() call of one draw reads.
#   CODE8_ARM=before|after Rscript dev/code8-ridge-trace2.R <j>
js <- as.integer(commandArgs(TRUE))
source("dev/code8-common.R")
options(width = 200)
set.seed(77)
dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(8, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
y0 <- dn$yn
form <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
ns <- asNamespace("frmtmb")
trace("se_tier3", tracer = quote({
  cat("  [se_tier3 in] p:", format(p, digits = 5), "\n   diagH:",
      format(diag(H), digits = 4), "\n   maxE:",
      format(apply(E, 1, max), digits = 3), " units:",
      format(fit$par_units %||% 1, digits = 3), " exact:", exact,
      " held:", held, " edge:", edge, "\n")
}), exit = quote({
  rv <- returnValue()
  cat("  [se_tier3 out] lost:", paste(names(rv$lost), rv$lost), "\n")
}), where = ns, print = FALSE)
for (j in js) {
  dn$yn <- y0 * (1 + 2e-15 * j)
  cat("== draw", j, "\n")
  r <- rfit(form, data = dn)
  cat(" final p:", format(r$f$opt$par, digits = 5), "\n")
}
cat("DONE\n")
