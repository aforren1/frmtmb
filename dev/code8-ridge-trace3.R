# Part C: what se_ridge_ritz() returns on a draw (lane build).
#   CODE8_ARM=after Rscript dev/code8-ridge-trace3.R <j...>
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
trace("se_ridge_ritz", tracer = quote(cat("  [ritz in] S:", S, " ev:",
  format(e$values, digits = 3), "\n")), exit = quote({
  rv <- returnValue()
  cat("  [ritz] values:", format(rv$values, digits = 3), " noise:",
      format(rv$noise[1], digits = 3), "\n")
  print(round(rv$vectors, 6))
}), where = ns, print = FALSE)
trace("se_tier3", exit = quote({
  rv <- returnValue()
  cat("  [se_tier3] lost:", paste(names(rv$lost), rv$lost), "\n")
}), where = ns, print = FALSE)
for (j in js) {
  dn$yn <- y0 * (1 + 2e-15 * j)
  cat("== draw", j, "\n")
  r <- rfit(form, data = dn)
  cat(" conds:", substr(cond_txt(r), 1, 150), "\n")
}
cat("DONE\n")
