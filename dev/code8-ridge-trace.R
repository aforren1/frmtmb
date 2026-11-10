# Part C: which step of se_tier3() decides each draw's lost set.
#   CODE8_ARM=before|after Rscript dev/code8-ridge-trace.R <j...>
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
for (fn in c("se_line_probe", "se_curvature_real", "se_span_remeasure",
             "se_ridge_curvature", "se_obj_curv", "conv_headroom",
             "grad_headroom", "se_boundary_names", "se_edge_sd")) {
  eval(bquote(trace(.(fn), exit = quote({
    rv <- returnValue()
    cat("  [", .(fn), "] ->", if (is.list(rv)) paste(names(rv), collapse = ",")
        else format(rv, digits = 4), "\n")
  }), where = ns, print = FALSE)))
}
trace("se_tier3", exit = quote({
  rv <- returnValue()
  cat("  [se_tier3] lost:", paste(names(rv$lost), rv$lost), " ev:",
      format(rv$ev, digits = 3), "\n")
}), where = ns, print = FALSE)
for (j in js) {
  dn$yn <- y0 * (1 + 2e-15 * j)
  cat("== draw", j, "\n")
  r <- rfit(form, data = dn)
  cat(" fit done; conds:", substr(cond_txt(r), 1, 200), "\n")
  sdr <- quiet(fns$sdr_of(r$f))
  cat(" final lost:", paste(names(sdr$se_lost), sdr$se_lost), "\n")
}
cat("DONE\n")
