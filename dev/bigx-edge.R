# Lane bigx: test-se-reparam.R's reparam_edge_fit() designs (raw powers
# of x on [1, 2], (1 | g), no group effect in the data) with the
# optimizer in the design's coordinates and in the natural ones, against
# lme4's ML maximum.
#   BIGX_ARM=before|after Rscript dev/bigx-edge.R
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/bigx-common.R")
suppressPackageStartupMessages(library(lme4))
edge <- function(deg, seed) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(300, 1, 2),
                  g = factor(rep(1:20, length.out = 300)))
  d$y <- sin(3 * d$x) + stats::rnorm(300, 0, 0.3)
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                                "+ (1 | g)"))
  list(d = d, fo = fo)
}
for (cs in list(c(5, 4), c(7, 4), c(6, 4), c(8, 1))) {
  e <- edge(cs[1], cs[2])
  l4 <- lmer(e$fo, data = e$d, REML = FALSE,
             control = lmerControl(check.conv.singular = "ignore"))
  arms <- if (bigx_arm == "after") c(TRUE, FALSE) else TRUE
  for (b in arms) {
    old <- options(frmtmb.opt_basis = b)
    x <- bigx_fit(e$fo, data = e$d)
    options(old)
    f <- x$f
    lost <- tryCatch(names(frmtmb:::sdr_of(f)$se_lost),
                     error = function(e) "?")
    cat(sprintf(paste0("deg %d seed %d basis %-5s ll - lmer %+.3e ",
                       "log sd %.3f code %d sd_escape %s lost [%s] ",
                       "warn %d %s\n"),
                cs[1], cs[2], b, as.numeric(logLik(f)) -
                  as.numeric(logLik(l4)),
                f$opt$par[["theta"]], f$opt$convergence,
                if (is.null(f$opt$sd_escape)) "-" else
                  format(f$opt$sd_escape[["gain"]], digits = 3),
                paste(lost, collapse = ","), length(x$w),
                if (length(x$w)) substr(x$w[1], 1, 60) else ""))
  }
  cat(sprintf("   lmer sd %.3g\n", attr(VarCorr(l4)$g, "stddev")))
}
