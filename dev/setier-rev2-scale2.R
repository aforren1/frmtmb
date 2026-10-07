# Reviewer of lane setier, re-check: why the boundary verdict changes with
# the units of the response (dev/setier-rev2-scale.R). Per seed and
# scale on the ri20 design: logLik against lme4 at the same scale, the
# group sd relative to sigma, the theta curvature, se_at_edge(), the
# verdict.
#   Rscript dev/setier-rev2-scale2.R <lib> [seeds]
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:40
suppressMessages({library(frmtmb); library(lme4)})
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), "\n")
X <- list()
for (s in seeds) for (sc in c(1e-3, 1, 1e3)) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- (1 + 0.5 * d$x + rnorm(100)) * sc
  m <- character()
  f <- withCallingHandlers(frm(y ~ x + (1 | g), data = d),
    warning = function(x) invokeRestart("muffleWarning"),
    message = function(x) {m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")})
  l4 <- suppressMessages(lmer(y ~ x + (1 | g), data = d, REML = FALSE))
  nm <- ns$outer_par_names(f)
  j <- match("theta_1", nm)
  hc <- f$cache$hessian_fixed
  sig <- exp(f$opt$par[match("sigma_(Intercept)", nm)])
  X[[length(X) + 1L]] <- data.frame(
    seed = s, scale = sc, code = f$opt$convergence,
    dll = as.numeric(logLik(f)) - as.numeric(logLik(l4)),
    rel_sd = exp(f$opt$par[j]) / sig,
    lme4_rel_sd = attr(VarCorr(l4)$g, "stddev") / sigma(l4),
    hd = if (!is.null(hc)) abs(hc$H[j, j]) else NA,
    edge = if (exists("se_at_edge", ns)) ns$se_at_edge(f, "theta_1", -1)
      else NA,
    flagged = any(grepl("^Boundary", m)), sing = isSingular(l4))
}
X <- do.call(rbind, X)
write.table(X, file.path("dev/setier-rev2-log", paste0("scale2-", basename(dirname(find.package("frmtmb"))), ".tsv")), sep = "\t", quote = FALSE,
            row.names = FALSE)
w <- reshape(X[, c("seed", "scale", "flagged")], idvar = "seed",
             timevar = "scale", direction = "wide")
diffs <- w$seed[apply(w[, -1], 1, function(r) length(unique(r)) > 1)]
cat("seeds whose verdict differs:", length(diffs), "\n")
cat("logLik - lme4, all fits: min", signif(min(X$dll), 3), "max",
    signif(max(X$dll), 3), "\n")
for (sc in unique(X$scale)) {
  x <- X[X$scale == sc, ]
  cat(sprintf(paste0("scale %-6g | below lme4 by > 1e-4: %d | lme4 ",
                     "singular %d, flagged %d (singular %d, other %d) | ",
                     "median rel sd of the singular fits %.2g\n"),
              sc, sum(x$dll < -1e-4), sum(x$sing), sum(x$flagged),
              sum(x$flagged & x$sing), sum(x$flagged & !x$sing),
              median(x$rel_sd[x$sing])))
}
print(X[X$seed %in% head(diffs, 6), ], row.names = FALSE, digits = 3)
