# Lane seflags, item 4: the degree-7 raw polynomial with (1 | g), seed 7
# (dev/setier-rev-deg7.R): the spectrum tier 3 reads, each direction's
# own noise, the gap and tau, and each parameter's projection onto the
# removed subspace and onto the subspace of the smallest m directions.
#   Rscript dev/seflags-deg7diag.R <lib> [deg] [seed]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
ns <- asNamespace("frmtmb")
deg <- if (length(a) > 1) as.integer(a[2]) else 7L
sd <- if (length(a) > 2) as.integer(a[3]) else deg
set.seed(sd)
d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
rhs <- paste0("x", 1:deg, collapse = " + ")
fo <- as.formula(paste("y ~", rhs, "+ (1 | g)"))
f <- suppressWarnings(frm(fo, data = d))
m <- suppressWarnings(lmer(fo, data = d, REML = FALSE))
cat("logLik frmtmb", format(-f$opt$objective, digits = 12), "lmer",
    format(as.numeric(logLik(m)), digits = 12), "\n")
h <- f$cache$hessian_fixed %||% ns$fit_outer_hessian(f)
H <- h$H; E <- h$E
D <- sqrt(abs(diag(H)))
e <- eigen(H / outer(D, D), symmetric = TRUE)
Es <- E / outer(D, D)
nk <- sqrt(colSums((Es %*% e$vectors)^2))
big <- max(abs(e$values))
cat("ev:", signif(e$values, 3), "\n")
cat("nk:", signif(nk, 3), "\n")
cat("thr:", signif(pmax(1e-9 * big, 3 * nk), 3), "\n")
cat("enorm", signif(sqrt(sum(Es^2)), 3), "\n")
nm <- ns$outer_par_names(f)
np <- length(nm)
for (mm in 1:4) {
  idx <- (np - mm + 1):np
  gap <- e$values[np - mm] - e$values[np - mm + 1]
  tau <- 10 * sqrt(sum(Es^2)) / gap
  pr <- sqrt(rowSums(e$vectors[, idx, drop = FALSE]^2))
  cat(sprintf("smallest %d: gap %.3g tau %.3g | proj %s\n", mm, gap, tau,
              paste(nm, signif(pr, 2), sep = "=", collapse = " ")))
}
sl <- sqrt(diag(as.matrix(vcov(m))))
se <- suppressWarnings(sqrt(diag(vcov(f))))
cat("frmtmb SE:", signif(se[seq_along(sl)], 4), "\n")
cat("lmer SE:  ", signif(sl, 4), "\n")
print(suppressWarnings(ns$sdr_of(f))$se_lost)
