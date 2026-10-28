# Lane reslot: a registered covariance structure (frmtmb.spatial's car())
# keeps class theta for the parameter only it reaches. Usage:
# Rscript dev/reslot-car.R before|after. Seed 7, 6 x 6 lattice.
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.spatial))
cat("frmtmb.spatial", format(packageVersion("frmtmb.spatial")),
    dirname(find.package("frmtmb.spatial")), "\n")
set.seed(7)
k <- 6
W <- matrix(0, k * k, k * k)
for (i in seq_len(k)) for (j in seq_len(k)) {
  a <- (i - 1) * k + j
  if (j < k) W[a, a + 1] <- W[a + 1, a] <- 1
  if (i < k) W[a, a + k] <- W[a + k, a] <- 1
}
rownames(W) <- colnames(W) <- paste0("s", seq_len(k * k))
loc <- factor(rep(rownames(W), each = 4), levels = rownames(W))
d <- data.frame(loc = loc, x = rnorm(length(loc)))
d$y <- 1 + 0.5 * d$x + rnorm(k * k)[as.integer(loc)] + rnorm(nrow(d))
f <- bf(y ~ x + car(W, gr = loc, type = "escar"))
tab <- default_prior(f, data = d, data2 = list(W = W))
print(as.data.frame(tab)[tab$class %in% c("sd", "theta"),
                         c("class", "coef", "group")], row.names = FALSE)
des <- frmtmb:::prior_design(f, d, NULL, list(W = W))
for (pl in list(set_prior("normal(0, 1)", class = "theta", coef = "theta_2"),
                set_prior("normal(0, 1)", class = "theta"),
                set_prior("normal(0, 1)", class = "theta", coef = "theta_1"),
                set_prior("exponential(1)", class = "sd", group = "loc"))) {
  r <- tryCatch(frmtmb:::resolve_priorlist(des, pl), error = function(e) e)
  cat("--", format(pl$class), format(pl$coef), ":",
      if (inherits(r, "error")) paste("ERROR", conditionMessage(r)) else
        paste(vapply(r$entries, function(e) paste(e$comp, e$idx), ""),
              collapse = ", "), "\n")
}
