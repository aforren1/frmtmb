# Lane predraw2: the excess distribution of S0.1 seed 14, whose 20-draw
# check medians exceed 1 in 12% of random subsets though its 1000-draw
# median is 0.11 (dev/predraw2-ncheck.R). Usage: Rscript <this> <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
source("dev/predraw2-cases.R")
fit <- case_fit("S0.1:14")$fit
set.seed(7)
dr <- fq(ns$predict_par_drawer(fit, TRUE, 2L))
ds <- ns$fit_draw_space(fit)
L <- chol(ds$V)
law <- fit$cache[["draw_edge_law"]]
if (!is.list(law)) law <- NULL
cat("law:", if (is.null(law)) "normal" else
  paste("edge", length(law$edge), "narrowed", length(law$narrowed)), "\n")
v0 <- ns$fit_outer_vector(fit, ds$map)
P <- length(v0)
print(round(cbind(est = v0, sd = sqrt(diag(ds$V))), 3))
dw <- list(L = L, law = law, v0 = v0, held = integer(0), off = 0L)
set.seed(7)
Z <- matrix(stats::rnorm(P * 1000), P, 1000)
D <- v0 + if (is.null(law)) crossprod(L, Z) else ns$draw_edge_apply(law, L, Z)
r <- ns$draw_logratio(fit, dw, D, Z)
cat("e quantiles .1 .25 .5 .6 .7 .75 .9:",
    format(stats::quantile(r$e, c(.1, .25, .5, .6, .7, .75, .9)), digits = 3),
    "\n")
j <- which(names(v0) == "theta_2" | seq_along(v0) == P - 1L)
cat("share e > 1:", mean(r$e > 1), "\n")
for (k in seq_len(P)) {
  cat(sprintf("param %d: corr(e > 1, draw > estimate) %.3f\n", k,
              stats::cor(r$e > 1, D[k, ] > v0[k])))
}
