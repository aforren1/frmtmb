# Lane parity2: two follow-ups to dev/parity2-brms-probe2.R on brms
# 2.23.0. (1) the draw-margin names of an ordinal posterior_epred() under
# ndraws and draw_ids; (2) where err sits in variables() of a latent
# ar(cov = FALSE) fit. fixed_param, 20 draws, seed 1.
#   Rscript dev/parity2-brms-probe3.R > dev/parity2-out/brms-probe3.txt 2>&1
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
options(mc.cores = 1)
here <- "C:/Users/adf44/source/r/frmtmb-wt-parity2/dev"
source(file.path(here, "parity2-data.R"))
bfix <- function(formula, data, ...) {
  suppressMessages(suppressWarnings(
    brm(formula, data = data, algorithm = "fixed_param", chains = 1,
        iter = 20, warmup = 0, init = 0, refresh = 0, seed = 1,
        silent = 2, ...)))
}
fo <- bfix(yo ~ x, dm, family = cumulative())
cat("ndraws = 5:", dimnames(posterior_epred(fo, ndraws = 5))[[1]], "\n")
cat("draw_ids = c(3, 7, 12):",
    dimnames(posterior_epred(fo, draw_ids = c(3, 7, 12)))[[1]], "\n")
dp <- data.frame(t = rep(1:8, 6), g = factor(rep(1:6, each = 8)),
                 x = dg$x)
set.seed(4)
dp$y <- rpois(48, exp(0.5 + 0.3 * dp$x))
fa <- bfix(y ~ x + ar(time = t, gr = g, cov = FALSE), dp,
           family = poisson())
v <- variables(fa)
cat("variables:", paste(v[!grepl("^err\\[([2-9]|[1-9][0-9])", v)],
                        collapse = " "), "\n")
cat("err positions", range(grep("^err\\[", v)), "of", length(v), "\n")
cat("done\n")
