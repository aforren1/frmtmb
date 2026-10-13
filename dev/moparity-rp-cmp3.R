# The lane's centered chain of dev/moparity-rp-cmp.R with the per-scale
# floor taken out (sd_floored() bound to exp() in the registry): if its
# draws are base's, the chains differ only where the floor acts.
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
ns <- asNamespace("frmtmb")
reg <- get("covstruct_registry", ns)
old_env <- new.env(parent = ns)
old_env$sd_floored <- exp
for (nm in names(reg)) {
  f <- reg[[nm]]$nll
  if (is.function(f) && identical(environment(f), ns)) {
    environment(f) <- old_env
    reg[[nm]]$nll <- f
  }
}
unlockBinding("covstruct_registry", ns)
assign("covstruct_registry", reg, envir = ns)
lockBinding("covstruct_registry", ns)
set.seed(31)
ng <- 12L
dd <- data.frame(x = stats::rnorm(120),
                 g = factor(rep(seq_len(ng), length.out = 120)))
u <- matrix(stats::rnorm(2 * ng), 2)
dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x + 0.7 * u[1, dd$g] +
                       0.4 * u[2, dd$g] * dd$x, 1)
cen <- suppressWarnings(suppressMessages(
  frm_sample(bf(y ~ x + (x | g)), family = gaussian(), data = dd,
             chains = 1, iter = 1000, refresh = 0, seed = 33,
             reparameterize = FALSE)))
a <- readRDS("dev/moparity-log/rp-cmp-base.rds")
b <- readRDS("dev/moparity-log/rp-cmp-lane.rds")
cat("lane without the floor vs base: identical", identical(cen$draws, a$cen),
    "| vs lane with it: identical", identical(cen$draws, b$cen), "\n")
th <- grep("^theta", colnames(cen$draws))
cat("smallest saved theta: base", min(a$cen[, th]), "lane", min(b$cen[, th]),
    "\n")
