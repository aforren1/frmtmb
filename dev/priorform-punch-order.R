# How the more-specific-wins rule changes frm() fits: the same four
# priors fitted on the base build (later wins) and the lane build (more
# specific wins), in both written orders.
#   Rscript dev/priorform-punch-order.R ref|lane   (writes an .rds)
#   Rscript dev/priorform-punch-order.R compare
args <- commandArgs(trailingOnly = TRUE)
which_lib <- args[1]
out_rds <- function(w) paste0("dev/priorform-punch-order-", w, ".rds")
if (identical(which_lib, "compare")) {
  a <- readRDS(out_rds("ref"))
  b <- readRDS(out_rds("lane"))
  cat("<!-- priorform-punch-order:begin -->\n")
  cat("| case | base logLik | lane logLik | fit identical |\n")
  cat("|---|---|---|---|\n")
  for (nm in names(a)) {
    cat(sprintf("| %s | %.6f | %.6f | %s |\n", nm, a[[nm]]$ll, b[[nm]]$ll,
                identical(a[[nm]], b[[nm]])))
  }
  cat("<!-- priorform-punch-order:end -->\n")
  quit(save = "no")
}
lib <- switch(which_lib,
  ref = "C:/Users/adf44/source/r/rellib-r3",
  lane = "C:/Users/adf44/source/r/priorform-lib")
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(20260916)
n <- 180
d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                g = factor(rep(1:15, 12)))
d$yb <- stats::plogis(0.3 * d$x - 0.5 + stats::rnorm(15, 0, 0.5)[d$g] +
                        stats::rnorm(n, 0, 0.3))
d$y <- 1 + 0.5 * d$x + stats::rnorm(15, 0, 0.5)[d$g] + stats::rnorm(n)
coef_x <- set_prior("normal(3, 0.01)", class = "b", coef = "x")
cls_b <- set_prior("normal(0, 0.01)", class = "b")
grp_sd <- set_prior("normal(0, 5)", class = "sd", group = "g")
cls_sd <- set_prior("normal(0, 0.05)", class = "sd")
beta_f <- bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta()
cases <- list(
  b_class_then_coef = list(y ~ x + z + (1 | g), cls_b + coef_x),
  b_coef_then_class = list(y ~ x + z + (1 | g), coef_x + cls_b),
  sd_class_then_group = list(beta_f, cls_sd + grp_sd),
  sd_group_then_class = list(beta_f, grp_sd + cls_sd))
res <- lapply(cases, function(cs) {
  f <- suppressWarnings(frm(cs[[1]], data = d, prior = cs[[2]]))
  list(ll = as.numeric(logLik(f)), par = f$opt$par)
})
saveRDS(res, out_rds(which_lib))
cat("saved", out_rds(which_lib), "\n")
