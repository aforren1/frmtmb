# Lane fixes, item 1: emmeans() on a transformed predictor against
# emmeans() on glm() / lm() fits of the same model (same MLE).
#   Rscript dev/fixes-emm2.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(emmeans)})
cat("LIB", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
set.seed(20260930)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n, 1, 2), time = runif(n, 1, 5),
                f = factor(sample(c("a", "b", "c"), n, TRUE)))
d$yc <- rpois(n, d$time * exp(0.3 + 0.4 * d$x + 0.2 * d$z))
d$yg <- 1 + 0.5 * d$z - 0.1 * d$z^2 + as.numeric(d$f) + rnorm(n)
cmp <- function(lab, fo, family, spec, ...) {
  r <- tryCatch({
    fit <- frm(bf(fo), data = d, family = family)
    ref <- if (family$family == "gaussian") lm(fo, data = d) else
      glm(fo, data = d, family = family)
    a <- summary(emmeans(fit, spec, ...))
    b <- summary(emmeans(ref, spec, ...))
    ra <- emmeans::ref_grid(fit, ...)@grid
    rb <- emmeans::ref_grid(ref, ...)@grid
    sprintf(paste("grid identical: %s | emmean max rel diff %.3g |",
                  "SE max rel diff %.3g | %s"),
            isTRUE(all.equal(ra[setdiff(names(ra), ".wgt.")],
                             rb[setdiff(names(rb), ".wgt.")])),
            max(abs(a$emmean - b$emmean) / abs(b$emmean)),
            max(abs(a$SE - b$SE) / b$SE),
            paste(format(a$emmean, digits = 7), collapse = " "))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-44s %s\n", lab, r))
}
pois <- poisson()
gaus <- gaussian()
cmp("pois poly(z,2)+f", yc ~ poly(z, 2) + f, pois, "f")
cmp("pois poly(z,2)+f+offset", yc ~ poly(z, 2) + f + offset(log(time)),
    pois, "f")
cmp("pois log(abs(z)+1)+f", yc ~ log(abs(z) + 1) + f, pois, "f")
cmp("pois scale(z)+f", yc ~ scale(z) + f, pois, "f")
cmp("pois scale(z)+f+offset", yc ~ scale(z) + f + offset(log(time)),
    pois, "f")
cmp("pois scale(z)+f at z=3", yc ~ scale(z) + f, pois, "f",
    at = list(z = 3))
cmp("pois poly(z,2)+f at z=c(-1,2)", yc ~ poly(z, 2) + f, pois, c("z", "f"),
    at = list(z = c(-1, 2)))
cmp("pois scale(z)*f by z at z=c(0,4)", yc ~ scale(z) * f, pois, "f",
    at = list(z = c(0, 4)))
cmp("gaus poly(z,2)+f", yg ~ poly(z, 2) + f, gaus, "f")
cmp("gaus scale(z)+f at z=-2", yg ~ scale(z) + f, gaus, "f",
    at = list(z = -2))
cmp("gaus log(abs(z)+1)*f at z=c(0,3)", yg ~ log(abs(z) + 1) * f, gaus,
    c("z", "f"), at = list(z = c(0, 3)))
# type = "response" on the poisson log link
cmp("pois poly(z,2)+f response", yc ~ poly(z, 2) + f, pois, "f",
    type = "response")
