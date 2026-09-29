# Reviewer: is frm_bootstrap() really untouched? Its default is
# re_formula = NA with redraw_smooths = TRUE, which takes a different
# branch of sim_re_plan(). Saved for identical() across the arms.
ARM <- if (identical(Sys.getenv("REVLIB"), "base")) "base" else "lane"
LIB <- if (ARM == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-resmooth-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("ARM:", ARM, "| frmtmb:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")

set.seed(5)
n <- 300L
d <- data.frame(x = stats::runif(n), g = factor(rep(1:10, each = 30L)))
d$y <- sin(2 * pi * d$x) + stats::rnorm(10, 0, 0.5)[d$g] * d$x +
  stats::rnorm(n, 0, 0.3)
out <- list()
for (lab in c("fs", "sx", "sx_g")) {
  fm <- switch(lab,
    fs = bf(y ~ s(x, g, bs = "fs", k = 5)),
    sx = bf(y ~ s(x)),
    sx_g = bf(y ~ s(x) + (1 | g)))
  fit <- suppressWarnings(frm(fm, family = gaussian(), data = d))
  pl <- get("sim_re_plan", envir = ns)
  a <- pl(fit, NA, smooths = TRUE)
  b <- pl(fit, NA)
  cat(sprintf("%-5s sim_re_plan(NA, smooths=TRUE): blocks %-8s every=%-5s ",
              lab, paste(a$blocks, collapse = ","), a$every))
  cat(sprintf("| sim_re_plan(NA): blocks %-8s every=%-5s | class(redraw) %s\n",
              paste(b$blocks, collapse = ","), b$every, class(b$redraw)))
  bt <- suppressWarnings(frm_bootstrap(fit, nsim = 20, seed = 9,
                                       FUN = function(f) {
                                         as.numeric(fixef(f)[, "Estimate"])
                                       }))
  out[[lab]] <- list(t0 = bt$t0, t = bt$t, converged = bt$converged)
  cat(sprintf("      frm_bootstrap t0 %s | sd(t[,1]) %.6f | converged %d\n",
              paste(sprintf("%.5f", bt$t0), collapse = " "),
              stats::sd(bt$t[, 1L]), sum(bt$converged)))
}
saveRDS(out, file.path("dev", paste0("resmooth-rev-boot-", ARM, ".rds")))
cat("DONE\n")
