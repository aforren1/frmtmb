# Reviewer, lane sampfix, script 09: the lane's frmtmb.sample on the
# 0.65.0 core (rellib-r3), which declares no inverse threshold map. The
# NEWS says the draws then keep the internal names. Seed 405 data,
# sampler seed 3.
LANE <- "C:/Users/adf44/source/r/wt-sampfix-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
.libPaths(c(REF, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(frmtmb, lib.loc = REF))
suppressMessages(library(frmtmb.sample, lib.loc = LANE))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
cat("frmtmb.sample from", dirname(find.package("frmtmb.sample")), "\n")
cat("core declares ord_thresholds_raw:",
    !is.null(frmtmb::cumulative()$post$ord_thresholds_raw), "\n")
set.seed(405L); n <- 300L
do <- data.frame(x = rnorm(n), fc = factor(sample(c("a", "b", "c"), n, TRUE)))
do$yo <- factor(cut(0.8 * do$x + rlogis(n), c(-Inf, -0.5, 0.7, Inf),
                    labels = FALSE), ordered = TRUE)
for (spec in list(list(bf(yo ~ x), cumulative()), list(bf(yo ~ x + cs(fc)), sratio()))) {
  fam <- spec[[2]]
  fit <- frm(spec[[1]], family = fam, data = do)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 3)))
  cat(fam$family, ": ", paste(colnames(ds$draws), collapse = " "), "\n", sep = "")
  cat("  fixef rows: ", paste(rownames(fixef(ds)), collapse = ","),
      "  epred finite: ", all(is.finite(posterior_epred(ds, ndraws = 5))), "\n", sep = "")
}
