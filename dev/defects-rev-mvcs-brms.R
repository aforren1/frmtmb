# Reviewer of lane defects, recheck of B2: brms 2.23.0's
# fitted(scale = "linear") on multivariate fits with ordinal cs()
# responses. Checks brms's own layer identity per draw and saves the
# posterior means for dev/defects-rev-mvcs-frmtmb.R.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
source("dev/defects-rev-mvcs-data.R")
d <- mvcs_data()
out <- list()
for (nm in names(mvcs_models)) {
  cat("=====", nm, "\n")
  f <- suppressWarnings(brm(eval(mvcs_models[[nm]]), data = d, chains = 2,
                            iter = 2000, refresh = 0, seed = 1))
  a <- fitted(f, scale = "linear")
  cat("dim:", dim(a), " layers:", dimnames(a)[[3]], "\n")
  pl <- posterior_linpred(f)
  dr <- as_draws_df(f)
  if (nm == "one_cs") {
    # brms's own identity per draw: eta_k = b_x x + bcs_k z
    h1 <- outer(dr$b_yo_x, d$x) + outer(dr$`bcs_yo_z[1]`, d$z)
    h2 <- outer(dr$b_yo_x, d$x) + outer(dr$`bcs_yo_z[2]`, d$z)
    cat("posterior_linpred dims:", dim(pl), "\n")
    lin <- posterior_linpred(f, resp = "yo")
    cat("max |brms eta1 - (b x + bcs1 z)| over draws:", max(abs(lin[, , 1] - h1)),
        "  eta2:", max(abs(lin[, , 2] - h2)), "\n")
  }
  out[[nm]] <- list(fitted = a, means = colMeans(as.matrix(dr[, grep("^b_|^bcs_|^Intercept", names(dr))])))
  print(round(out[[nm]]$means, 3))
}
saveRDS(out, "dev/defects-rev-log/mvcs-brms.rds")
