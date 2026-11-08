# Review of lane mfxse, attack 8: does interop_vcov_re() work on lane
# latperf's sparse joint covariance (frm_joint_sparse)?
#   Rscript dev/review-mfxse-latperf.R
# Reads wt-latperf-lib (lane latperf's build, read-only) and evaluates this
# lane's interop functions from R/interop.R in a child of its namespace. Data
# seed 9501: poisson (1 | g) + ar(cov = TRUE), 400 rows.
.libPaths(c("C:/Users/adf44/source/r/wt-latperf-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
cat("latperf's get_joint_cov formals:", names(formals(ns$get_joint_cov)), "\n")
ex <- parse("C:/Users/adf44/source/r/frmtmb-wt-mfxse/R/interop.R", keep.source = FALSE)
want <- c("interop_re_part", "interop_has_b", "block_reads_theta",
          "interop_re_names", "interop_joint_map", "interop_vcov_re",
          "mo_gap_cum_jacobian", "interop_vcov", "interop_coef_vector",
          "interop_coef_names", "interop_extra_comps", "interop_mo_gaps",
          "interop_mo_const")
env <- new.env(parent = ns)
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]]) %in% want) eval(e, env)
}
cat("loaded:", sort(ls(env)), "\n")
set.seed(9501)
J <- 20; m <- 20
d <- data.frame(g = factor(rep(1:J, each = m)), t = rep(1:m, J), x = rnorm(J * m))
d$y <- rpois(J * m, exp(0.3 + 0.4 * d$x + rnorm(J, 0, 0.6)[d$g]))
fit <- suppressWarnings(frm(y ~ x + (1 | g) + ar(time = t, gr = g, cov = TRUE),
                            family = poisson(), data = d))
jc <- ns$get_joint_cov(fit)
cat("joint covariance class:", class(jc$V), "\n")
V0 <- env$interop_vcov(fit)
t1 <- system.time(Vs <- tryCatch(env$interop_vcov_re(fit, V0),
                                 error = function(e) conditionMessage(e)))
if (is.character(Vs)) cat("sparse route ERROR:", Vs, "\n") else {
  cat("sparse route: dim", dim(Vs), "held:", attr(Vs, "held") %||% "none",
      "time", t1[["elapsed"]], "s\n")
  # the same with the joint covariance made dense
  fit$cache$Vjoint <- list(V = as.matrix(jc$V), names = jc$names)
  Vd <- env$interop_vcov_re(fit, V0)
  cat("max |sparse - dense| / max|dense|:", format(max(abs(Vs - Vd)) / max(abs(Vd)), digits = 3), "\n")
}
