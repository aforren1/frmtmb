f <- "R/compat.R"
x <- readLines(f)
old <- "      \"Refused under both forms. brms refuses se() with cov = FALSE ('Please set cov = TRUE in ARMA structures when including known standard errors'), and the cov = TRUE form here has no per-row density for a known standard error to add to.\")"
new <- "      \"Refused under both forms. brms refuses se() with cov = FALSE ('Please set cov = TRUE in ARMA structures when including known standard errors'). Under cov = TRUE brms 2.23.0 fits it, adding se^2 to the diagonal of each group's covariance (normal_time_hom_se_lpdf() and its relatives); frmtmb refuses it there, a parity gap.\")"
i <- which(x == old)
stopifnot(length(i) == 1L)
x[i] <- new
writeLines(x, f)
