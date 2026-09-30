# Punch round 1, m3: export arma_cond_fill_epred() for frmtmb.sample.
# Record of the edit of R/sampling-api.R.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/R/sampling-api.R"
x <- readLines(p)
i <- which(x == "#' @aliases arma_cond_fill_dpars"); stopifnot(length(i) == 1L)
x <- append(x, "#' @aliases arma_cond_fill_epred", i)
j <- which(x == "#'   arma_cond_resp, arma_cond_dpars, arma_cond_fill_dpars,")
stopifnot(length(j) == 1L)
x[j] <- "#'   arma_cond_resp, arma_cond_dpars, arma_cond_fill_dpars,"
x[j + 1L] <- "#'   arma_cond_fill_epred, response_codes_newdata, subset_resp_check,"
x <- append(x, "#'   subset_newdata)", j + 1L)
k <- which(x == "#' response it is `dpars_fn(fit)`."); stopifnot(length(k) == 1L)
x <- append(x, c("#'",
  "#' `arma_cond_fill_epred(fit, rspec, newdata, re_formula, dpar)` is the",
  "#' expected response of one draw on `newdata` with brms's fill: each",
  "#' missing `cov = FALSE` response is a draw at its shifted mean, so",
  "#' `posterior_epred()` carries the spread of the unobserved past, as",
  "#' brms's does. It returns `NULL` when nothing needs filling, and the",
  "#' caller keeps its own route."), k)
con <- file(p, "wb"); writeLines(x, con, sep = "\r\n"); close(con)
