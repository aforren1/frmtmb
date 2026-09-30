# One-off edit of R/sampling-api.R: export response_codes_newdata for
# frmtmb.sample. Kept as the record of the edit.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/R/sampling-api.R"
x <- readLines(p)
i <- which(x == "#' @aliases arma_cond_fill_dpars")
stopifnot(length(i) == 1L)
x <- append(x, "#' @aliases response_codes_newdata", i)
j <- which(x == "#'   arma_cond_resp, arma_cond_dpars, arma_cond_fill_dpars,")
stopifnot(length(j) == 1L)
x[j + 1L] <- "#'   response_codes_newdata, subset_resp_check, subset_newdata)"
k <- which(x == "#' response it is `dpars_fn(fit)`.")
stopifnot(length(k) == 1L)
x <- append(x, c("#'",
  "#' `response_codes_newdata(rspec, y, what)` codes a response read from",
  "#' newdata as the fit coded its own: a bernoulli response takes the 0",
  "#' and 1 of the fit's two values (brms's level order), and is refused,",
  "#' naming `what`, when it holds a third value. Any other response is",
  "#' returned unchanged."), k)
con <- file(p, "wb"); writeLines(x, con, sep = "\r\n"); close(con)
