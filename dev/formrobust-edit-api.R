# One-off edit of R/sampling-api.R: export arma_cond_fill_dpars for
# frmtmb.sample. Kept as the record of the edit.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/R/sampling-api.R"
x <- readLines(p)
i <- which(x == "#' @aliases arma_cond_dpars")
stopifnot(length(i) == 1L)
x <- append(x, "#' @aliases arma_cond_fill_dpars", i)
j <- which(x == "#'   arma_cond_resp, arma_cond_dpars, subset_resp_check, subset_newdata)")
stopifnot(length(j) == 1L)
x[j] <- "#'   arma_cond_resp, arma_cond_dpars, arma_cond_fill_dpars,"
x <- append(x, "#'   subset_resp_check, subset_newdata)", j)
con <- file(p, "wb"); writeLines(x, con, sep = "\r\n"); close(con)
