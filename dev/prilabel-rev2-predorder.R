# Re-check of prilabel fix round 1, finding 4: the order ACROSS
# predictors (dpars written out of brms's order, nlpars), brms's
# frame_re ids against the lane's blocks. Data seed 71.
.libPaths(c("C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
})
sink(file.path(wt, "dev/prilabel-rev-log/rev2-predorder.txt"), split = TRUE)
set.seed(71)
n <- 300
d <- data.frame(x = runif(n), g = factor(sample(1:8, n, TRUE)),
                h = factor(sample(1:6, n, TRUE)),
                ga = factor(sample(1:5, n, TRUE)))
d$y <- 2 * exp(0.5 * d$x) + rnorm(8, 0, 0.3)[d$g] + rnorm(n, 0, 0.3)
d$cnt <- rpois(n, 2)
try_ <- function(expr) tryCatch(expr, error = function(e) {
  paste("ERROR:", conditionMessage(e))
})
brms_ids <- function(f, fam) {
  re <- as.data.frame(brms:::frame_re(brms:::brmsterms(
    brms::brmsformula(f[[1]], flist = f[-1], family = fam,
                      nl = isTRUE(attr(f, "nl")))), data = d))
  paste(re$dpar, re$nlpar, re$group, re$coef, sep = "|")
}
frm_ids <- function(f, fam) {
  bff <- do.call(frmtmb::bf, c(f, list(nl = isTRUE(attr(f, "nl")))))
  des <- frmtmb:::prior_design(bff, d, fam, list(), check_response = FALSE)
  out <- character(0)
  for (bk in des$frame$re_blocks) {
    cp <- frmtmb:::block_col_prefix(des$spec, des$frame, bk)
    out <- c(out, paste(cp$dpar, cp$nlpar, bk$group_name,
                        frmtmb:::block_col_coefs(bk), sep = "|"))
  }
  out
}
nl <- function(...) structure(list(...), nl = TRUE)
cases <- list(
  list(f = list(y ~ (1 | h), alpha ~ (1 | g), sigma ~ (1 | ga)),
       fb = brms::skew_normal(), ff = frmtmb::skew_normal()),
  list(f = list(cnt ~ (1 | h), zi ~ (1 | g)),
       fb = brms::zero_inflated_poisson(), ff = frmtmb::zero_inflated_poisson()),
  list(f = nl(y ~ b * exp(a * x), a ~ 1 + (1 | h), b ~ 1 + (1 | g) + (1 | ga)),
       fb = gaussian(), ff = gaussian()),
  list(f = nl(y ~ b * exp(a * x), b ~ 1 + (1 | h), a ~ 1 + (1 | ga)),
       fb = gaussian(), ff = gaussian())
)
for (cs in cases) {
  b <- try_(brms_ids(cs$f, cs$fb))
  r <- try_(frm_ids(cs$f, cs$ff))
  cat(sprintf("\n%s %s\n", if (identical(b, r)) "SAME" else "DIFF",
              paste(vapply(cs$f, deparse1, ""), collapse = " ; ")))
  cat(" brms:  ", paste(b, collapse = "  "), "\n")
  if (!identical(b, r)) cat(" frmtmb:", paste(r, collapse = "  "), "\n")
}
sink()
