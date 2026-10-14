p <- "extensions/frmtmb.eam/R/wiener-family.R"
x <- paste(readLines(p), collapse = "\n")
old <- "#' Under an identity link nothing keeps `bs` above zero, `bias` inside
#' (0, 1) or `ndt` above zero. The density is `NaN` outside those
#' ranges, which the optimizer treats as a rejected step."
new <- "#' Under an identity link nothing keeps `bs` above zero, `bias` inside
#' (0, 1) or `ndt` above zero. The density is `NaN` outside those
#' ranges, on the tape the optimizer evaluates as well as on doubles, and
#' the optimizer treats it as a rejected step. On data with a bias of
#' 0.97 (seeds 20261101 to 20261105, n = 80), `link_bias = \"identity\"`
#' fits reach the logit fits' maxima on all five.
#'
#' A brms family object is read with its links. `brms::wiener()` and
#' `brms::brmsfamily(\"wiener\")` carry brms's `link_ndt = \"log\"`, so
#' `family = brms::wiener()` fits brms's log link for `ndt`, while
#' `family = \"wiener\"` and this package's `wiener()` fit the bounded
#' default. With both packages attached, a bare `wiener()` is the one
#' attached last: brms's gives the log link, this package's the bounded
#' one. On a scalar `ndt` the two are the same model with the same
#' maximum; under `ndt ~ x` the coefficients are on different scales, and
#' the log link refuses `max_ndt`, `ndt_group()`, `variability` and
#' `contaminant = TRUE`. Write `frmtmb.eam::wiener()` to be sure which
#' you have."
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
writeLines(strsplit(x, "\n")[[1]], p)
