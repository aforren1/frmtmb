# Reviewer: size of the non-identical numeric outputs (base vs lane).
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
b <- readRDS(file.path(wt, "dev/ordinal-rev-bitwise-base.rds"))
l <- readRDS(file.path(wt, "dev/ordinal-rev-bitwise-lane.rds"))
for (mk in list(c("cum_cloglog_gr", "fixef"), c("cum_cloglog_gr", "vcov"),
                c("cum_cloglog_gr", "coef"), c("mv_ord_gauss", "coef"),
                c("hurdle", "draws"), c("hurdle", "draws_epred"),
                c("hurdle", "draws_vars"))) {
  x <- b[[mk[1]]][[mk[2]]]; y <- l[[mk[1]]][[mk[2]]]
  cat("==", mk, " class", class(x), "\n")
  if (is.numeric(x) && is.numeric(y) && length(x) == length(y)) {
    cat(" dims equal:", identical(dim(x), dim(y)),
        " dimnames equal:", identical(dimnames(x), dimnames(y)),
        " names equal:", identical(names(x), names(y)),
        " max rel diff:", max(abs(x - y) / pmax(abs(x), 1e-300)), "\n")
    if (!identical(dimnames(x), dimnames(y))) {
      print(setdiff(unlist(dimnames(x)), unlist(dimnames(y))))
      print(setdiff(unlist(dimnames(y)), unlist(dimnames(x))))
    }
  } else if (is.list(x)) {
    for (k in names(x)) cat(" ", k, identical(x[[k]], y[[k]]),
                            if (is.numeric(x[[k]])) max(abs(x[[k]] - y[[k]])), "\n")
  } else { print(setdiff(x, y)); print(setdiff(y, x)) }
}
