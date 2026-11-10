# Lane hollow: why a model frame rebuilt from the raw data with the
# frame's terms does not reproduce a frame with poly() and offset().
.libPaths(c("C:/Users/adf44/source/r/wt-hollow-lib",
            "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(4691)
dd <- data.frame(x = stats::rnorm(40), age = stats::runif(40, 1, 3),
                 g = factor(rep(1:5, 8)))
dd$y <- stats::rpois(40, exp(0.2 + 0.3 * dd$x + log(dd$age)))
fit <- frm(bf(y ~ poly(x, 2) + offset(log(age)) + (1 | g)),
           family = poisson(), data = dd)
mf <- fit$data
tt <- attr(mf, "terms")
print(tt)
re <- stats::model.frame(tt, data = dd, na.action = stats::na.pass)
cat("frame names:", names(mf), "\n")
cat("rebuilt names:", names(re), "\n")
for (v in intersect(names(mf), names(re))) {
  cat(v, ": all.equal ->",
      format(all.equal(unclass(re[[v]]), unclass(mf[[v]]), tolerance = 0,
                       check.attributes = FALSE)), "\n")
}
str(fit$frame[["raw_vars"]]); str(fit$frame[["data_vars"]])
