source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
library(RTMB)

cat("RTMB exports:\n")
print(sort(getNamespaceExports("RTMB")))

cat("\nS3 methods for advector:\n")
print(sort(as.character(utils::methods(class = "advector"))))

f <- MakeTape(function(x) {
  for (nm in c("as.numeric", "as.double", "as.vector", "unclass")) {
    v <- try(do.call(nm, list(x)), silent = TRUE)
    cat(" ", nm, ": ")
    if (inherits(v, "try-error")) cat("ERROR:", conditionMessage(attr(v, "condition")), "\n")
    else cat(class(v), " -> ", paste(format(utils::head(v, 4)), collapse = " "), "\n")
  }
  sum(x)
}, c(2, 5))
invisible(f(c(2, 5)))
