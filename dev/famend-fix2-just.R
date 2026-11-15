# Lane famend, fix round 2: the rows of famend-short-gevclamp.R's "JUST
# INSIDE the clamp" block. That block calls a row inside when xi * x is
# above the clamp on doubles; the tape forms xi through log1p and expm1.
# Each row's value, gradient and Hessian against 0.82.0's gev_parts(),
# split by which side of the clamp the tape's z is on.
#   Rscript dev/famend-fix2-just.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "\n")
src <- readLines("dev/reviews/famend-short-out/fe-base.R")
i0 <- grep("^gev_parts <- function", src)
i1 <- i0 - 1L + grep("^}", src[i0:length(src)])[1]
old_parts <- eval(parse(text = src[i0:i1]))
environment(old_parts) <- ns
fam <- gen_extreme_value()
with_parts <- function(parts, expr) {
  orig <- get("gev_parts", ns)
  unlockBinding("gev_parts", ns)
  assign("gev_parts", parts, ns)
  on.exit({ assign("gev_parts", orig, ns); lockBinding("gev_parts", ns) })
  expr
}
mu0 <- 0.3
ls0 <- -0.2
eval_slot <- function(slot, y, xi) {
  p <- c(mu0, ls0, log1p(xi))
  tp <- RTMB::MakeTape(function(p) {
    fam[[slot]](y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                        .eta_sigma = p[2], .eta_xi = p[3]), list())
  }, p)
  list(v = tp(p), g = tp$jacobian(p), h = tp$jacfun()$jacobian(p))
}
cl <- -1 + 1e-15
tab <- matrix(0, 2, 2, dimnames = list(c("tape inside", "tape outside"),
                                       c("identical", "differ")))
for (xi in c(0.06, 0.3, 1, 3, 10, -0.5)) {
  for (k in 1:8) {
    x <- (cl + k * .Machine$double.eps / 2) / xi
    y <- mu0 + exp(ls0) * x
    xx <- (y - mu0) / exp(ls0)
    if (!(xi * xx > cl)) next
    side <- if (expm1(log1p(xi)) * xx > cl) 1 else 2
    for (slot in c("lpdf", "lcdf", "lccdf")) {
      same <- identical(with_parts(old_parts, eval_slot(slot, y, xi)),
                        eval_slot(slot, y, xi))
      tab[side, 2 - same] <- tab[side, 2 - same] + 1
    }
  }
}
print(tab)
