# Reviewer: is the by-split failure of quad_passes() its pass order?
# Same arm as dev/fuzzfix-rev-quad-shapes.R by_split_inner4, with the
# passes ordered inner blocks first (this process only).
#   Rscript dev/fuzzfix-rev-order.R <lib> <ng> <lane|inner_first>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
ng <- as.integer(a[2])
if (a[3] == "inner_first") {
  qp <- function(blocks) {
    outer <- vapply(blocks, function(bk) bk$group_name == "ga", NA)
    left <- sort(unlist(lapply(blocks, `[[`, "b_idx")))
    lapply(blocks[order(outer)], function(bk) {
      w <- match(bk[["b_idx"]], left)
      left <<- left[-w]
      list(b = structure(list(dim = 1, adaptive = FALSE, debug = FALSE,
                              method = "marginal_gk", which = w),
                         class = "GK"))
    })
  }
  utils::assignInNamespace("quad_passes", qp, "frmtmb")
}
set.seed(31)
m <- 3
d <- expand.grid(r = seq_len(m), gb = factor(c("i", "ii")),
                 ga = factor(seq_len(ng)))
n <- nrow(d)
d$x <- rnorm(n)
d$cell <- interaction(d$ga, d$gb, drop = TRUE)
d$cellu <- factor(paste0("c", as.integer(d$cell)))
d$f4 <- factor(as.integer(d$cell) %% 4)
eta <- 0.3 + 0.5 * d$x + rnorm(ng, 0, 0.8)[d$ga] +
  rnorm(nlevels(d$cell), 0, 0.6)[d$cell]
d$y <- rbinom(n, 1, plogis(eta))
t0 <- proc.time()[[3]]
f <- tryCatch(suppressWarnings(suppressMessages(
  frm(bf(y ~ x + (1 | ga) + (1 | gr(cellu, by = f4))), family = bernoulli(),
      data = d, quadrature = TRUE))), error = function(e) e)
el <- proc.time()[[3]] - t0
if (inherits(f, "error")) {
  cat(sprintf("ORDER %s ng %d rows %d ERROR (%.1fs) %s\n", a[3], ng, n, el,
              substr(conditionMessage(f), 1, 120)))
} else {
  fl <- suppressWarnings(suppressMessages(
    frm(bf(y ~ x + (1 | ga) + (1 | gr(cellu, by = f4))), family = bernoulli(),
        data = d)))
  cat(sprintf(paste0("ORDER %s ng %d rows %d ok (%.1fs) left_random %d ",
                     "logLik quad %.8f laplace %.8f\n"), a[3], ng, n, el,
              length(f$obj$env$random), as.numeric(logLik(f)),
              as.numeric(logLik(fl))))
}
