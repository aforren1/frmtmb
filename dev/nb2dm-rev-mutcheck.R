# Reviewer: did the noeta and nopost mutations of dev/nb2dm-rev-mut.R
# take effect? And the case where .eta_shape matters (log sigma < -709).
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
ns <- asNamespace("frmtmb")
src <- deparse(get("fam_reparam", ns))
cat("nopost pattern found in deparse:",
    any(grepl("for (k in names(v)) v[[k]] <- reparam_slot(v[[k]], to_inner)",
              src, fixed = TRUE)), "\n")
print(grep("for [(]k", src, value = TRUE))
nb2 <- brmsfamily("negbinomial2")
tape_at <- function(fam, ls) {
  y <- c(0, 2, 9)
  p0 <- c(rep(log(3), 3), rep(ls, 3))
  tp <- RTMB::MakeTape(function(p) {
    frmtmb:::row_lpdf(fam, y, y, list(mu = exp(p[1:3]), sigma = exp(p[4:6]),
                                      .eta_mu = p[1:3], .eta_sigma = p[4:6]),
                      list(), NULL)
  }, p0)
  as.numeric(tp(p0))
}
g <- get("fam_negbinomial2", ns)
body(g) <- do.call(substitute, list(body(g), list(
  robust_logmu = quote(function(...) NULL))))
environment(g) <- ns
nb2m <- g()
for (ls in c(-40, -700, -720, -800, 720)) {
  cat("log sigma", ls, ": lane", format(tape_at(nb2, ls), digits = 10),
      "| noeta", format(tape_at(nb2m, ls), digits = 10), "\n")
}
cat("dpois(y, 3, log):", format(dpois(c(0, 2, 9), 3, log = TRUE),
                                digits = 10), "\n")
