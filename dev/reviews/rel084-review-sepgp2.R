# rel084 release review, finding (a) continued: why no gp message beside the
# separation warning, and a quasi-separated case where the field stays
# (bernoulli, one factor level all 1, field from ls 0.7 at n 200).
# Prints per fit: code, white flag, sdgp, max |field|, max |eta|, and the
# result of gp_basis_message() called on the fit, with conv_classes as set
# and as if a convergence warning had been recorded.
#   REVLIB=r23 Rscript dev/reviews/rel084-review-sepgp2.R
source("dev/reviews/rel084-review-common.R")
options(width = 170)
ns <- asNamespace("frmtmb")
field <- function(seed, n) {
  set.seed(seed)
  x <- sort(runif(n, 0, 10))
  K <- exp(-outer(x, x, "-")^2 / (2 * 0.7^2)) + diag(1e-8, n)
  list(x = x, f = as.vector(t(chol(K)) %*% rnorm(n)))
}
probe <- function(tag, r) {
  f <- r$v
  if (is.null(f)) { cat(tag, "ERROR\n"); return() }
  th <- f$estimates$theta
  wh <- NA; sdg <- NA
  for (bk in f$frame$re_blocks) if (identical(bk$covstruct, "hsgp")) {
    wh <- ns$gp_basis_white(bk, th[bk$theta_idx])
    sdg <- exp(th[bk$theta_idx][1])
  }
  eta <- tryCatch(max(abs(predict(f, type = "link")[, 1])), error = function(e) NA)
  gm <- tryCatch(ns$gp_basis_message(f), error = function(e) conditionMessage(e))
  old <- f$cache$conv_classes
  f$cache$conv_classes <- "frmtmb_conv_code"
  gm2 <- tryCatch(ns$gp_basis_message(f), error = function(e) conditionMessage(e))
  f$cache$conv_classes <- old
  said <- paste(vapply(r$said, function(s) paste0(s$k, ":", s$cls[[1]],
    if ("frmtmb_boundary_fit" %in% s$cls) "(bnd)" else ""), ""), collapse = " ")
  cat(sprintf("%s | code %s | white %s sdgp %.3g max|eta| %.3g | gp msg now: %s | with conv class: %s | %s\n",
              tag, f$opt$convergence, wh, sdg, eta,
              if (length(gm)) paste(vapply(gm, function(m) paste(m$class, collapse = "/"), ""), collapse = ";") else "none",
              if (length(gm2)) paste(vapply(gm2, function(m) paste(m$class, collapse = "/"), ""), collapse = ";") else "none",
              said))
}
for (s in 1:3) {
  fl <- field(9600 + s, 60)
  d <- data.frame(x = fl$x, y = fl$f + 0.3 * rnorm(60), cc = "right")
  probe(sprintf("SEPA seed %d", s),
        conds_of(frm(y | cens(cc) ~ gp(x, k = 5), family = asym_laplace(), data = d)))
  d$cc <- "none"
  probe(sprintf("NOCENS-gauss seed %d", s),
        conds_of(frm(y ~ gp(x, k = 5), data = d)))
}
for (s in 1:8) {
  fl <- field(9800 + s, 200)
  g <- factor(sample(c("a", "b", "c"), 200, TRUE))
  p <- plogis(1.5 * fl$f)
  y <- rbinom(200, 1, p)
  y[g == "c"] <- 1L
  d <- data.frame(x = fl$x, g = g, y = y)
  probe(sprintf("QSEP seed %d", s),
        conds_of(frm(y ~ g + gp(x, k = 5), family = bernoulli(), data = d)))
}
cat("DONE\n")
