# Fix round, coordinator's question: on the tight-sd poly(x, 2) fit, are
# theta_2 and theta_4's SEs at the clamp (0.497, 0.660) or the reviewer's
# Richardson references at the stop (0.017, 0.203) the right ones? The
# profile of each (every other parameter re-optimized, theta_5 and
# theta_6 free) at steps of 0.25, 0.5 and 1 times each candidate SE.
#   CODE8_ARM=after Rscript dev/code8-tightsd-prof.R
source("dev/code8-common.R")
env <- new.env(parent = asNamespace("frmtmb"))
for (ex in parse("tests/testthat/test-brms-names-specials.R")) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-"))) eval(ex, env)
}
d <- env$bns_prior_data()
form <- bf(y ~ poly(x, 2) + (1 + poly(x, 2) | g))
p <- set_prior("normal(2, 0.001)", class = "sd", coef = "polyx22",
               group = "g")
f <- quiet(frm(form, data = d, prior = p))
nm <- fns$outer_par_names(f)
pp <- f$opt$par
f0 <- f$opt$objective
obj <- f$obj
prof <- function(j, t) {
  q0 <- pp; q0[j] <- pp[j] + t
  fn <- function(q) { z <- q0; z[-j] <- q; as.numeric(obj$fn(z)) }
  gr <- function(q) { z <- q0; z[-j] <- q; drop(obj$gr(z))[-j] }
  o <- tryCatch(quiet(nlminb(pp[-j], fn, gr, control = list(
    iter.max = 2000, eval.max = 2000, rel.tol = 1e-12))),
    error = function(e) list(objective = NA))
  o$objective - f0
}
cand <- list(theta_2 = c(clamp = 0.497, richardson = 0.0172),
             theta_4 = c(clamp = 0.660, richardson = 0.203),
             theta_1 = c(clamp = 0.241, richardson = 0.2145))
for (nmj in names(cand)) {
  j <- match(nmj, nm)
  for (k in names(cand[[nmj]])) {
    se <- cand[[nmj]][[k]]
    l <- sapply(c(-1, -0.5, -0.25, 0.25, 0.5, 1), function(m) prof(j, m * se))
    cat(sprintf("%-8s %-10s SE %.4g: loss at -1,-.5,-.25,.25,.5,1 SE: %s\n",
                nmj, k, se, paste(format(l, digits = 3), collapse = " ")))
  }
}
cat("(a quadratic profile with that SE loses 0.5, 0.125, 0.031 there)\n")
