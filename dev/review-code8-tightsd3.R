# Review of lane code8: are (Intercept) and poly(x, 2)1 of the poly
# tight-sd fit determined? Profile each (others re-optimized from the
# stop) at +-0.5 and +-1 of the reference SE.
#   RV_ARM=after Rscript dev/review-code8-tightsd3.R
source("dev/review-code8-common.R")
env <- new.env(parent = asNamespace("frmtmb"))
for (ex in parse("tests/testthat/test-brms-names-specials.R")) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-"))) eval(ex, env)
}
d <- env$bns_prior_data()
form <- bf(y ~ poly(x, 2) + (1 + poly(x, 2) | g))
p <- set_prior("normal(2, 0.001)", class = "sd", coef = "polyx22", group = "g")
f <- quiet(frm(form, data = d, prior = p))
x <- f$opt$par
nm <- ns$outer_par_names(f)
obj <- f$obj
f0 <- obj$fn(x)
for (pn in c("(Intercept)", "poly(x, 2)1", "poly(x, 2)2", "theta_1")) {
  i <- match(pn, nm)
  sc <- c("(Intercept)" = 0.2846, "poly(x, 2)1" = 1.453,
          "poly(x, 2)2" = 1.297, theta_1 = 0.2145)[[pn]]
  out <- sapply(c(-2, -1, -0.5, 0.5, 1, 2), function(m) {
    st <- x; st[i] <- x[i] + m * sc
    oth <- setdiff(seq_along(x), i)
    o <- nlminb(x[oth], function(q) { z <- st; z[oth] <- q; obj$fn(z) },
                function(q) { z <- st; z[oth] <- q; drop(obj$gr(z))[oth] },
                control = list(iter.max = 500, eval.max = 1000))
    o$objective - f0
  })
  cat(sprintf("%-12s profile loss at -2,-1,-0.5,0.5,1,2 ref SE: %s\n", pn,
              paste(signif(out, 3), collapse = " ")))
}
cat("DONE\n")
