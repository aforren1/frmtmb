# Where the frozen nested tape goes below its own optimum: fix beta_x,
# minimize the frozen tape over the rest, then re-anchor there.
source("C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev/review-nestquad-common.R")
set.seed(20261107)
ns_ <- 14
d <- do.call(rbind, lapply(seq_len(ns_), function(s) {
  do.call(rbind, lapply(seq_len(sample(1:4, 1)), function(c) {
    data.frame(school = s, cls = c, x = rnorm(sample(2:5, 1)))
  }))
}))
d$school <- factor(d$school)
d$cls_uni <- factor(paste0(d$school, "_", d$cls))
ua <- rnorm(ns_, 0, 0.9); ui <- rnorm(nlevels(d$cls_uni), 0, 1)
d$y <- rbinom(nrow(d), 1, plogis(0.2 + 0.7 * d$x + ua[d$school] +
                                   ui[d$cls_uni]))
fq <- frm(bf(y ~ x + (1 | school/cls_uni)) + bernoulli(), data = d,
          quadrature = TRUE)
par <- fq$opt$par
lp <- function(eta, r) sum(dbinom(d$y[r], 1, plogis(eta), log = TRUE))
ref <- function(p) nest_ref(p[1] + p[2] * d$x, rep(1, nrow(d)), rep(1, nrow(d)),
                            as.integer(d$school), as.integer(d$cls_uni),
                            exp(p[3]), exp(p[4]), lp, 40)
for (bx in c(0.69, 0.9, 1.17)) {
  o <- nlminb(par[-2], function(q) { p <- par; p[-2] <- q; p[2] <- bx; fq$obj$fn(p) })
  p <- par; p[-2] <- o$par; p[2] <- bx
  cat(sprintf("bx %.2f frozen-min %.6f at %s | reanchored %.6f ref %.6f\n", bx,
              o$objective, paste(format(o$par, digits = 4), collapse = " "),
              -nest_at(fq, p, c(11, 15)), -ref(p)))
}
# the reference profile at the same beta_x, maximized over the rest
for (bx in c(0.69, 0.9)) {
  o <- optim(par[-2], function(q) { p <- par; p[-2] <- q; p[2] <- bx; -ref(p) },
             method = "BFGS")
  cat(sprintf("REFPROFILE bx %.2f min %.6f (fit %.6f) 2*drop %.4f\n", bx,
              o$value, -ref(par), 2 * (o$value + ref(par))))
}
# confint(method = "uniroot") reads the same frozen tape
cu <- try(confint(fq, method = "uniroot", parm = "x"))
print(cu)
cat("Wald:", confint(fq)[2, ], "\n")
