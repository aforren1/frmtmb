# Fix round 1, M1 and m5: every route that reads a nested fit's tape away
# from its estimate refuses by class; Wald and bootstrap still run.
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
f1 <- frm(bf(y ~ x + (1 | school)) + bernoulli(), data = d, quadrature = TRUE)
show <- function(tag, expr) {
  r <- try_msg(expr)
  cat(sprintf("ROUTE %-40s %s\n", tag,
              if (inherits(r, "error")) paste("REFUSED", paste(class(r), collapse = "/"),
                                              substr(conditionMessage(r), 1, 90))
              else paste("RAN", class(r)[1])))
}
show("confint profile", confint(fq, method = "profile", parm = "x"))
show("confint uniroot", confint(fq, method = "uniroot", parm = "x"))
show("profile()", profile(fq))
show("hypothesis profile", hypothesis(fq, "x > 0", method = "profile"))
show("conditional_effects band profile",
     conditional_effects(fq, effects = "x", band = "profile"))
show("ranef condVar", ranef(fq, condVar = TRUE))
show("confint wald", confint(fq))
show("confint boot nsim 20", confint(fq, method = "boot", nsim = 20, seed = 1))
show("hypothesis wald", hypothesis(fq, "x > 0"))
show("conditional_effects default", conditional_effects(fq, effects = "x"))
# the one-block route keeps its profile
show("one-block confint profile", confint(f1, method = "profile", parm = "x"))
show("one-block ranef condVar", ranef(f1, condVar = TRUE))
cat("DONE\n")
show("profile(parm = 'x')", profile(fq, parm = "x"))
show("profile(parm = 2)", profile(fq, parm = 2))
