# Attack 4: post-fit methods on a nested quadrature fit. Each is checked
# against the quadrature objective or an independent quantity.
#   Rscript dev/review-nestquad-postfit.R
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
cat("rows", nrow(d), "\n")
fo <- bf(y ~ x + (1 | school/cls_uni)) + bernoulli()
fq <- try_msg(frm(fo, data = d, quadrature = TRUE))
cat("fit warns:", attr(fq, "warns"), "\n")
fl <- frm(fo, data = d)
f1 <- try_msg(frm(bf(y ~ x + (1 | school)) + bernoulli(), data = d,
                  quadrature = TRUE))
par <- fq$opt$par
lq <- as.numeric(logLik(fq))
cat(sprintf("logLik %.10f df %s  -objective %.10f  AIC %.10f  -2LL+2df %.10f\n",
            lq, attr(logLik(fq), "df"), -fq$opt$objective, AIC(fq),
            -2 * lq + 2 * attr(logLik(fq), "df")))
cat(sprintf("BIC %.6f nobs %d\n", BIC(fq), nobs(fq)))

lp <- function(eta, r) sum(dbinom(d$y[r], 1, plogis(eta), log = TRUE))
ref <- function(p, k = 40) {
  pl <- ns$imp_par_list(fq$frame$par_template, p, fq$frame$map)
  s <- exp(pl$theta)
  nest_ref(pl$beta[1] + pl$beta[2] * d$x, rep(1, nrow(d)), rep(1, nrow(d)),
           as.integer(d$school), as.integer(d$cls_uni), s[1], s[2], lp, k)
}
# confint(method = "root") reads the same frozen tape
cr <- try_msg(confint(fq, method = "root", parm = "x"))
print(cr); cat("root warns:", attr(cr, "warns"), "\n")
for (nm in c("fq", "f1", "fl")) {
  f <- get(nm)
  r <- try_msg(ranef(f, condVar = TRUE))
  if (inherits(r, "error")) cat("ranef condVar", nm, "ERROR", conditionMessage(r), "\n") else {
    cat("ranef condVar", nm, ":\n"); str(r[[1]], max.level = 1)
    pv <- attr(r[[1]], "postVar") %||% attr(r[[1]], "condVar")
    if (!is.null(pv)) cat("  first condVar sd:", format(sqrt(as.numeric(pv)[1:3]), digits = 5), "\n")
  }
  if (length(attr(r, "warns"))) cat("  warns:", attr(r, "warns"), "\n")
}
# the Laplace fit's conditional sd at the quad optimum, for comparison
# predict at new levels, simulate, update, anova
nd <- data.frame(x = c(0, 1), school = factor(c("1", "99")),
                 cls_uni = factor(c("1_1", "99_1")))
pr <- try_msg(predict(fq, newdata = nd, allow_new_levels = TRUE))
print(pr); cat("predict warns:", attr(pr, "warns"), "\n")
fi <- try_msg(fitted(fq, newdata = nd, allow_new_levels = TRUE))
print(fi); cat("fitted warns:", attr(fi, "warns"), "\n")
sm <- try_msg(simulate(fq, nsim = 2, seed = 1))
cat("simulate:", class(sm), dim(sm), "\n")
fu <- try_msg(update(fq, control = frmtmb_control(quad_nodes = c(21, 21))))
cat(sprintf("update quad_nodes 21: logLik %.10f (default %.10f) ref %.10f\n",
            as.numeric(logLik(fu)), lq, ref(fu$opt$par)))
fu2 <- try_msg(update(fq, quadrature = FALSE))
cat("update quadrature=FALSE quadrature flag:", fu2$quadrature,
    sprintf("logLik %.8f", as.numeric(logLik(fu2))), "\n")
fu3 <- try_msg(update(fq, . ~ . + I(x^2)))
cat("update formula: quadrature flag", fu3$quadrature, "random",
    length(fu3$obj$env$random), "\n")
an <- try_msg(anova(f1, fq))
print(an); cat("anova warns:", attr(an, "warns"), "\n")
an2 <- try_msg(anova(fl, fq))
print(an2); cat("anova laplace-vs-quad warns:", attr(an2, "warns"), "\n")
sq <- try_msg(summary(fq)); cat("summary warns:", attr(sq, "warns"), "\n")
dg <- try_msg(diagnose(fq)); cat("diagnose:", class(dg), attr(dg, "warns"), "\n")
cat("DONE\n")
