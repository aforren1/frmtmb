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

# independent reference and its Hessian at the optimum
lp <- function(eta, r) sum(dbinom(d$y[r], 1, plogis(eta), log = TRUE))
ref <- function(p, k = 40) {
  pl <- ns$imp_par_list(fq$frame$par_template, p, fq$frame$map)
  s <- exp(pl$theta)
  nest_ref(pl$beta[1] + pl$beta[2] * d$x, rep(1, nrow(d)), rep(1, nrow(d)),
           as.integer(d$school), as.integer(d$cls_uni), s[1], s[2], lp, k)
}
r0 <- ref(par)
cat(sprintf("REF logLik err %.3e (laplace at par err %.3e)\n", abs(lq - r0),
            abs(lap_at(fq, par) - r0)))
h <- 2e-3
p <- length(par)
H <- matrix(0, p, p)
for (i in seq_len(p)) for (j in i:p) {
  e <- function(a, b) { q <- par; q[i] <- q[i] + a; q[j] <- q[j] + b; ref(q) }
  H[i, j] <- H[j, i] <- if (i == j) {
    (e(h, 0) - 2 * r0 + e(-h, 0)) / h^2
  } else {
    (e(h, h) - e(h, -h) - e(-h, h) + e(-h, -h)) / (4 * h^2)
  }
}
se_ref <- sqrt(diag(solve(-H)))
cat("SE ref:", format(se_ref, digits = 6), "\n")
v <- try(vcov(fq), silent = TRUE)
print(v)
sdr <- ns$sdr_of(fq)
se_sdr <- sqrt(diag(sdr$cov.fixed))
cat("SE sdr:", format(se_sdr, digits = 6), "\n")
cat("SE ratio sdr/ref:", format(se_sdr / se_ref, digits = 6), "\n")
se_lap <- sqrt(diag(solve(optimHess(par, fl$obj$fn, fl$obj$gr))))
cat("SE laplace-objective-at-quad-par / ref:", format(se_lap / se_ref, digits = 6), "\n")

# confint: Wald and profile; count which objective the profile refits use
cw <- try(confint(fq), silent = TRUE)
print(cw)
calls <- new.env(); calls$nest <- 0; calls$gk <- 0
trace("quad_nest_fit", quote(calls$nest <- calls$nest + 1),
      print = FALSE, where = ns)
trace("quad_fit", quote(calls$gk <- calls$gk + 1), print = FALSE,
      where = ns)
cp <- try_msg(confint(fq, method = "profile", parm = "x"))
untrace("quad_nest_fit", where = ns); untrace("quad_fit", where = ns)
print(cp); cat("profile warns:", attr(cp, "warns"), "\n")
cat("profile refits via quad_nest_fit:", calls$nest, " via quad_fit:", calls$gk, "\n")
# the profile endpoint: the reference likelihood with beta_x fixed there,
# maximized over the rest, against qchisq(0.95, 1) / 2
if (!inherits(cp, "error")) {
  ix <- which(names(par) == "beta")[2]
  for (bnd in as.numeric(cp[1, 1:2])) {
    obj <- function(q) { pp <- par; pp[-ix] <- q; pp[ix] <- bnd; -ref(pp, 30) }
    o <- optim(par[-ix], obj, method = "BFGS", control = list(reltol = 1e-10))
    cat(sprintf("PROFILE endpoint %.5f 2*drop %.4f (target %.4f)\n", bnd,
                2 * (o$value + ref(par, 30)), qchisq(0.95, 1)))
  }
}

# ranef: under quadrature the modes come from the Laplace inner solve at
# the quadrature optimum; compare with the one-block route's convention
rq <- ranef(fq)
str(rq, max.level = 3)
inner <- ns$solved_par_list(RTMB::MakeADFun(ns$build_objective(fq$frame),
                                            fq$frame$par_template, random = "b",
                                            map = fq$frame$map, silent = TRUE),
                            par)
cat("ranef == Laplace modes at quad par:",
    isTRUE(all.equal(unname(inner$b), unname(fq$estimates$b))), "\n")
# the posterior mean under the exact integral, for one outer level, by
# quadrature in plain R, against the reported mode
pl <- ns$imp_par_list(fq$frame$par_template, par, fq$frame$map)
s <- exp(pl$theta)
eta0 <- pl$beta[1] + pl$beta[2] * d$x
a1 <- which(d$school == "1")
g <- gh_phys(60)
fa <- function(u) {
  dnorm(u, 0, s[1], log = TRUE) + sum(vapply(unique(d$cls_uni[a1]), function(c) {
    rc <- a1[d$cls_uni[a1] == c]
    adapt_int(function(v) dnorm(v, 0, s[2], log = TRUE) + lp(eta0[rc] + u + v, rc), 40)
  }, 1))
}
uu <- seq(-8, 8, length.out = 801) * s[1]
lw <- vapply(uu, fa, 1)
pm <- sum(uu * exp(lw - max(lw))) / sum(exp(lw - max(lw)))
psd <- sqrt(sum((uu - pm)^2 * exp(lw - max(lw))) / sum(exp(lw - max(lw))))
cat(sprintf("school 1: reported mode %.5f, exact posterior mean %.5f sd %.5f\n",
            fq$estimates$b[fq$frame$re_blocks[[1]]$b_idx[1]], pm, psd))
rq1 <- try(ranef(f1$frame %||% f1), silent = TRUE)
# Est.Error columns, nested and one-block quadrature and Laplace
print(head(rq[[1]][, , 1])); print(head(ranef(fl)[[1]][, , 1]))
if (!inherits(f1, "error")) print(head(ranef(f1)[[1]][, , 1]))

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
