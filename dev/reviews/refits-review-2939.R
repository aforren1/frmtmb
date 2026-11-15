# Review of lane refits, item 2939: (… || gr(g, cov = A)) and gr(prec =)
# before and after; names, logLik, re_formula, by =, plain ||.
# Seed 11, 12 levels, A[i, j] = 0.5^|i - j|, n = 120.
source("C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews/refits-review-lib.R")
set.seed(11)
ng <- 12
A <- 0.5 ^ abs(outer(1:ng, 1:ng, "-"))
dimnames(A) <- list(paste0("g", 1:ng), paste0("g", 1:ng))
Q <- solve(A)
d <- data.frame(g = factor(rep(paste0("g", 1:ng), each = 10),
                           levels = paste0("g", 1:ng)),
                f = factor(rep(c("a", "b", "c"), 40)),
                f2 = factor(rep(c("a", "b"), 60)),
                f4 = factor(rep(c("a", "b", "c", "d"), 30)),
                h = factor(rep(c("u", "v"), each = 60)),
                x = rnorm(120), z = rnorm(120))
u <- MASS::mvrnorm(3, rep(0, ng), A)
d$y <- 1 + 0.5 * d$x + u[1, as.integer(d$g)] +
  (d$f == "b") * u[2, as.integer(d$g)] + rnorm(120, 0, 0.7)
show <- function(fit) {
  if (is.character(fit)) return(fit)
  vc <- VarCorr(fit)
  list(logLik = sprintf("%.10f", as.numeric(logLik(fit))),
       npar = attr(logLik(fit), "df"),
       vc_names = names(vc),
       vc_sd = lapply(vc, function(v) round(v$sd[, 1], 5)),
       re_names = names(ranef(fit)),
       re_cols = lapply(ranef(fit), function(r) dimnames(r)[[3]]),
       sd_vars = grep("^sd_|^cor_", variables(fit), value = TRUE),
       r_vars_head = head(grep("^r_", variables(fit), value = TRUE), 4),
       n_r_vars = length(grep("^r_", variables(fit), value = TRUE)))
}
forms <- list(
  num_cov = y ~ x + (1 + x || gr(g, cov = A)),
  num_cov_cor = y ~ x + (1 + x | gr(g, cov = A, cor = FALSE)),
  num_prec = y ~ x + (1 + x || gr(g, prec = Q)),
  f_cov = y ~ x + (1 + f || gr(g, cov = A)),
  f0_cov = y ~ x + (0 + f || gr(g, cov = A)),
  f2_cov = y ~ x + (1 + f2 || gr(g, cov = A)),
  f4_0_cov = y ~ x + (0 + f4 || gr(g, cov = A)),
  mix_cov = y ~ x + (1 + z + f || gr(g, cov = A)),
  f_prec = y ~ x + (0 + f || gr(g, prec = Q)),
  plain_f = y ~ x + (1 + f || g),
  plain_num = y ~ x + (1 + x || g),
  by_num_cov = y ~ x + (1 + x || gr(g, cov = A, by = h)),
  by_f_cov = y ~ x + (0 + f2 || gr(g, cov = A, by = h)),
  by_num_plain = y ~ x + (1 + x || gr(g, by = h))
)
fits <- list()
for (nm in names(forms)) {
  fits[[nm]] <- go(paste(nm, deparse1(forms[[nm]])),
                   frm(forms[[nm]], data = d, data2 = list(A = A, Q = Q)))
  go(paste(nm, "summary"), show(fits[[nm]]))
}
cat("\n== re_formula naming one piece of the numeric || term ==\n")
fit <- fits[["num_cov"]]
if (!is.character(fit)) {
  nd <- d[1:4, ]
  go("predict re_formula = ~ (1 | gr(g, cov = A))",
     predict(fit, newdata = nd, re_formula = ~ (1 | gr(g, cov = A))))
  go("predict re_formula = ~ (1 + x || gr(g, cov = A))",
     predict(fit, newdata = nd, re_formula = ~ (1 + x || gr(g, cov = A))))
  go("predict re_formula = NULL", predict(fit, newdata = nd))
  go("predict re_formula = NA", predict(fit, newdata = nd, re_formula = NA))
  go("prior sd coef x group g",
     { f2 <- frm(forms$num_cov, data = d, data2 = list(A = A),
                 prior = set_prior("normal(0, 0.01)", class = "sd",
                                   group = "g", coef = "x"))
       round(VarCorr(f2)[[1]]$sd, 5) })
  go("prior_summary sd rows",
     { ps <- prior_summary(fit); ps[ps$class == "sd", c("class", "coef", "group")] })
}
cat("\n== new levels ==\n")
fit <- fits[["f_cov"]]
if (!is.character(fit)) {
  nd <- d[1:3, ]
  go("fitted(newdata) f_cov", round(fitted(fit, newdata = nd)[, 1], 6))
}
cat("\n== plain || bit-identity values ==\n")
for (nm in c("plain_f", "plain_num", "by_num_plain")) {
  if (!is.character(fits[[nm]])) {
    cat(nm, "fixef:", sprintf("%.15g", fixef(fits[[nm]], flatten = TRUE)), "\n")
    cat(nm, "logLik:", sprintf("%.15g", as.numeric(logLik(fits[[nm]]))), "\n")
  }
}
cat("\n== logLik full precision ==\n")
for (nm in names(fits)) if (!is.character(fits[[nm]])) {
  cat(nm, sprintf("%.12f", as.numeric(logLik(fits[[nm]]))), "\n")
}
saveRDS(lapply(fits, function(f) if (is.character(f)) f else
  list(ll = as.numeric(logLik(f)), fe = fixef(f, flatten = TRUE))),
  paste0("C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews/refits-review-2939-", arm, ".rds"))
