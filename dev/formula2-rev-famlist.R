# Reviewer: try to falsify the family-list claims.
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(expr)),
           error = function(e) structure(conditionMessage(e), class = "ERR"))
}
show <- function(x) {
  if (inherits(x, "ERR")) paste("ERR:", substr(x, 1, 170)) else x
}
fams <- function(fit) {
  if (inherits(fit, "ERR")) return(show(fit))
  paste(vapply(fit$spec$responses, function(r) r$family$family, ""),
        collapse = ", ")
}
set.seed(8)
n <- 200
d <- data.frame(x = rnorm(n))
d$y <- 1 + d$x + rnorm(n)
d$cnt <- rpois(n, exp(0.5 + 0.3 * d$x))
d$cat <- factor(sample(c("a", "b", "c"), n, TRUE))
d$yb <- rbinom(n, 1, 0.4)
k <- rbinom(n, 1, 0.5)
d$ymix <- ifelse(k == 1, rnorm(n, 3), rnorm(n, -1))
mv <- bf(y ~ x) + bf(cnt ~ x) + set_rescor(FALSE)

a <- q(frm(mv, data = d, family = list(gaussian(), poisson())))
b <- q(frm(bf(y ~ x) + gaussian() + bf(cnt ~ x) + poisson() +
             set_rescor(FALSE), data = d))
cat("order y, cnt:", fams(a), "| logLik identical to + family:",
    identical(logLik(a), logLik(b)), "\n")
a2 <- q(frm(bf(cnt ~ x) + bf(y ~ x) + set_rescor(FALSE), data = d,
            family = list(poisson(), gaussian())))
cat("order cnt, y:", fams(a2), "\n")
cat("named list in the wrong order (positional, as brms):",
    fams(q(frm(mv, data = d, family = list(cnt = poisson(),
                                             y = gaussian())))), "\n")
cat("shorter list:", show(q(frm(mv, data = d, family = list(gaussian())))),
    "\n")
cat("longer list:", show(q(frm(mv, data = d, family = list(gaussian(),
                                                         poisson(),
                                                         gaussian())))),
    "\n")
cat("univariate with a list:", show(q(frm(y ~ x, data = d,
                                         family = list(gaussian())))), "\n")
cat("strings:", fams(q(frm(mv, data = d,
                          family = list("gaussian", "poisson")))), "\n")
cat("functions:", fams(q(frm(mv, data = d,
                            family = list(gaussian, poisson)))), "\n")
cat("with + family on cnt (list fills only y):",
    fams(q(frm(bf(y ~ x) + bf(cnt ~ x) + poisson() + set_rescor(FALSE),
               data = d, family = list(student(), gaussian())))), "\n")
cat("mixture in a list:",
    fams(q(frm(bf(ymix ~ x) + bf(cnt ~ x) + set_rescor(FALSE), data = d,
               family = list(mixture(gaussian(), gaussian()),
                             poisson())))), "\n")
cat("categorical (deferred) in a list:",
    fams(q(frm(bf(cat ~ x) + bf(y ~ x) + set_rescor(FALSE), data = d,
               family = list(categorical(), gaussian())))), "\n")
cat("bernoulli + hurdle_poisson:",
    fams(q(frm(bf(yb ~ x) + bf(cnt ~ x) + set_rescor(FALSE), data = d,
               family = list(bernoulli(), hurdle_poisson())))), "\n")
cat("mvbind with a list:",
    fams(q(frm(bf(mvbind(y, cnt) ~ x) + set_rescor(FALSE), data = d,
               family = list(gaussian(), poisson())))), "\n")
cat("rescor with two gaussians:",
    fams(q(frm(bf(mvbind(y, ymix) ~ x) + set_rescor(TRUE), data = d,
               family = list(gaussian(), gaussian())))), "\n")
cat("list with a NULL entry:",
    show(q(fams(frm(mv, data = d, family = list(NULL, poisson()))))), "\n")
cat("list with a bad entry:",
    show(q(fams(frm(mv, data = d, family = list(gaussian(), "nope"))))), "\n")
cat("brmsfamily in a list:",
    fams(q(frm(mv, data = d, family = list(brmsfamily("gaussian"),
                                           brmsfamily("poisson"))))), "\n")
cat("== other entry points\n")
gp1 <- q(get_prior(mv, data = d, family = list(gaussian(), poisson())))
gp2 <- q(get_prior(bf(y ~ x) + gaussian() + bf(cnt ~ x) + poisson() +
                     set_rescor(FALSE), data = d))
cat("get_prior identical:", identical(gp1, gp2), "\n")
pt1 <- q(par_template(mv, data = d, family = list(gaussian(), poisson())))
pt2 <- q(par_template(bf(y ~ x) + gaussian() + bf(cnt ~ x) + poisson() +
                        set_rescor(FALSE), data = d))
cat("par_template identical:", identical(pt1, pt2), "\n")
cat("frm_simulate:", show(q(class(frm_simulate(mv, data = d,
  family = list(gaussian(), poisson()), newparams = pt1, seed = 1)))), "\n")
cat("update(a, family = list(student(), poisson())):",
    fams(q(update(a, family = list(student(), poisson())))), "\n")
cat("update(b, newdata = d[-1, ]):", fams(q(update(b, newdata = d[-1, ]))),
    "\n")
cat("update(a, newdata = d[-1, ]):", fams(q(update(a, newdata = d[-1, ]))),
    "\n")
cat("DONE\n")
cat("fixed construction: bf(y ~ x) + (bf(cnt ~ x) + poisson()), list(student(), gaussian()):",
    fams(q(frm(bf(y ~ x) + (bf(cnt ~ x) + poisson()) + set_rescor(FALSE),
               data = d, family = list(student(), gaussian())))), "\n")
cat("brms on the same:",
    paste(vapply(brms:::validate_formula(brms::bf(y ~ x) +
      (brms::bf(cnt ~ x) + brms::brmsfamily("poisson")) + brms::set_rescor(FALSE),
      data = d, family = list(brms::student(), brms::brmsfamily("gaussian")))$forms,
      function(f) f$family$family, ""), collapse = ", "), "\n")
