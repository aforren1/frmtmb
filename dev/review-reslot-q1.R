# Reviewer short check, item 1: the named list is gone; no spelling left
# reaches a position cov_class_coefs() lists; lognormal on sd equals
# normal on its log with the Jacobian. Usage: Rscript
# dev/review-reslot-q1.R after. Seed 1.
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
q <- function(expr) suppressWarnings(suppressMessages(expr))
msg <- function(expr) tryCatch({ force(expr); "ACCEPTED" },
                               error = function(e) conditionMessage(e))
dd <- re_data(1L)
f3 <- bf(y ~ x + (1 + x + z | g))
fit <- q(frm(f3, data = dd))
fitd <- q(frm(f3, data = dd, dry_run = "objective"))
show <- function(lab, m) cat(sprintf("-- %s\n   %s\n", lab, substr(m, 1, 400)))

cat("== former spellings\n")
show("frm(prior = list(theta_4 = ))",
     msg(q(frm(f3, data = dd, prior = list(theta_4 = prior_normal(0, 1))))))
show("frm_sample resolve, list(theta_4 = )",
     msg(q(frmtmb.sample:::sample_resolve_priors(
       fitd, list(theta_4 = prior_normal(0, 1))))))
show("frm(prior = list(beta = )) whole component",
     msg(q(frm(f3, data = dd, prior = list(beta = prior_normal(0, 5))))))
show("frm(prior = plain data.frame)",
     msg(q(frm(f3, data = dd, prior = data.frame(prior = "normal(0, 1)",
                                                 class = "b")))))
show("mixed: list(set_prior(...), theta_4 = prior_normal())",
     msg(q(frm(f3, data = dd, prior = list(set_prior("normal(0, 1)", class = "b"),
                                           theta_4 = prior_normal(0, 1))))))
show("mixed: c(set_prior(...), list(theta_4 = ))",
     msg(q(frm(f3, data = dd, prior = c(set_prior("normal(0, 1)", class = "b"),
                                        list(theta_4 = prior_normal(0, 1)))))))
show("update(fit, prior = list(theta_4 = ))",
     msg(q(update(fit, prior = list(theta_4 = prior_normal(0, 1))))))
fitp <- q(frm(f3, data = dd, prior = set_prior("normal(0, 1)", class = "b")))
show("update(MAP fit, prior = list(theta_1 = ))",
     msg(q(update(fitp, prior = list(theta_1 = prior_normal(0, 1))))))
show("frm_sample(fit, prior = list(theta_1 = )) on the sampler",
     msg(q(frm_sample(fit, prior = list(theta_1 = prior_normal(0, 1)),
                      chains = 1, iter = 20, refresh = 0))))
show("frm(prior = 5)", msg(q(frm(f3, data = dd, prior = 5))))

cat("\n== no spelling reaches a class-reached position\n")
des <- frmtmb:::prior_design(f3, dd, NULL, list())
cc <- frmtmb:::cov_class_coefs(des$spec, des$frame)
for (nm in cc$name) {
  m1 <- msg(frmtmb:::resolve_priorlist(des, set_prior("normal(0, 1)",
                                                      class = "theta",
                                                      coef = nm)))
  # as.brmsprior() of a table with a theta row for it
  tb <- data.frame(prior = "normal(0, 1)", class = "theta", coef = nm,
                   group = "", resp = "", dpar = "", nlpar = "",
                   lb = NA, ub = NA, stringsAsFactors = FALSE)
  m2 <- msg(q(frm(f3, data = dd, prior = as.brmsprior(tb))))
  m3 <- msg(q(frmtmb.sample:::sample_resolve_priors(fitd, as.brmsprior(tb))))
  cat(sprintf("  %-8s set_prior: %-9s as.brmsprior frm: %-9s frm_sample: %s\n",
              nm, if (grepl("class of its own", m1)) "refused" else m1,
              if (grepl("class of its own", m2)) "refused" else substr(m2, 1, 60),
              if (grepl("class of its own", m3)) "refused" else substr(m3, 1, 60)))
}
show("class-wide theta",
     msg(frmtmb:::resolve_priorlist(des, set_prior("normal(0, 1)",
                                                   class = "theta"))))
# update() of a fit's prior: start from a class prior, update with theta
show("update(fitp, prior = theta_4 set_prior)",
     msg(q(update(fitp, prior = set_prior("normal(0, 1)", class = "theta",
                                          coef = "theta_4")))))
# a fit made with a class-wide theta on toep, updated: positions stay
# the unreached ones
dd$tf <- factor((dd$t - 1) %% 3)
ft <- bf(y ~ x + toep(0 + tf | g))
fitt <- q(frm(ft, data = dd, prior = set_prior("normal(0, 2)", class = "theta")))
r <- frmtmb:::resolve_prior_input(fitt, fitt$prior)
cat("  toep fit's class-wide theta reaches:",
    vapply(r$entries, function(e) paste(e$comp, paste(e$idx, collapse = ",")), ""),
    "\n")
fitt2 <- q(update(fitt, prior = set_prior("normal(0, 1)", class = "b")))
r2 <- frmtmb:::resolve_prior_input(fitt2, fitt2$prior)
cat("  after update(prior = b):",
    vapply(r2$entries, function(e) paste(e$comp, paste(e$idx, collapse = ","),
                                         e$dist$kind), ""), "\n")

cat("\n== lognormal(0, s) on class sd vs normal(0, s) on log sd\n")
f1 <- bf(y ~ x + (1 | g))
fit1 <- q(frm(f1, data = dd, dry_run = "objective"))
for (s in c(0.5, 2)) {
  ri <- frmtmb:::resolve_prior_input(
    fit1, set_prior(sprintf("lognormal(0, %s)", s), class = "sd"))
  e <- Filter(function(e) identical(e$comp, "theta"), ri$entries)
  th <- c(-3, -0.4, 0, 1.3, 4)
  got <- vapply(th, function(t) frmtmb:::prior_logdens(t, e[[1]]$dist,
                                                       e[[1]]$scale), 0)
  ref <- dnorm(th, 0, s, log = TRUE)
  cat(sprintf("  s = %s: scale %s, max |got - dnorm(theta)| / max|ref| = %.3g\n",
              s, e[[1]]$scale, max(abs(got - ref)) / max(abs(ref))))
}
