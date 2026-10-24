# Fix round 1 check of F2. Usage: Rscript review-splitsurv-saved2.R <save|absent>
# save:   lane lib, frmtmb.survival attached; fits saved to RDS
# absent: rev lib holding only a copy of the lane's frmtmb ahead of
#         rellib-r11 (which has no frmtmb.survival)
arm <- commandArgs(trailingOnly = TRUE)[1]
first <- if (arm == "absent") "C:/Users/adf44/source/r/rev-splitsurv-lib" else
  "C:/Users/adf44/source/r/wt-splitsurv-lib"
.libPaths(c(first, "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
O <- "C:/Users/adf44/source/r/frmtmb-wt-splitsurv/dev/review-splitsurv-out"
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "; survival installed:",
    nzchar(system.file(package = "frmtmb.survival")), "\n")
if (arm == "save") {
  suppressMessages(library(frmtmb.survival))
  set.seed(1)
  dd <- data.frame(x = rnorm(200), g = factor(rep(1:2, 100)))
  dd$time <- rexp(200, exp(-0.5 + 0.7 * dd$x))
  dd$out <- rbinom(200, 1, 0.3)
  fc <- frm(bf(time | cens(out) ~ x), family = cox(), data = dd)
  fr <- frm(bf(time | cens(out) ~ x), family = royston_parmar(df = 2),
            data = dd)
  # a mixture of a core and an extension component
  fm <- tryCatch(suppressWarnings(
    frm(bf(time ~ x), family = mixture(lognormal(), royston_parmar(df = 1)),
        data = dd)), error = function(e) conditionMessage(e))
  cat("mixture fit:", if (is.character(fm)) fm else "ok", "\n")
  if (!is.character(fm)) {
    cat("mixture component packages:",
        attr(family(fm), "frmtmb_component_packages"), "\n")
  }
  # a family a user builds at top level, from the global environment
  ufam <- frmtmb_family(
    "ugauss", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      -0.5 * ((y - dpars$mu) / dpars$sigma)^2 - log(dpars$sigma) -
        0.5 * log(2 * pi)
    },
    init_dpars = list(mu = function(y, aterms) mean(y),
                      sigma = function(y, aterms) sd(y)),
    post = list(mean_fn = function(dpars, aterms) dpars$mu))
  cat("user family package attr:",
      deparse(attr(ufam, "frmtmb_package", exact = TRUE)),
      "; frm_family_package:", frm_family_package(ufam), "\n")
  dd$y <- 1 + 0.5 * dd$x + rnorm(200)
  fu <- frm(y ~ x, family = ufam, data = dd)
  saveRDS(list(fc = fc, fr = fr, fm = fm, fu = fu, dd = dd),
          file.path(O, "saved2-fits.rds"))
}
s <- readRDS(file.path(O, "saved2-fits.rds"))
dd <- s$dd
nd <- data.frame(x = c(-1, 0, 1), g = factor(1:2)[c(1, 2, 1)],
                 time = c(1, 1, 1), out = 0L, y = 0)
cat("survival ns loaded after readRDS:", isNamespaceLoaded("frmtmb.survival"),
    "\n")
show <- function(lab, expr) {
  r <- tryCatch({
    v <- suppressWarnings(suppressMessages(expr))
    paste0("OK ", paste(class(v), collapse = "/"))
  }, error = function(e) {
    paste0(if (inherits(e, "frmtmb_extension_missing")) "REFUSED[ext] "
           else "ERROR ", substr(conditionMessage(e), 1, 110))
  })
  cat(sprintf("  %-26s %s\n", lab, r))
}
methods <- function(f, nm, rp = FALSE) {
  cat("==", nm, "\n")
  show("print", capture.output(print(f)))
  show("summary", summary(f))
  show("logLik", logLik(f))
  show("AIC", AIC(f))
  show("fixef", fixef(f))
  show("vcov", vcov(f))
  show("confint", confint(f))
  show("family", family(f))
  show("nobs", nobs(f))
  show("frm_linpred(link)", frm_linpred(f, type = "link"))
  show("frm_linpred(newdata)", frm_linpred(f, newdata = nd, type = "link"))
  show("frm_linpred(newdata,se)", frm_linpred(f, newdata = nd, type = "link",
                                              se.fit = TRUE))
  show("predict()", predict(f))
  show("predict(newdata=)", predict(f, newdata = nd))
  show("fitted", fitted(f))
  show("residuals", residuals(f))
  show("log_lik", log_lik(f))
  show("simulate", simulate(f, nsim = 1))
  show("conditional_effects", conditional_effects(f, effects = "x"))
  show("anova(f)", anova(f))
  show("anova(f, f)", anova(f, f))
  show("drop1", drop1(f))
  show("update(.~1)", update(f, . ~ 1, data = dd))
  show("frm_curve", frm_curve(f, newdata = data.frame(
    x = c(-1, 0, 1), time = 1), simultaneous = FALSE))
}
methods(s$fc, "cox")
methods(s$fr, "royston_parmar")
if (!is.character(s$fm)) methods(s$fm, "mixture(lognormal, royston_parmar)")
methods(s$fu, "user family (global env)")
