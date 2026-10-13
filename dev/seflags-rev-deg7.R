# Reviewer of lane seflags, claim 4: is the degree-7 raw polynomial with
# (1 | g) (seed 7, test-se-flags.R) identified, and do the lane's NaN
# standard errors come from the fit not being at the optimum? lm(),
# lmer() (ML and REML), frmtmb default, frmtmb refit from lmer's
# estimates, frmtmb REML (coefficient SEs from the joint precision).
#   Rscript dev/seflags-rev-deg7.R <lane|base> [deg] [seed]
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev-common.R")
a <- commandArgs(trailingOnly = TRUE)
deg <- if (length(a) > 1) as.integer(a[2]) else 7L
sd0 <- if (length(a) > 2) as.integer(a[3]) else 7L
set.seed(sd0)
d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
  stats::rnorm(200, 0, 0.2)
for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
rhs <- paste0("x", 1:deg, collapse = " + ")
fo <- stats::as.formula(paste("y ~", rhs, "+ (1 | g)"))
cn <- c("(Intercept)", paste0("x", 1:deg))
X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), d)
Xs <- sweep(X, 2, sqrt(colSums(X^2)), "/")
sv <- svd(Xs)$d
cat(sprintf("degree %d seed %d | scaled X: kappa %.3g, smallest sv^2 / largest %.3g\n",
            deg, sd0, max(sv) / min(sv), (min(sv) / max(sv))^2))
lmf <- stats::lm(stats::as.formula(paste("y ~", rhs)), data = d)
cat("lm SE:   ", format(sqrt(diag(vcov(lmf))), digits = 4), "\n")
m <- allow <- suppressWarnings(suppressMessages(
  lme4::lmer(fo, data = d, REML = FALSE,
             control = lme4::lmerControl(check.conv.singular = "ignore"))))
m2 <- suppressWarnings(lme4::lmer(fo, data = d, REML = FALSE,
  control = lme4::lmerControl(optimizer = "bobyqa",
                              optCtrl = list(maxfun = 1e5))))
ref <- sqrt(diag(as.matrix(vcov(m))))
cat(sprintf("lmer ML logLik %.6f (bobyqa %.6f) | SE: %s\n",
            as.numeric(logLik(m)), as.numeric(logLik(m2)),
            paste(format(ref, digits = 4), collapse = " ")))
cat("lmer bobyqa SE:", format(sqrt(diag(as.matrix(vcov(m2)))), digits = 4),
    "\n")
mr <- suppressWarnings(lme4::lmer(fo, data = d, REML = TRUE))
cat("lmer REML SE:", format(sqrt(diag(as.matrix(vcov(mr)))), digits = 4),
    "\n")
show <- function(lab, r) {
  f <- r$v
  if (inherits(f, "rev_err")) { cat(lab, "ERROR", f, "\n"); return() }
  se <- se_outer(f)
  s <- se[match(cn, names(se))]
  cat(sprintf("%-22s code %d logLik %.6f | SE/lmer %s\n   lost %s | warnings %d\n",
              lab, f$opt$convergence, -f$opt$objective,
              paste(format(s / ref, digits = 3), collapse = " "), lost_of(f),
              length(r$w)))
  for (w in r$w) cat("   W:", substr(w, 1, 230), "\n")
  invisible(f)
}
f0 <- show("frmtmb default", cap(frm(fo, data = d)))
tm <- par_template(f0)
cat("template slots:", paste(names(tm), collapse = " "), "\n")
st <- list(beta = unname(lme4::fixef(m2)),
           theta = log(attr(lme4::VarCorr(m2)$g, "stddev")[[1]]),
           betad = log(stats::sigma(m2)))
r1 <- cap(frm(fo, data = d, start = st))
f1 <- show("frmtmb from lmer est", r1)
# the objective at lmer's estimate, in frmtmb's own parameters
p <- f0$opt$par
on <- ns$outer_par_names(f0)
q <- p
q[match(cn, on)] <- lme4::fixef(m2)
q[on == "theta_1"] <- st$theta
q[grepl("^sigma", on)] <- st$betad
cat(sprintf("frmtmb objective at its optimum %.6f, at lmer's %.6f\n",
            f0$obj$fn(p), f0$obj$fn(q)))
f0$obj$fn(p)
# SEs at lmer's point from the frmtmb objective: FD Hessian of the AD
# gradient there, then a plain solve()
h <- ns$fd_hessian(q, f0$obj$gr)
f0$obj$fn(p)
V <- tryCatch(solve(h$H), error = function(e) NULL)
if (!is.null(V)) {
  cat("solve(FD Hessian at lmer's point) SE/lmer:",
      format(sqrt(diag(V))[match(cn, on)] / ref, digits = 3), "\n")
}
show("frmtmb REML", cap(frm(fo, data = d, REML = TRUE)))
rr <- sqrt(diag(as.matrix(vcov(mr))))
fr <- suppressWarnings(frm(fo, data = d, REML = TRUE))
sr <- suppressWarnings(fixef(fr)[, "Est.Error"])
cat("frmtmb REML fixef SE / lmer REML:", format(sr / rr, digits = 4), "\n")
