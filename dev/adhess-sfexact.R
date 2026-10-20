# Lane adhess: the battery rows where the lane's largest SE error is
# higher than the base's (dev/adhess-log/sf-ref/sum.txt). The battery's
# reference there is glm()'s or lmer()'s vcov(), which is not exact.
# This reruns the battery's "lin" and "poly" groups (same data, same
# seeds: the loop of dev/adhess-sf-rev3-falarm.R, evaluated with its
# own code) and records, per fit, the largest |SE / ref - 1| against
#   batt   the battery's own reference (reproduces its max_rel)
#   exact  without random effects: the inverse of X' W X formed and
#          inverted in 256 bits at frmtmb's own estimate; with random
#          effects: a Richardson Jacobian of the AD gradient in SE units
#          (numDeriv, six halvings, base step 0.01, as
#          dev/adhess-restep.R), inverted in double
#   Rscript dev/adhess-sfexact.R <lane|base> <group> <out.rds> [seeds]
a <- commandArgs(TRUE)
Sys.setenv(ADHESS_ARM = a[1])
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
grp <- a[2]
out <- a[3]
seeds <- if (length(a) > 3) eval(parse(text = a[4])) else 1:20
# the battery's loop assigns `m` (a glm) in the global environment, so
# the 256-bit helpers get an environment of their own
hold <- new.env(parent = globalenv())
for (nm in c("m", "mM", "num", "col1", "mp_inv", "glm_info_exact",
             "glm_cov_exact")) {
  f <- get(nm)
  if (is.function(f)) environment(f) <- hold
  assign(nm, f, envir = hold)
}
glm_cov_exact <- hold$glm_cov_exact
se_outer <- function(fit) {
  v <- suppressWarnings(vcov(fit, full = TRUE))
  stats::setNames(sqrt(diag(as.matrix(v))), ns$outer_par_names(fit))
}
ref_hessian <- function(obj, p, D, d0) {
  g0 <- function(dz) as.numeric(obj$gr(p + dz / D)) / D
  J <- numDeriv::jacobian(g0, rep(0, length(p)), method = "Richardson",
                          method.args = list(eps = d0, d = d0, zero.tol = 1,
                                             r = 6, v = 2))
  J <- (J + t(J)) / 2
  J * outer(D, D)
}
lmer_ref <- function(fo, d) {
  m <- suppressWarnings(suppressMessages(lme4::lmer(fo, data = d,
                                                    REML = FALSE)))
  sqrt(diag(as.matrix(stats::vcov(m))))
}
rows <- list()
rec <- function(model, seed, expr, ref = NULL, coef_names = NULL) {
  d <- get("d", parent.frame())
  r <- cap(expr)
  fit <- r$v
  if (is_err(fit)) return(invisible())
  se <- se_outer(fit)
  on <- names(se)
  nmz <- if (is.null(coef_names)) on else coef_names
  s1 <- se[match(nmz, on)]
  batt <- if (is.null(ref)) NA else
    max(abs(s1 / ref[seq_along(nmz)] - 1))
  re <- length(fit$obj$env$random) > 0
  ex <- NA
  if (!re && grepl("^pois_cor|^logit_unc", model)) {
    fam <- if (grepl("^pois", model)) "poisson" else "bernoulli"
    X <- if (fam == "poisson") cbind(1, d$x1, d$x2) else cbind(1, d$x)
    b <- unname(fit$opt$par[seq_len(ncol(X))])
    s0 <- sqrt(diag(glm_cov_exact(X, b, fam)))
    ex <- max(abs(s1 / s0 - 1))
  } else if (re) {
    obj <- fit$obj
    p <- ns$sdr_outer_point(obj)
    st <- ns$obj_state_save(obj)
    H0 <- ns$fd_hessian(p, function(q) obj$gr(q))$H
    R2 <- ref_hessian(obj, p, sqrt(abs(diag(H0))), 1e-2)
    ns$obj_state_restore(obj, st)
    s0 <- sqrt(diag(solve(R2)))[match(nmz, on)]
    ex <- max(abs(s1 / s0 - 1))
  } else if (!re) {
    he <- fit$obj$he(fit$opt$par)
    s0 <- sqrt(diag(solve((he + t(he)) / 2)))[match(nmz, on)]
    ex <- max(abs(s1 / s0 - 1))
  }
  rows[[length(rows) + 1L]] <<- data.frame(model = model, seed = seed,
    batt = batt, exact = ex, warn = wtag(r$w), lost = lost_of(fit))
}
# the battery's loop, evaluated as it is written there
src <- parse("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-sf-rev3-falarm.R")
loop <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("for")),
               as.list(src))
stopifnot(length(loop) == 1L)
eval(loop[[1]])
res <- do.call(rbind, rows)
saveRDS(res, out)
cat("fits", nrow(res), "\n")
