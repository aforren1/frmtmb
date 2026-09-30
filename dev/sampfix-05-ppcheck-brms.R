# Lane sampfix, script 05: pp_check()'s loo_* types against brms 2.23.0.
#
#   Rscript dev/sampfix-05-ppcheck-brms.R
#
# One compiled brms fit of y ~ x + (1 | g) (data as in
# dev/arcovsample-rev-07-ppcheck.R, seed 31; sampler seed 31, 2 chains of
# 500 post-warmup draws). A frmtmb draws object is built from brms's OWN
# draws, with brms's stanfit, so both pp_check()s see the same parameter
# draws and the same chains. bayesplot's ppc_loo_* functions are traced
# to capture what each pp_check() hands them: y, lw, psis_object. The
# pictures are not compared.

LANE <- "C:/Users/adf44/source/r/wt-sampfix-lib"
.libPaths(c(LANE, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
suppressMessages(library(brms))
cat("brms ", format(packageVersion("brms")), " bayesplot ",
    format(packageVersion("bayesplot")), " loo ",
    format(packageVersion("loo")), "\n", sep = "")

set.seed(31L)
dd <- data.frame(g = factor(rep(1:5, each = 6L)), t = rep(1:6, 5L))
n <- nrow(dd)
dd$x <- rnorm(n)
dd$y <- 0.6 + 0.5 * dd$x + rnorm(n, 0, 0.8)

rds <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-log/05-brmsfit.rds"
bfit <- if (file.exists(rds)) readRDS(rds) else {
  b <- brm(y ~ x + (1 | g), data = dd, family = gaussian(),
           chains = 2, iter = 1000, seed = 31, refresh = 0,
           backend = "rstan")
  saveRDS(b, rds)
  b
}

fit <- frm(frmtmb::bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
lab <- c(frmtmb::brms_par_labels(fit), "lp__")
bm <- as.matrix(bfit)
cat("frmtmb labels: ", paste(lab, collapse = " "), "\n", sep = "")
M <- matrix(NA_real_, nrow(bm), length(lab), dimnames = list(NULL, lab))
for (nm in intersect(lab, colnames(bm))) M[, nm] <- bm[, nm]
# the one covariance parameter is stored as the log sd; log_lik() does not
# read it, VarCorr() does, and the identity is checked below
M[, "theta_1"] <- log(bm[, "sd_g__Intercept"])
stopifnot(!anyNA(M))
ds <- structure(list(stanfit = bfit$fit, draws = M, fit = fit),
                class = "frmtmb_draws")
cat("nchains: brms ", nchains(bfit), ", frmtmb ", nchains(ds), "\n", sep = "")
cat("VarCorr sd identity, max abs diff: ",
    format(max(abs(VarCorr(ds, summary = FALSE)$g$sd[, 1] -
                     bm[, "sd_g__Intercept"]))), "\n", sep = "")
llb <- log_lik(bfit)
llf <- log_lik(ds)
cat("log_lik max abs diff: ", format(max(abs(llb - llf))),
    " (scale ", format(max(abs(llb))), ")\n", sep = "")

captured <- new.env()
types <- c("loo_pit_overlay", "loo_pit_qq", "loo_pit", "loo_intervals",
           "loo_ribbon")
for (ty in types) {
  fn <- paste0("ppc_", ty)
  have <- intersect(c("y", "lw", "psis_object"),
                    names(formals(get(fn, asNamespace("bayesplot")))))
  body <- lapply(have, function(a) {
    bquote(if (!missing(.(as.name(a)))) .(as.name(a)))
  })
  names(body) <- have
  trace(fn, where = asNamespace("bayesplot"), print = FALSE,
        tracer = bquote({
          captured[[.(fn)]] <- .(as.call(c(as.name("list"), body)))
        }))
}
grab <- function(obj, ty, ...) {
  rm(list = ls(captured), envir = captured)
  pdf(NULL)
  on.exit(dev.off())
  r <- tryCatch(suppressMessages(suppressWarnings(
    pp_check(obj, type = ty, ...))), error = function(e) e)
  if (inherits(r, "error")) {
    cat("    pp_check error: ", conditionMessage(r), "\n", sep = "")
  }
  captured[[paste0("ppc_", ty)]]
}
lwmat <- function(cap) {
  if (!is.null(cap$lw)) cap$lw else weights(cap$psis_object, log = TRUE)
}
for (ids in list(NULL, seq(1L, 1000L, by = 7L))) {
  cat("\n-- draws: ", if (is.null(ids)) "all 1000" else
    paste0(length(ids), " (draw_ids = seq(1, 1000, by = 7))"), "\n",
    sep = "")
  for (ty in types) {
    cb <- grab(bfit, ty, draw_ids = ids)
    cf <- grab(ds, ty, draw_ids = ids)
    if (is.null(cb) || is.null(cf)) {
      cat(sprintf("  %-16s brms %s, frmtmb %s\n", ty,
                  if (is.null(cb)) "NOT CALLED" else "called",
                  if (is.null(cf)) "NOT CALLED" else "called"))
      next
    }
    which_arg <- if (!is.null(cb$lw)) "lw" else "psis_object"
    same_arg <- identical(is.null(cb$lw), is.null(cf$lw))
    lb <- lwmat(cb)
    lf <- lwmat(cf)
    kb <- if (!is.null(cb$psis_object)) loo::pareto_k_values(cb$psis_object)
    kf <- if (!is.null(cf$psis_object)) loo::pareto_k_values(cf$psis_object)
    cat(sprintf(paste0("  %-16s arg %-11s same arg %s  dim %s  ",
                       "y identical %s  max |lw diff| %.3g%s\n"),
                ty, which_arg, same_arg, paste(dim(lf), collapse = "x"),
                identical(as.numeric(cb$y), as.numeric(cf$y)),
                max(abs(lb - lf)),
                if (!is.null(kb)) sprintf("  max |k diff| %.3g",
                                          max(abs(kb - kf))) else ""))
  }
}
cat("\nDONE\n")
