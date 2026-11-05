# Review of lane gpgap, item 1: the five pairs at the marginal maximum, read
# independently: brms's nugget in a covariance conditioned by sigma^2, at
# frm()'s estimate, and each model's own maximum. For student and nl there is
# no closed form, so the non-centered tape at brms's line stands in.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-pairs.R pairs.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
requireNamespace("brms")
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
ns <- asNamespace("frmtmb")
orig <- get("gp_corr", envir = ns)
put <- function(val) {
  unlockBinding("gp_corr", ns); assign("gp_corr", val, envir = ns)
  lockBinding("gp_corr", ns)
}
abs_corr <- function(theta, blk) {
  Q <- blk[["aux_D2"]][[1]] / (2 * exp(2 * theta[2]))
  exp(-Q) + diag(1, nrow(Q)) * (1e-12 * exp(-2 * theta[1]))
}
keys <- c(p022 = "gp_pred() x mvbf", p023 = "gp_pred() x prior", p094 = "lognormal x gp",
          p109 = "student x gp", p136 = "gp x nl")
for (id in names(keys)) {
  spec <- E$eq_specs[[keys[[id]]]]()
  spec$prior <- NULL
  fit <- suppressWarnings(E$eq_frm_fit(spec))
  bk <- Filter(function(b) identical(b$covstruct, "gp"), fit$frame$re_blocks)[[1]]
  th <- fit$estimates$theta[bk$theta_idx]
  fam <- spec$family %||% "gaussian()"
  ll <- -fit$opt$objective
  warned <- !is.null(gp_nugget_message <- get("gp_nugget_message", ns)(fit))
  cat(sprintf("PAIR %s %s: frm ll %.6f sdgp %.5f rho %.5f; L4 tolerance 1e-4*|ll| %.2e; nugget warning %s\n",
              id, keys[[id]], ll, exp(th[[1]]), exp(th[[2]]), 1e-4 * abs(ll), warned))
  if (!grepl("student", fam) && !grepl("nl = TRUE", spec$bf)) {
    d <- spec$data
    y <- d[[if (grepl("y1", spec$bf)) "y1" else "y"]]
    if (grepl("lognormal", fam)) y <- log(y)
    u <- sort(unique(d$x)); J <- match(d$x, u); D2 <- outer(u, u, "-")^2
    nll <- function(p, nug) {
      s2 <- exp(2 * p[2])
      K <- s2 * exp(-D2 / (2 * exp(2 * p[3]))); diag(K) <- diag(K) + nug(s2)
      V <- K[J, J]; diag(V) <- diag(V) + exp(2 * p[4])
      R <- chol(V)
      r <- backsolve(R, y - p[1], transpose = TRUE)
      0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * length(y) * log(2 * pi)
    }
    st <- c(mean(y), th[[1]], th[[2]], -1)
    a <- nlminb(st, nll, nug = function(s2) 1e-6 * s2,
                control = list(iter.max = 800, eval.max = 1200, rel.tol = 1e-12))
    b <- nlminb(st, nll, nug = function(s2) 1e-12,
                control = list(iter.max = 800, eval.max = 1200, rel.tol = 1e-12))
    cat(sprintf(paste0("     closed form: frm nugget ll %.8f sigma %.6f sdgp %.6f | brms nugget ll %.8f ",
                       "sigma %.6f sdgp %.6f | maxima differ %.3e | brms model at frm max loses %.3e\n"),
                -a$objective, exp(a$par[4]), exp(a$par[2]), -b$objective, exp(b$par[4]), exp(b$par[2]),
                a$objective - b$objective, nll(a$par, function(s2) 1e-12) - b$objective))
  }
  frn <- fit$frame
  frn$ncp_blocks <- which(vapply(frn$re_blocks, function(b) identical(b$covstruct, "gp"), NA))
  put(abs_corr)
  r <- tryCatch({
    fn <- frmtmb::build_objective(frn)
    tpl <- frn$par_template
    tpl$b[] <- 0
    o <- RTMB::MakeADFun(fn, tpl, random = "b", map = frn$map, silent = TRUE)
    op <- suppressWarnings(nlminb(o$par, o$fn, o$gr, control = list(iter.max = 500, eval.max = 800)))
    h <- op$par[names(op$par) == "theta"][bk$theta_idx]
    sprintf("conv %d ll %.6f (gap to frm %.3e) sdgp %.5f rho %.5f", op$convergence, -op$objective,
            -op$objective - ll, exp(h[[1]]), exp(h[[2]]))
  }, error = function(e) paste("ERROR", conditionMessage(e)))
  put(orig)
  cat("     ncp at brms's line:", r, "\n")
}
