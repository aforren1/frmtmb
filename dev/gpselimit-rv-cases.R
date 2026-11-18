# Lane gpselimit, fix round 1: a copy of dev/reviews/gpselimit-review-cases.R
# whose output goes to dev/gpselimit-out/rv/ (the review's logs stay).
# Review of lane gpselimit, checks 2, 4 and 5: cases built to break the
# white-limit certificate and the message, run on both builds.
#  A  plateau: resolved sweep rows started inside the white plateau
#  B  local: the 13 white sweep rows refitted from a resolved start
#  C  priors: weak lognormal on lscale, brms's default-like prior on sdgp
#  D  2-D fields, iso and aniso, exp_quad, matern32, matern52
#  E  gp(by = f), F  sigma ~ gp(), G  two gp() terms, H  REML, I  poisson
# Per fit: conditions, the SE of every outer parameter, lost reasons and
# kinds, each hsgp block's ls/spacing and the objective along its ridges
# (t > 0 toward lscale 0), and for each hypothesis() string the message
# names: its estimate and SE beside a numerical Hessian's SE of the same
# quantity in (psi, ridge) coordinates.
# Usage: Rscript gpselimit-review-cases.R <lane|base> [groups]
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-rv-lib.R")
a <- commandArgs(TRUE)
arm <- a[1]
groups <- if (length(a) >= 2) strsplit(a[2], "")[[1]] else LETTERS[1:9]
rv_libs(arm)
ns <- asNamespace("frmtmb")
ts <- c(-6, -4, -2, -1, 1, 2, 12)
res <- list()
psi_se <- function(fit, blocks) {
  om <- ns$outer_par_map(fit)
  thp <- which(om$comp == "theta")
  p <- fit$opt$par
  H <- stats::optimHess(p, fit$obj$fn, fit$obj$gr)
  fit$obj$fn(p)
  A <- diag(length(p))
  drop <- integer(0); psi <- integer(0)
  for (bk in blocks) {
    pos <- thp[bk$theta_idx]
    dirs <- rv_dirs(bk)
    A[, pos] <- 0
    A[pos[1], pos[1]] <- 0.5
    for (j in seq_along(dirs)) A[pos, pos[1 + j]] <- dirs[[j]]
    psi <- c(psi, pos[1]); drop <- c(drop, pos[-1])
  }
  Hn <- t(A) %*% H %*% A
  keep <- setdiff(seq_along(p), drop)
  V <- tryCatch(solve(Hn[keep, keep]), error = function(e) NULL)
  if (is.null(V)) return(rep(NA, length(psi)))
  sqrt(diag(V)[match(psi, keep)])
}
report <- function(lab, q) {
  fit <- q$fit
  if (inherits(fit, "error")) {
    cat("CASE", lab, arm, "ERROR", conditionMessage(fit), "\n")
    return(invisible(NULL))
  }
  ss <- rv_ses(fit)
  kind <- fit$cache$se_boundary_kind
  cat(sprintf(paste0("CASE %s %s | ll %.6f conv %d | se %s | lost [%s] | ",
                     "kind [%s] | conds [%s]\n"),
              lab, arm, as.numeric(logLik(fit)), fit$opt$convergence,
              paste(paste0(names(ss$se), "=", format(ss$se, digits = 4)),
                    collapse = ","),
              paste(paste0(names(ss$lost), "=", ss$lost), collapse = ","),
              paste(paste0(names(kind), "=", kind), collapse = ","),
              rv_cstr(q$conds)))
  white <- list()
  for (bk in rv_hsgp(fit)) {
    th <- fit$estimates$theta[bk$theta_idx]
    w <- ns$gp_basis_white(bk, th)
    cert <- if (exists("gp_white_ridge", ns)) {
      tryCatch(ns$gp_white_ridge(fit, bk), error = function(e) NA)
    } else NA
    rg <- rv_ridge2(fit, bk, ts)
    cat(sprintf("   block %s | sd %.4g ls/sp %s | white %d cert %s | %s\n",
                bk$term_label %||% "?", exp(th[1]),
                paste(format(exp(th[-1]) / bk$gp_spacing, digits = 3),
                      collapse = ","), w, cert,
                paste(vapply(rg, function(r) {
                  paste0("ridge(t=", paste(ts, collapse = ","), ") ",
                         paste(format(r, digits = 3), collapse = ","))
                }, ""), collapse = " ; ")))
    if (w) white[[length(white) + 1L]] <- bk
  }
  for (c in q$conds) {
    if (!"frmtmb_gp_basis" %in% c$class) next
    cat("   text:", gsub("\n", " ", c$text), "\n")
  }
  m <- Filter(function(c) "frmtmb_gp_basis" %in% c$class &&
                c$type == "message", q$conds)
  if (length(m) && length(white)) {
    hs <- gsub("\"", "", regmatches(m[[1]]$text,
                                    gregexpr("\"[^\"]*\"", m[[1]]$text))[[1]])
    hse <- psi_se(fit, white)
    for (i in seq_along(hs)) {
      r <- tryCatch(rv_q(hypothesis(fit, hs[i], class = ""))$hypothesis,
                    error = function(e) conditionMessage(e))
      if (is.character(r)) {
        cat("   hyp", hs[i], "ERROR", r, "\n")
      } else {
        cat(sprintf("   hyp %s | est %.6f se %.6f | hessian psi se %.6f\n",
                    hs[i], r$Estimate, r$Est.Error, hse[i]))
      }
    }
  }
  res[[lab]] <<- list(est = fit$estimates, par = fit$opt$par,
                      ll = as.numeric(logLik(fit)), se = ss$se,
                      lost = ss$lost, conds = q$conds)
  invisible(fit)
}
gfit <- function(...) rv_fit(frm(...))

if ("A" %in% groups) {
  # resolved rows started inside the white plateau: the ridge toward 0 is
  # level there, while the resolved maximum is higher
  for (row in c(6, 12, 18, 24, 30)) {
    s <- rv_sweep_row(row)
    f0 <- report(sprintf("A%02d-default", row),
                 gfit(bf(y ~ gp(x, k = s$k)), data = s$d,
                      family = gaussian()))
    sp <- rv_hsgp(f0)[[1]]$gp_spacing
    for (fr in c(3e-5, 1e-3)) {
      report(sprintf("A%02d-start%.0e", row, fr),
             gfit(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian(),
                  start = list(theta = c(log(30), log(fr * sp)))))
    }
  }
}
if ("B" %in% groups) {
  for (row in c(1:5, 9, 10, 14, 21, 22, 25:27)) {
    s <- rv_sweep_row(row)
    f0 <- gs <- rv_q(frm(bf(y ~ gp(x, k = s$k)), data = s$d,
                         family = gaussian()))
    sp <- rv_hsgp(f0)[[1]]$gp_spacing
    l0 <- as.numeric(logLik(f0))
    for (fr in c(0.3, 1, 3)) {
      q <- gfit(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian(),
                start = list(theta = c(log(1), log(fr * sp))))
      if (inherits(q$fit, "error")) next
      bk <- rv_hsgp(q$fit)[[1]]
      cat(sprintf("LOCAL B%02d %s start %.1f sp | white ll %.6f | refit ll %.6f diff %.3g ls/sp %.3g conds [%s]\n",
                  row, arm, fr, l0, as.numeric(logLik(q$fit)),
                  as.numeric(logLik(q$fit)) - l0,
                  exp(q$fit$estimates$theta[bk$theta_idx[2]]) / sp,
                  rv_cstr(q$conds)))
    }
  }
}
if ("C" %in% groups) {
  s <- rv_sweep_row(3)
  f0 <- rv_q(frm(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian()))
  sp <- rv_hsgp(f0)[[1]]$gp_spacing
  # gp_brms_values() reads lscale on brms's scaled inputs
  bv <- ns$gp_brms_values(f0, f0$estimates$theta)
  ratio <- bv[["lscale_gpx"]] / exp(f0$estimates$theta[2])
  for (sdl in c(1, 3)) {
    pr <- set_prior(sprintf("lognormal(%.10g, %g)", log(0.01 * sp * ratio),
                            sdl), class = "lscale")
    report(sprintf("C-lscale-lognormal-sd%g", sdl),
           gfit(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian(),
                prior = pr))
  }
  report("C-sdgp-student", gfit(bf(y ~ gp(x, k = s$k)), data = s$d,
                                family = gaussian(),
                                prior = set_prior("student_t(3, 0, 2.5)",
                                                  class = "sdgp")))
  report("C-sdgp-exp0.001", gfit(bf(y ~ gp(x, k = s$k)), data = s$d,
                                 family = gaussian(),
                                 prior = set_prior("exponential(0.001)",
                                                   class = "sdgp")))
}
if ("D" %in% groups) {
  set.seed(5150)
  n <- 150
  d2 <- data.frame(x1 = runif(n), x2 = runif(n), y = rnorm(n))
  for (cv in c("exp_quad", "matern32", "matern52")) {
    for (iso in c(TRUE, FALSE)) {
      form <- bf(y ~ gp(x1, x2, k = 4, iso = iso, cov = cv))
      f0 <- rv_q(frm(form, data = d2, family = gaussian()))
      bk <- rv_hsgp(f0)[[1]]
      lp <- f0$frame$linpreds[[ns$gp_block_lp(bk, f0$frame)]]
      Z <- as.matrix(ns$lp_Z(f0, lp))[, bk$c_idx, drop = FALSE]
      set.seed(5161)
      dd <- d2
      dd$y <- as.vector(Z %*% rnorm(ncol(Z))) + 0.1 * rnorm(n)
      report(sprintf("D-%s-iso%d", cv, iso),
             gfit(form, data = dd, family = gaussian()))
    }
  }
  # anisotropic, one dimension white only: x2 carries a smooth signal
  set.seed(5170)
  dd <- d2
  dd$y <- sin(3 * d2$x2) + rnorm(n, 0, 0.3)
  report("D-aniso-onewhite",
         gfit(bf(y ~ gp(x1, x2, k = 4, iso = FALSE)), data = dd,
              family = gaussian()))
}
if ("E" %in% groups) {
  set.seed(9)
  d <- data.frame(x = stats::runif(60), f = factor(rep(c("a", "b"), 30)))
  d$y <- sin(3 * d$x) + stats::rnorm(60, 0, 0.3)
  report("E-by-f-k5", gfit(y ~ gp(x, by = f, k = 5), data = d))
  for (row in c(3, 21)) {
    s <- rv_sweep_row(row)
    s$d$f <- factor(rep(c("a", "b"), length.out = nrow(s$d)))
    report(sprintf("E-by-f-row%d", row),
           gfit(y ~ gp(x, by = f, k = s$k), data = s$d))
  }
}
if ("F" %in% groups) {
  for (row in c(3, 21)) {
    s <- rv_sweep_row(row)
    set.seed(700 + row)
    s$d$ys <- rnorm(nrow(s$d), 0, exp(0.5 * s$d$y))
    for (k in c(5, 10)) {
      report(sprintf("F-sigma-row%d-k%d", row, k),
             gfit(bf(ys ~ 1, sigma ~ gp(x, k = k)), data = s$d,
                  family = gaussian()))
    }
  }
}
if ("G" %in% groups) {
  s3 <- rv_sweep_row(3)
  s6 <- rv_sweep_row(21)
  set.seed(4401)
  d <- data.frame(x1 = s3$d$x, x2 = sample(s6$d$x, 60),
                  y = s3$d$y + rnorm(60, 0, 0.1))
  d$y <- d$y + sin(d$x2)
  report("G-two-k10-k30", gfit(bf(y ~ gp(x1, k = 10) + gp(x2, k = 30)),
                               data = d, family = gaussian()))
  report("G-two-k10-k10", gfit(bf(y ~ gp(x1, k = 10) + gp(x2, k = 10)),
                               data = d, family = gaussian()))
  d$y2 <- s3$d$y + rnorm(60, 0, 0.2)
  report("G-two-white", gfit(bf(y2 ~ gp(x1, k = 10) + gp(x1, k = 5)),
                             data = d, family = gaussian()))
}
if ("H" %in% groups) {
  for (row in c(3, 9, 6)) {
    s <- rv_sweep_row(row)
    report(sprintf("H-REML-row%d", row),
           gfit(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian(),
                REML = TRUE))
  }
}
if ("I" %in% groups) {
  for (row in c(3, 9, 21)) {
    s <- rv_sweep_row(row)
    set.seed(800 + row)
    s$d$cnt <- rpois(nrow(s$d), exp(1 + 0.5 * s$d$y))
    report(sprintf("I-poisson-row%d", row),
           gfit(bf(cnt ~ gp(x, k = s$k)), data = s$d, family = poisson()))
  }
}
saveRDS(res, file.path(RV_OUT, sprintf("cases-%s-%s.rds", arm,
                                       paste(groups, collapse = ""))))
