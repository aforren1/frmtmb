# Short check of lane drawdir's fix round 1 (option (b)), part 1, lane
# library. Checks 1 to 4 of the short brief:
#   A  FIX 1 on slpw:12 and on constructed objectives: no Inf cached or
#      printed; a direction read back at t is judged on t's scale.
#   B  option (b): the components of every firing direction on the 12
#      fires (review's 8, lane's 4) and on test-gr-by.R's fit; "all"
#      against "any".
#   C  frmtmb_draw_wide_tail: how many per call, for each entry point.
#   D  the new-level caveat: with and without a new level.
#   E  slpw:12's reweighted arm (a non-variance fire) at 20000 draws,
#      seed 1, against the review's first-build line.
# Usage: Rscript dev/reviews/drawdir-short-1.R [A B C D E]
args <- commandArgs(TRUE)
if (!length(args)) args <- c("A", "B", "C", "D", "E")
.libPaths(c("C:/Users/adf44/source/r/wt-drawdir-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), "\n")
ex <- parse("dev/reviews/drawdir-review-survey.R")
keep <- vapply(ex, function(e) {
  is.call(e) && identical(e[[1]], as.name("<-")) &&
    as.character(e[[2]]) %in% c("fq", "ilogit", "cum_mean", "cum_draw",
                                "design", "seeds_of")
}, TRUE)
for (e in ex[keep]) eval(e, globalenv())
`%||%` <- function(a, b) if (is.null(a)) b else a

# S designs (lane's 4), written from the test's generator
s_fit <- function(kind, seed) {
  sds <- as.numeric(substring(kind, 2))
  set.seed(2000 + seed)
  d <- data.frame(g = factor(rep(1:15, each = 8)), x = stats::rnorm(120))
  ug <- stats::rnorm(15, 0, 0.8)
  sg <- stats::rnorm(15, 0, sds)
  d$y <- 1 + 0.5 * d$x + ug[d$g] + sg[d$g] * d$x + stats::rnorm(120)
  fq(frm(bf(y ~ x + (1 + x | g)), data = d))
}
fit_of <- function(id) {
  p <- strsplit(id, ":")[[1]]
  if (p[1] %in% c("S0.1", "S0.2")) return(s_fit(p[1], as.integer(p[2])))
  design(p[1], as.integer(p[2]))$fit
}

# every direction's stat, from the installed function with its last
# expression replaced (the reading itself is unchanged)
dirs_all <- local({
  f <- ns$draw_check_dirs
  b <- as.list(body(f))
  at <- which(vapply(b, function(e) {
    identical(deparse(e)[1], "i <- which.max(stat)")
  }, TRUE))
  body(f) <- as.call(c(b[seq_len(at - 1L)],
                       quote(list(stat = stat, R = R, k = k, at = at,
                                  eA = eA, n = n))))
  f
})
cap <- new.env()
grab_dw <- function(fit) {
  real <- ns$draw_reweight
  assignInNamespace("draw_reweight", function(object, dw) {
    dw$off <- length(dw$v0) - length(object$opt$par)
    cap$dw <- dw
    dw$D
  }, ns = "frmtmb")
  on.exit(assignInNamespace("draw_reweight", real, ns = "frmtmb"))
  invisible(fq(ns$predict_par_drawer(fit, TRUE, 2L)))
  cap$dw
}
# conditions said by an expression, by class
said <- function(expr) {
  out <- list()
  val <- withCallingHandlers(expr,
    warning = function(w) {
      out[[length(out) + 1L]] <<- list(cls = class(w)[1],
                                       msg = conditionMessage(w))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      out[[length(out) + 1L]] <<- list(cls = class(m)[1],
                                       msg = conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  list(val = val, said = out)
}
n_cls <- function(s, cls) sum(vapply(s$said, function(x) x$cls == cls, NA))

if ("A" %in% args) {
  cat("\n#### A. FIX 1\n")
  fit <- fit_of("slpw:12")
  dw <- grab_dw(fit)
  lab <- dw$labels
  comp <- dw$comps
  # independent reading of the formerly non-finite ends: the review's
  # wall script's way, with t scaled deviates through the side law
  off <- dw$off
  np <- length(fit$opt$par)
  pos <- off + seq_len(np)
  Lf <- dw$law$Lc %||% dw$L
  sv <- colSums(Lf^2)
  R <- setdiff(pos[sv[pos] > 0], c(dw$held, dw$law$edge))
  k <- length(R)
  U <- Lf[, R, drop = FALSE] / rep(sqrt(sv[R]), each = nrow(Lf))
  Z <- cbind(2 * U, -2 * U)
  et <- function(t, cols) {
    Zt <- t * Z[, cols, drop = FALSE]
    Dt <- dw$v0 + if (is.null(dw$law)) crossprod(dw$L, Zt) else
      ns$draw_edge_apply(dw$law, dw$L, Zt)
    if (length(dw$held)) Dt[dw$held, ] <- dw$v0[dw$held]
    ns$draw_logratio(fit, dw, Dt, Zt)$e
  }
  e1 <- et(1, seq_len(2 * k))
  da <- dirs_all(fit, dw)
  cat(sprintf("slpw:12: %d directions; non-finite at t = 1: %d\n",
              2 * k, sum(!is.finite(e1))))
  for (i in seq_len(2 * k)) {
    j <- (i - 1L) %% k + 1L
    nm <- sprintf("%s%s (%s)", if (i > k) "-" else "+", lab[R[j]],
                  comp[R[j]])
    line <- sprintf("  %-26s e(1) %9.4g", nm, e1[i])
    if (!is.finite(e1[i])) {
      e75 <- et(0.75, i)
      e50 <- et(0.5, i)
      line <- paste0(line, sprintf(
        " | e(0.75) %9.4g [thr 6*t^2 = %.3f; /t^2 %9.4g] e(0.5) %9.4g [thr %.3f; /t^2 %9.4g]",
        e75, 6 * 0.75^2, e75 / 0.75^2, e50, 6 * 0.25, e50 / 0.25))
    }
    line <- paste0(line, sprintf(" | lane at %.2f eA %9.4g stat %9.4g",
                                 da$at[i], da$eA[i], da$stat[i]))
    cat(line, "\n")
  }
  cat("any stat infinite:", any(is.infinite(da$stat) & da$stat > 0), "\n")
  set.seed(1)
  s <- said(predict(fit, ndraws = 200, summary = FALSE))
  chk <- fit$cache[["draw_check"]]
  cat(sprintf("cached: fire %s dir %s dir_par %s dir_n %d warn_only %s\n",
              chk$fire, format(chk$dir), chk$dir_par, chk$dir_n,
              chk$warn_only))
  cat("cached fields finite:", all(is.finite(unlist(chk[c("dir")]))), "\n")
  for (x in s$said) cat("  said", x$cls, ":", substr(x$msg, 1, 400), "\n")
  cat("any condition text contains 'Inf':",
      any(vapply(s$said, function(x) grepl("Inf", x$msg), NA)), "\n")

  # constructed: 2-D, correlation 0.9, a quadratic misfit a/2 (u'p)^2
  # along u = (1, 1)/sqrt(2) (no axis), and the objective NaN (or +Inf)
  # where |u'p| > 1.6: every direction's end (u'p = 2.687 at t = 1, 2.015
  # at 0.75) fails, t = 0.5 (1.34) reads. The misfit is quadratic, so
  # the excess at t over t^2 is the full one exactly; read unscaled
  # (e(0.5) = 0.25 x full = 3.42) it would not pass 6
  V <- matrix(c(1, 0.9, 0.9, 1), 2)
  Vi <- solve(V)
  u <- c(1, 1) / sqrt(2)
  a <- 4
  quad <- function(p) sum(p * (Vi %*% p)) / 2
  csd <- 1 / sqrt(diag(Vi))
  full_dir <- 2 * V[, 1] / sqrt(V[1, 1])
  want <- a / 2 * sum(u * full_dir)^2 - a / 2 * (u[1] * 2 * csd[1])^2
  dwc <- list(L = chol(V), law = NULL, v0 = c(0, 0), held = integer(0),
              off = 0L, labels = c("a", "b"), comps = c("beta", "theta"))
  one <- function(fn) {
    object <- list(obj = list(fn = fn, env = new.env()),
                   opt = list(par = c(0, 0)))
    ns$draw_check_dirs(object, dwc)
  }
  for (bad in c(NaN, Inf)) {
    r <- one(function(p) {
      if (abs(sum(u * p)) > 1.6) bad else quad(p) + a / 2 * sum(u * p)^2
    })
    cat(sprintf(paste0("constructed (end %s): stat %.6f want %.6f ",
                       "(unscaled e(0.5) %.4f) n %d par %s vc_only %s ",
                       "fires %s\n"),
                format(bad), r$stat, want, 0.25 * a / 2 *
                  sum(u * full_dir)^2 - 0.25 * a / 2 *
                  (u[1] * 2 * csd[1])^2,
                r$n, r$par, r$vc_only, r$stat > ns$draw_dir_tol))
  }
  # a hard wall: +Inf from 0.4 of the way out on every direction
  r <- one(function(p) if (abs(sum(u * p)) > 1.0) Inf else quad(p))
  cat("constructed hard wall at every read t: result",
      if (is.null(r)) "NULL (left unread, no fire)" else format(r$stat), "\n")
}

if ("B" %in% args) {
  cat("\n#### B. option (b): components of the firing directions\n")
  cat("draw_vc_comps:", ns$draw_vc_comps, "\n")
  ids <- c("slpw:12", "slpw:14", "slpw:15", "slpw:18", "slp10:10", "prsp:1",
           "prsp:10", "ordrs:9", "S0.1:14", "S0.1:16", "S0.1:24", "S0.2:8")
  tol <- ns$draw_dir_tol
  for (id in ids) {
    fit <- fit_of(id)
    dw <- grab_dw(fit)
    da <- dirs_all(fit, dw)
    fi <- which(da$stat > tol)
    pp <- da$R[(fi - 1L) %% da$k + 1L]
    cmp <- dw$comps[pp]
    s <- said(predict(fit, ndraws = 50, summary = FALSE))
    chk <- fit$cache[["draw_check"]]
    cat(sprintf("%-9s fired: %s | all-in-set %s any-in-set %s | cached warn_only %s | classes %s\n",
                id, paste(sprintf("%s%s(%s) %.1f",
                                  ifelse(fi > da$k, "-", "+"),
                                  dw$labels[pp], cmp, da$stat[fi]),
                          collapse = ", "),
                all(cmp %in% ns$draw_vc_comps),
                any(cmp %in% ns$draw_vc_comps), chk$warn_only,
                paste(unique(vapply(s$said, `[[`, "", "cls")),
                      collapse = ",")))
  }
  # test-gr-by.R's fit
  ex2 <- parse("tests/testthat/test-gr-by.R")
  for (e in ex2) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        as.character(e[[2]]) %in% c("grby_data")) eval(e, globalenv())
  }
  fit <- fq(frm(bf(y ~ x + (1 + x | gr(g, by = f))) + gaussian(),
                data = grby_data()))
  dw <- grab_dw(fit)
  da <- dirs_all(fit, dw)
  fi <- which(da$stat > tol)
  pp <- da$R[(fi - 1L) %% da$k + 1L]
  s <- said(predict(fit, ndraws = 50, summary = FALSE))
  chk <- fit$cache[["draw_check"]]
  cat(sprintf("gr-by fit fired: %s | warn_only %s | classes %s\n",
              paste(sprintf("%s%s(%s) %.1f", ifelse(fi > da$k, "-", "+"),
                            dw$labels[pp], dw$comps[pp], da$stat[fi]),
                    collapse = ", "), chk$warn_only,
              paste(unique(vapply(s$said, `[[`, "", "cls")),
                    collapse = ",")))
}

if ("C" %in% args || "D" %in% args) {
  fit <- s_fit("S0.1", 14)
  d <- fit$data
  nd_new <- data.frame(x = 2, g = factor("new"))
  nd_old <- data.frame(x = 2, g = factor("3", levels = levels(d$g)))
  W <- "frmtmb_draw_wide_tail"
  cnt <- function(lab, expr) {
    t0 <- proc.time()[[3]]
    s <- said(expr)
    msgs <- vapply(s$said, function(x) if (x$cls == W) x$msg else NA, "")
    msgs <- msgs[!is.na(msgs)]
    cat(sprintf("  %-58s wide_tail %d caveat %d (%.1f s) other: %s\n", lab,
                length(msgs),
                sum(grepl("cannot yet predict a grouping level", msgs,
                          fixed = TRUE)),
                proc.time()[[3]] - t0,
                paste(unique(setdiff(vapply(s$said, `[[`, "", "cls"), W)),
                      collapse = ",")))
    invisible(s$val)
  }
  if ("C" %in% args) {
    cat("\n#### C. warnings per call (S0.1:14)\n")
    cnt("predict() first call (check computed)", predict(fit, ndraws = 100))
    cat("  cached:", !is.null(fit$cache[["draw_check"]]), "warn_only",
        fit$cache[["draw_check"]]$warn_only, "\n")
    cnt("predict() second call", predict(fit, ndraws = 100))
    cnt("predict(summary = FALSE)", predict(fit, ndraws = 100,
                                            summary = FALSE))
    cnt("predict() in a loop of 3", for (i in 1:3) predict(fit, ndraws = 20))
    cnt("predict(propagate_error = FALSE)",
        predict(fit, ndraws = 20, propagate_error = FALSE))
    cnt("fitted()", fitted(fit))
    cnt("fitted(newdata = new level, allow_new_levels)",
        fitted(fit, newdata = nd_new, allow_new_levels = TRUE))
    cnt("conditional_effects()", conditional_effects(fit))
    cnt("conditional_effects(method = 'posterior_predict')",
        conditional_effects(fit, method = "posterior_predict"))
    cnt("pp_check()", pp_check(fit, ndraws = 10))
    cnt("loo_predictive_interval()",
        tryCatch(loo_predictive_interval(fit), error = function(e) {
          cat("    error:", conditionMessage(e), "\n")
        }))
    cnt("simulate()", simulate(fit, nsim = 2))
    # a fresh copy of the fit: the verdict is recomputed once
    fit2 <- fit
    fit2$cache <- new.env()
    cnt("predict() on a fresh cache", predict(fit2, ndraws = 100))
  }
  if ("D" %in% args) {
    cat("\n#### D. the new-level caveat (S0.1:14)\n")
    cnt("predict(new level, allow_new_levels = TRUE)",
        predict(fit, newdata = nd_new, allow_new_levels = TRUE, ndraws = 50))
    cnt("predict(old level, allow_new_levels = TRUE)",
        predict(fit, newdata = nd_old, allow_new_levels = TRUE, ndraws = 50))
    cnt("predict(old level)", predict(fit, newdata = nd_old, ndraws = 50))
    cnt("predict() in-sample", predict(fit, ndraws = 50))
    cnt("predict(new level, re_formula = NA, allow_new_levels = TRUE)",
        predict(fit, newdata = nd_new, re_formula = NA,
                allow_new_levels = TRUE, ndraws = 50))
    cnt("predict(new level, sample_new_levels = 'gaussian')",
        predict(fit, newdata = nd_new, allow_new_levels = TRUE,
                sample_new_levels = "gaussian", ndraws = 50))
    s <- said(predict(fit, newdata = nd_new, allow_new_levels = TRUE,
                      ndraws = 50))
    for (x in s$said) cat("  said", x$cls, ":", x$msg, "\n")
  }
}

if ("E" %in% args) {
  cat("\n#### E. slpw:12 reweighted arm (betad direction fires)\n")
  r <- design("slpw", 12L)
  fit <- r$fit
  set.seed(99)
  o <- r$onew(2e5)
  cover <- function(pn, o) {
    q90 <- stats::quantile(pn, c(0.05, 0.95), names = FALSE)
    q98 <- stats::quantile(pn, c(0.01, 0.99), names = FALSE)
    sprintf("cov90 %.3f w90 %8.2f cov98 %.3f w98 %8.2f",
            mean(o >= q90[1] & o <= q90[2]), q90[2] - q90[1],
            mean(o >= q98[1] & o <= q98[2]), q98[2] - q98[1])
  }
  set.seed(1)
  s <- said(predict(fit, ndraws = 20000, summary = FALSE))
  pin <- s$val
  k <- vapply(s$said, function(x) {
    m <- regmatches(x$msg, regexpr("Pareto k (of )?[0-9.]+", x$msg))
    if (length(m)) m else x$cls
  }, "")
  set.seed(1)
  pn <- fq(predict(fit, newdata = r$nd, allow_new_levels = TRUE,
                   ndraws = 20000, summary = FALSE))[, 1]
  chk <- fit$cache[["draw_check"]]
  cat(sprintf("  lane  mae %.4f | %s | %s | fire %s warn_only %s dir %.3f (%s)\n",
              mean(abs(colMeans(pin) - r$mu)), paste(k, collapse = ","),
              cover(pn, o), chk$fire, chk$warn_only, chk$dir, chk$dir_par))
  cat("  review first build: mae 7.4062 | k=0.8 | cov90 0.839 w90 114.78 cov98 0.957 w98 165.38\n")
}
cat("DONE\n")
