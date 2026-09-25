# Summarise the recovery runs from the RDS files alone. Counts come
# from the files that exist, never from the launcher, and the seed grid
# is checked for gaps.
#
# Output: dev/phase3b-log/summary-<arm>.txt, pasted verbatim into
# dev/phase3b-findings.md.
.libPaths(c("C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
arms <- commandArgs(trailingOnly = TRUE)
dir <- "dev/phase3b-log/recov"

wilson <- function(k, n) {
  z <- 1.959964
  p <- k / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n),
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n))
}

row_of <- function(name, truth, est, lwr, upr) {
  ok <- is.finite(est) & is.finite(lwr) & is.finite(upr)
  k <- sum(lwr[ok] < truth[ok] & truth[ok] < upr[ok])
  n <- sum(ok)
  w <- wilson(k, n)
  sprintf("%-22s %9.4f %9.4f %8.4f %4d/%-4d %6.1f%%  [%4.1f, %4.1f]",
          name, mean(truth), mean(est[ok]), sd(est[ok]) / sqrt(n), k, n,
          100 * k / n, 100 * w[1], 100 * w[2])
}

header <- sprintf("%-22s %9s %9s %8s %9s %7s  %s", "quantity", "truth",
                  "mean", "mcse", "covered", "rate", "Wilson 95")

theta_rows <- function(ci) grep("^theta_", rownames(ci), value = TRUE)

summ_eam <- function(arm) {
  f <- list.files(dir, sprintf("^%s-[0-9]+[.]rds$", arm), full.names = TRUE)
  x <- lapply(f, readRDS)
  # The rerun pass (dev/phase3b-eam-recovery3.R) refitted every seed
  # that stopped with an error, sequentially and on the same build,
  # plus two seeds that had succeeded as a control. A rerun replaces
  # an ERRORED replicate only; the control says whether a rerun
  # reproduces a success.
  rr <- list.files("dev/phase3b-log/recov-rerun",
                   sprintf("^%s-[0-9]+[.]rds$", arm), full.names = TRUE)
  rerun_note <- character(0)
  if (length(rr)) {
    y <- lapply(rr, readRDS)
    ys <- vapply(y, `[[`, 0, "seed")
    xs <- vapply(x, `[[`, 0, "seed")
    was_err <- vapply(x, function(r) !is.null(r$error), TRUE)
    ctrl <- ys %in% xs[!was_err]
    for (k in which(ctrl)) {
      a <- x[[match(ys[k], xs)]]
      b <- y[[k]]
      rerun_note <- c(rerun_note, sprintf(
        "control seed %d: logLik first pass %.10f, rerun %.10f, identical: %s",
        ys[k], a$loglik, b$loglik, identical(a$fixef, b$fixef)))
    }
    fixed <- 0L
    for (k in which(!ctrl)) {
      j <- match(ys[k], xs)
      if (!is.na(j) && was_err[j] && is.null(y[[k]]$error)) {
        x[[j]] <- y[[k]]
        fixed <- fixed + 1L
      }
    }
    rerun_note <- c(sprintf(paste0("first-pass errors: %d; refitted without ",
                                   "error on rerun: %d"),
                            sum(was_err), fixed), rerun_note)
  }
  seeds <- vapply(x, `[[`, 0, "seed")
  err <- vapply(x, function(r) !is.null(r$error), TRUE)
  out <- c(sprintf("arm %s: %d files, seeds %d to %d, gaps: %s", arm,
                   length(x), min(seeds), max(seeds),
                   paste(setdiff(seq(min(seeds), max(seeds)), seeds),
                         collapse = " ")),
           rerun_note,
           sprintf("errors after the rerun: %d", sum(err)))
  x <- x[!err]
  conv <- vapply(x, function(r) r$conv, 0)
  # diagnose() prints three standard lines and then either its verdict
  # or its findings; with a gradient between its two thresholds it
  # prints neither. A replicate where diagnose() itself failed carries
  # the error text instead.
  dfail <- vapply(x, function(r) {
    !any(grepl("^Optimizer convergence code", r$diagnose))
  }, TRUE)
  flag <- vapply(x, function(r) {
    length(r$diagnose) > 3L &&
      !any(grepl("No convergence problems detected", r$diagnose))
  }, TRUE) & !dfail
  out <- c(out, sprintf(paste0("convergence code 0: %d of %d; diagnose() ",
                               "reported a finding: %d; diagnose() failed: %d"),
                        sum(conv == 0), length(x), sum(flag), sum(dfail)))
  # what diagnose() said, by the first words of each finding after its
  # three standard lines
  why <- unlist(lapply(x, function(r) {
    d <- r$diagnose[-(1:3)]
    d <- d[nzchar(trimws(d)) & !grepl("^ ", d)]
    unique(substr(d, 1, 48))
  }))
  if (length(why)) {
    tb <- sort(table(why), decreasing = TRUE)
    out <- c(out, sprintf("  %3d  %s", as.integer(tb), names(tb)))
  }
  mg <- vapply(x, function(r) {
    l <- grep("^Max [|]gradient[|]", r$diagnose, value = TRUE)
    if (length(l)) as.numeric(sub("^Max [|]gradient[|]: ([^ ]+) .*", "\\1", l[1]))
    else NA_real_
  }, 0)
  out <- c(out, sprintf("max |gradient| at the optimum: median %.2g, max %.2g",
                        median(mg, na.rm = TRUE), max(mg, na.rm = TRUE)))
  ci_ok <- vapply(x, function(r) !is.null(r$confint), TRUE)
  out <- c(out, sprintf("fits with a confint(): %d of %d", sum(ci_ok), length(x)))
  get <- function(row, col) vapply(x, function(r) {
    ci <- r$confint
    if (row %in% rownames(ci)) ci[row, col] else NA_real_
  }, 0)
  # the non-decision time is a scaled logit onto (0, ub), ub the bound
  # the fit used, so its truth on the link scale is per fit
  ub <- vapply(x, function(r) r$ndt_ub, 0)
  ndt_truth <- qlogis(0.25 / ub)
  ndt_truth[0.25 >= ub] <- NA
  out <- c(out, "", header,
           row_of("mu intercept", rep(0.4, length(x)),
                  get("(Intercept)", "est"), get("(Intercept)", "lwr"),
                  get("(Intercept)", "upr")),
           row_of("mu condition", rep(0.9, length(x)), get("cond", "est"),
                  get("cond", "lwr"), get("cond", "upr")),
           row_of("log bs", rep(log(1.4), length(x)),
                  get("bs_(Intercept)", "est"), get("bs_(Intercept)", "lwr"),
                  get("bs_(Intercept)", "upr")))
  ok_ndt <- is.finite(ndt_truth)
  out <- c(out, sprintf("ndt: the true 0.25 is below the fitted bound on %d of %d",
                        sum(ok_ndt), length(x)))
  if (any(ok_ndt)) {
    out <- c(out, row_of("ndt (link scale)", ndt_truth[ok_ndt],
                         get("ndt_(Intercept)", "est")[ok_ndt],
                         get("ndt_(Intercept)", "lwr")[ok_ndt],
                         get("ndt_(Intercept)", "upr")[ok_ndt]))
  }
  out <- c(out, sprintf("ndt natural scale: mean %.4f (truth 0.25)",
                        mean(plogis(get("ndt_(Intercept)", "est")) * ub, na.rm = TRUE)))
  # the two standard deviations: theta rows in the order VarCorr lists
  th <- theta_rows(x[[1]]$confint)
  sds <- c(log(0.35), log(0.20))
  for (j in seq_along(th)) {
    out <- c(out, row_of(sprintf("log sd %d (%s)", j,
                                 c("mu | s", "log bs | s")[j]),
                         rep(sds[j], length(x)), get(th[j], "est"),
                         get(th[j], "lwr"), get(th[j], "upr")))
  }
  if (arm %in% c("cont", "contdef", "contfix", "collapse")) {
    le <- get("lambda_(Intercept)", "est")
    if (arm != "collapse") {
      out <- c(out, row_of("lambda (logit)", rep(qlogis(0.05), length(x)),
                           le, get("lambda_(Intercept)", "lwr"),
                           get("lambda_(Intercept)", "upr")))
      out <- c(out, sprintf("lambda natural: mean %.4f, sd %.4f (truth 0.05)",
                            mean(plogis(le), na.rm = TRUE),
                            sd(plogis(le), na.rm = TRUE)))
      nc <- vapply(x, function(r) r$n_cont, 0)
      out <- c(out, sprintf("realized contaminant share: mean %.4f",
                            mean(nc / 12000)))
    } else {
      dll <- vapply(x, function(r) r$loglik - r$loglik_plain, 0)
      edge <- le < -15
      lflag <- vapply(x, function(r) {
        any(grepl("end of its link: lambda", r$diagnose))
      }, TRUE)
      out <- c(out,
               sprintf("lambda at the link edge (logit below -15): %d of %d",
                       sum(edge), length(x)),
               sprintf("diagnose() names lambda at the end of its link: %d; of those at the edge: %d",
                       sum(lflag), sum(lflag & edge)),
               sprintf("interior lambda: median %.2g, min %.2g, max %.2g",
                       median(plogis(le[!edge])), min(plogis(le[!edge])),
                       max(plogis(le[!edge]))),
               sprintf("logLik(contaminant) - logLik(plain): min %.3g, median %.3g, max %.3g",
                       min(dll, na.rm = TRUE), median(dll, na.rm = TRUE),
                       max(dll, na.rm = TRUE)),
               sprintf("  at the edge: max |difference| %.3g",
                       max(abs(dll[edge]), na.rm = TRUE)))
    }
  }
  if (arm == "cens") {
    nc <- vapply(x, function(r) r$n_cens, 0)
    out <- c(out, sprintf("censored share: mean %.4f, range %.4f to %.4f",
                          mean(nc / 12000), min(nc / 12000),
                          max(nc / 12000)))
  }
  el <- vapply(x, function(r) r$elapsed, 0)
  out <- c(out, sprintf("seconds per replicate: median %.0f", median(el)))
  out
}

summ_session <- function() {
  f <- list.files(dir, "^session-[0-9]+[.]rds$", full.names = TRUE)
  x <- lapply(f, readRDS)
  seeds <- vapply(x, `[[`, 0, "seed")
  out <- sprintf("session: %d files, seeds %d to %d, gaps: %s", length(x),
                 min(seeds), max(seeds),
                 paste(setdiff(seq(min(seeds), max(seeds)), seeds),
                       collapse = " "))
  for (which in c("with_session", "without")) {
    y <- lapply(x, `[[`, which)
    err <- vapply(y, function(r) !is.null(r$error), TRUE)
    y <- y[!err]
    get <- function(row, col) vapply(y, function(r) r$confint[row, col], 0)
    conv <- vapply(y, function(r) r$conv, 0)
    out <- c(out, "", sprintf("== %s: %d fits, %d errors, code 0 on %d", which,
                              length(y), sum(err), sum(conv == 0)),
             header,
             row_of("logit alpha", rep(qlogis(0.35), length(y)),
                    get("alpha_(Intercept)", "est"),
                    get("alpha_(Intercept)", "lwr"),
                    get("alpha_(Intercept)", "upr")),
             row_of("log tau", rep(log(3), length(y)),
                    get("tau_(Intercept)", "est"),
                    get("tau_(Intercept)", "lwr"),
                    get("tau_(Intercept)", "upr")),
             row_of("log sd(alpha | id)", rep(log(0.5), length(y)),
                    get("theta_1", "est"), get("theta_1", "lwr"),
                    get("theta_1", "upr")))
  }
  out
}

for (arm in arms) {
  res <- if (arm == "session") summ_session() else summ_eam(arm)
  writeLines(res, sprintf("dev/phase3b-log/summary-%s.txt", arm))
  cat(res, sep = "\n")
  cat("\n")
}
