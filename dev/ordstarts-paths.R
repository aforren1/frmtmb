# Copy of dev/ordstarts-refits.R (which was running when this was made)
# with B with starts = 0 added to the influence() variants.
# Lane ordstarts (backlog 6136), items 2 and 3 on the seeds whose answer
# the starts moved (dev/ordstarts-sum.R writes the list). Per seed, the
# data and model of dev/ordstarts-run.R; A is the default fit, B the
# fit with frmtmb_control(starts = 10).
#
# real:     is B's point a maximum? Four refits from B's estimates, each
#           outer coordinate moved by N(0, 0.1) (seed 7000 + k), with no
#           starts; a code-0 end above B by more than 1e-3 says it is
#           not. The same about A, for the comparison.
# refit:    refit(A, y) and refit(B, y) on the same response, and B with
#           starts = 0: whether starts run (the starts record) and where
#           each ends.
# boot:     frm_bootstrap(nsim = 5, seed = 1) of A as it is and of A with
#           starts = 10 in its control, which simulate the same
#           responses; and of B (starts inherited) and of B with
#           starts = 0. Per replicate: logLik, code, the starts record.
# influence (n = 300 only): influence(groups = "g") of A, of A with
#           starts = 10 and of B with starts = 0, each deletion refit
#           read off
#           interop_coef_vector(), which influence() calls on it.
# allfit:   frm_allfit() of A and of B: each optimizer's logLik and
#           whether any ran starts.
#   Rscript dev/ordstarts-refits.R <list.tsv> <rows> <out.tsv> [checks]
#     rows: R expression for the list's rows; checks: comma list of
#     real,refit,boot,influence,allfit (default all)
args <- commandArgs(trailingOnly = TRUE)
lst <- utils::read.delim(args[1], stringsAsFactors = FALSE)
rows <- eval(parse(text = args[2]))
out <- args[3]
checks <- if (length(args) > 3) strsplit(args[4], ",")[[1]] else
  c("real", "refit", "boot", "influence", "allfit")
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
# the data of dev/ordstarts-run.R
mkg <- function(seed, nobs, G, re) {
  d <- mk(seed, n = nobs, mix = TRUE)
  set.seed(seed)
  x <- rnorm(nobs)
  z <- rnorm(nobs)
  cls <- rbinom(nobs, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(nobs)
  stopifnot(identical(x, d$x), identical(z, d$z))
  g <- rep_len(seq_len(G), nobs)
  u0 <- rnorm(G, 0, 0.5)
  lat <- lat + u0[g]
  if (re %in% c("slope", "slopeu")) {
    u1 <- rnorm(G, 0, 0.3)
    lat <- lat + u1[g] * x
  }
  d$y <- 1L + (lat > -1 + 0.1 * z) + (lat > 0) + (lat > 1 - 0.1 * z)
  d$g <- factor(g)
  d
}
forms <- list(int = bf(y ~ x + z + (1 | g)),
              slope = bf(y ~ x + z + (1 + x | g)),
              slopeu = bf(y ~ x + z + (1 + x || g)))
fams <- list(cum_cum = function() mixture(cumulative(), cumulative()),
             cum_sratio = function() mixture(cumulative(), sratio()))
quiet <- function(e) suppressWarnings(suppressMessages(e))
tm <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  v <- tryCatch(quiet(expr), error = function(e) e)
  list(v = v, el = proc.time()[["elapsed"]] - t0)
}
rec <- function(f) {
  if (inherits(f, "error") || is.null(f)) {
    return(c(ll = NA, code = NA, ms_starts = NA, ms_reached = NA,
             ms_gain = NA, ms_tier = NA))
  }
  ms <- f$opt$mixture_starts
  c(ll = as.numeric(logLik(f)), code = f$opt$convergence,
    ms_starts = if (is.null(ms)) 0 else ms[["starts"]],
    ms_reached = if (is.null(ms)) NA else ms[["reached"]],
    ms_gain = if (is.null(ms)) NA else ms[["gain"]],
    ms_tier = if (is.null(ms)) NA else ms[["tier"]])
}
first <- !file.exists(out)
emit <- function(key, check, variant, k, r, el, note = "") {
  row <- data.frame(key, check = check, variant = variant, k = k,
                    t(r), elapsed = el, note = note)
  utils::write.table(row, out, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = first, append = !first)
  first <<- FALSE
}
with_starts <- function(f, n) {
  f$control$starts <- as.integer(n)
  f
}
for (i in rows) {
  L <- lst[i, ]
  key <- data.frame(family = L$family, re = L$re, n = L$n, G = L$G,
                    seed = L$seed)
  d <- mkg(L$seed, L$n, L$G, L$re)
  fit <- function(ctl) {
    frm(forms[[L$re]], family = fams[[L$family]](), data = d, control = ctl)
  }
  a <- tm(fit(frmtmb_control()))
  b <- tm(fit(frmtmb_control(starts = 10)))
  emit(key, "fit", "A", 0, rec(a$v), a$el)
  emit(key, "fit", "B", 0, rec(b$v), b$el)
  A <- a$v
  B <- b$v
  if (inherits(A, "error") || inherits(B, "error")) next
  if ("real" %in% checks) {
    # which parameters carry the largest standard errors at each point
    # (the forced fit-end covariance), and where they sit
    for (who in c("A", "B")) {
      F <- if (who == "A") A else B
      an <- tryCatch(quiet(ns$se_fit_analysis(F, force = TRUE)),
                     error = function(e) NULL)
      note <- "no covariance"
      if (is.matrix(an$V)) {
        s <- sqrt(pmax(diag(an$V), 0))
        nm <- tryCatch(ns$outer_par_map(F)$names, error = function(e) NULL)
        if (length(nm) != length(s)) nm <- names(F$opt$par)
        top <- utils::head(order(-s), 3L)
        note <- paste0(nm[top], " se ", signif(s[top], 3), " at ",
                       signif(F$opt$par[top], 3), collapse = "; ")
      }
      emit(key, "se", who, 0, rec(F), NA, note)
    }
    for (who in c("A", "B")) {
      F0 <- with_starts(if (who == "A") A else B, 0)
      est <- F0$estimates
      for (k in 1:4) {
        set.seed(7000 + k)
        st <- est
        for (nm in setdiff(names(st), "b")) {
          st[[nm]] <- st[[nm]] + stats::rnorm(length(st[[nm]]), 0, 0.1)
        }
        r <- tm(refit(F0, d$y, start = st))
        emit(key, "real", who, k, rec(r$v), r$el,
             if (inherits(r$v, "error")) substr(conditionMessage(r$v), 1,
                                                120) else "")
      }
    }
  }
  if ("refit" %in% checks) {
    r <- tm(refit(A, d$y))
    emit(key, "refit", "A", 0, rec(r$v), r$el)
    r <- tm(refit(B, d$y))
    emit(key, "refit", "B", 0, rec(r$v), r$el)
    r <- tm(refit(with_starts(B, 0), d$y))
    emit(key, "refit", "B_starts0", 0, rec(r$v), r$el)
  }
  if ("boot" %in% checks) {
    FUN <- function(f) rec(f)
    for (v in c("A", "A_starts10", "B", "B_starts0")) {
      F <- switch(v, A = A, A_starts10 = with_starts(A, 10), B = B,
                  B_starts0 = with_starts(B, 0))
      r <- tm(frm_bootstrap(F, FUN = FUN, nsim = 5, seed = 1))
      if (inherits(r$v, "error")) {
        emit(key, "boot", v, 0, rec(NULL), r$el,
             substr(conditionMessage(r$v), 1, 120))
        next
      }
      for (k in 1:5) emit(key, "boot", v, k, r$v$t[k, ], r$el / 5)
    }
  }
  if ("influence" %in% checks && L$n <= 300) {
    for (v in c("A", "A_starts10", "B_starts0")) {
      F <- switch(v, A = A, A_starts10 = with_starts(A, 10),
                  B_starts0 = with_starts(B, 0))
      cap <- list()
      .GlobalEnv$.ordstarts_cap <- function(m) {
        cap[[length(cap) + 1L]] <<- rec(m)
      }
      suppressMessages(trace("interop_coef_vector", where = ns,
                             tracer = quote(.ordstarts_cap(model)),
                             print = FALSE))
      r <- tm(stats::influence(F, groups = "g"))
      suppressMessages(untrace("interop_coef_vector", where = ns))
      # the first call reads the full fit; one per deletion after it
      el1 <- r$el / max(length(cap) - 1L, 1L)
      for (k in seq_along(cap)) {
        emit(key, "influence", v, k - 1L, cap[[k]], el1,
             if (inherits(r$v, "error")) substr(conditionMessage(r$v), 1,
                                                120) else "")
      }
    }
  }
  if ("allfit" %in% checks) {
    for (v in c("A", "B")) {
      F <- if (v == "A") A else B
      r <- tm(frm_allfit(F))
      if (inherits(r$v, "error")) {
        emit(key, "allfit", v, 0, rec(NULL), r$el,
             substr(conditionMessage(r$v), 1, 120))
        next
      }
      for (k in seq_along(r$v$fits)) {
        emit(key, "allfit", v, k, rec(r$v$fits[[k]]), r$v$times[k],
             names(r$v$fits)[k])
      }
    }
  }
  cat("done", L$family, L$re, L$n, L$seed, "\n")
}
cat("DONE refits\n")
