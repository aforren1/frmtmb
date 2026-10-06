# Defect 8, measured over seeds of brms_monotonic's data code.
# Per seed and fit (mo(income) * age, and mo(income) alone): whether
# sdreport's SEs are all non-finite, the simplex values, how many
# simplex components sit at 0 (below 1e-8), and what three inverses of
# the same outer Hessian give:
#   plain    solve(H), what sdreport does
#   scaled   solve(D^-1 H D^-1) with D = sqrt|diag H|, mapped back
#   held     the scaled inverse over the parameters off the null space
#            (eigenvalue of the scaled H below 1e-12 of the largest),
#            the null-space parameters held
# Output: one TSV row per seed and fit.
#
#   Rscript dev/nanse-mo-sweep.R [lib] [seeds] [out]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r5"
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:200
out <- if (length(args) > 2) args[3] else "dev/nanse-log/mo-sweep.tsv"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
mk <- function(s) {
  set.seed(s)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}
inv_scaled <- function(H) {
  D <- sqrt(abs(diag(H)))
  V <- try(solve(H / outer(D, D)), silent = TRUE)
  if (inherits(V, "try-error")) return(H * NaN)
  V / outer(D, D)
}
inv_held <- function(H, tol = 1e-12) {
  D <- sqrt(abs(diag(H)))
  held <- !(D > 0)
  pos <- which(!held)
  S <- H / outer(D, D)
  e <- eigen(S[pos, pos, drop = FALSE], symmetric = TRUE)
  null <- abs(e$values) < tol * max(abs(e$values))
  if (any(null)) {
    held[pos] <- apply(abs(e$vectors[, null, drop = FALSE]) > 0.1, 1, any)
  }
  V <- H * NaN
  k <- !held
  Vk <- try(solve(S[k, k, drop = FALSE]), silent = TRUE)
  if (!inherits(Vk, "try-error")) V[k, k] <- Vk / outer(D[k], D[k])
  list(V = V, held = held)
}
fmt <- function(x) paste(signif(x, 4), collapse = ",")
rows <- list()
for (s in seeds) {
  d <- mk(s)
  for (form in c("int", "main")) {
    fml <- if (form == "int") ls ~ mo(income) * age else ls ~ mo(income)
    w <- character()
    f <- withCallingHandlers(frm(fml, data = d), warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
    sdr <- frmtmb:::sdr_of(f)
    nm <- frmtmb:::outer_par_names(f)
    se0 <- sqrt(diag(sdr$cov.fixed))
    p <- f$opt$par
    H <- optimHess(p, f$obj$fn, f$obj$gr)
    Vs <- inv_scaled(H)
    hv <- inv_held(H)
    sese <- suppressWarnings(sqrt(diag(Vs)))
    seh <- suppressWarnings(sqrt(diag(hv$V)))
    zt <- f$estimates[grepl("^zeta", names(f$estimates))]
    simp <- lapply(zt, function(z) {
      x <- exp(c(0, z))
      x / sum(x)
    })
    smin <- vapply(simp, min, 0)
    dh <- sqrt(abs(diag(H)))
    ok <- dh > 0
    ev <- eigen((H / outer(dh, dh))[ok, ok], TRUE, only.values = TRUE)$values
    b <- nm %in% c("(Intercept)", "age", "moincome", "moincome:age",
                   "sigma_(Intercept)")
    rows[[length(rows) + 1L]] <- data.frame(
      seed = s, form = form, code = f$opt$convergence,
      n_warn = length(w), pdHess = isTRUE(sdr$pdHess),
      se_plain_finite = sum(is.finite(se0)), n_par = length(se0),
      se_scaled_finite = sum(is.finite(sese)),
      se_held_finite = sum(is.finite(seh)),
      n_held = sum(hv$held), held = paste(nm[hv$held], collapse = ","),
      n_zero_diag = sum(!ok),
      n_sat = sum(unlist(simp) < 1e-8),
      simplex_min = fmt(smin),
      rcond_H = signif(rcond(H), 3),
      mineig_scaled = signif(min(ev) / max(ev), 3),
      maxgrad = signif(max(abs(f$obj$gr(p))), 3),
      se_b_plain = fmt(se0[b]), se_b_scaled = fmt(sese[b]),
      se_b_held = fmt(seh[b]),
      est_b = fmt(p[b]),
      zeta = fmt(unlist(zt)),
      stringsAsFactors = FALSE)
  }
  if (s %% 20 == 0) cat("seed", s, "\n")
}
X <- do.call(rbind, rows)
utils::write.table(X, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", nrow(X), "rows to", out, "\n")
