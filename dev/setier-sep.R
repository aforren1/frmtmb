# Lane setier: separation named at fit time (check_convergence()'s
# separation_check()), against an independent reference, on separated
# and on ordinary designs.
#   Rscript dev/setier-sep.R <lib or "base"> <out .tsv> [seeds]
# The reference is Konis's linear program for separation (brglm2's
# detect_separation(), which is not installed here), solved as a
# quadratic program (quadprog, the projection of sum(A) onto the cone): the data
# are separated iff some b with s_i x_i'b >= 0 for every
# observation (= 0 where both outcomes occur) has x'b != 0.
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
if (identical(args[1], "src")) {
  suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-setier",
                                     quiet = TRUE))
} else suppressPackageStartupMessages(library(frmtmb))
out <- args[2]
S <- if (length(args) >= 3) as.integer(args[3]) else 50L
cat("lib:", find.package("frmtmb"), " BLAS probe (s):",
    system.time({m <- matrix(1, 600, 600); m %*% m})[["elapsed"]], "\n")
conds <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  }, message = function(c) invokeRestart("muffleMessage"))
  list(value = val, warnings = w)
}
lp_separated <- function(X, y, size) {
  s <- ifelse(y >= size, 1, ifelse(y <= 0, -1, 0))
  X <- sweep(X, 2L, pmax(sqrt(colSums(X^2)), 1e-12), `/`)
  A <- X * s
  p <- ncol(X)
  eq <- s == 0
  # equalities as an orthonormal basis of their row space, so that a
  # redundant one does not make solve.QP() call them inconsistent
  E <- matrix(0, 0, p)
  if (any(eq)) {
    q <- qr(t(X[eq, , drop = FALSE]))
    if (q$rank) E <- t(qr.Q(q)[, seq_len(q$rank), drop = FALSE])
  }
  Amat <- t(rbind(E, A[!eq, , drop = FALSE]))
  bvec <- rep(0, nrow(E) + sum(!eq))
  dvec <- colSums(A[!eq, , drop = FALSE])
  # min |b|^2 / 2 - d'b over the cone: b* is the projection of d onto
  # it, nonzero (with A b* != 0) exactly when some b in the cone has
  # d'b = sum(A b) > 0, that is when the data are separated
  sol <- tryCatch(quadprog::solve.QP(diag(1, p), dvec, Amat, bvec,
                                     meq = nrow(E)),
                  error = function(e) NULL)
  if (is.null(sol)) return(NA)
  m <- drop(X %*% sol$solution)
  max(abs(m)) > 1e-6
}
gen <- function(design, seed) {
  set.seed(seed)
  size <- 1
  if (design == "plain") {
    d <- data.frame(x = rnorm(100)); d$y <- rbinom(100, 1, plogis(d$x))
    ff <- y ~ x
  } else if (design == "strong") {
    d <- data.frame(x = rnorm(400))
    d$y <- rbinom(400, 1, plogis(-1 + 4 * d$x)); ff <- y ~ x
  } else if (design == "rare") {
    d <- data.frame(x = rnorm(1000))
    d$y <- rbinom(1000, 1, plogis(-5 + d$x)); ff <- y ~ x
  } else if (design == "cells") {
    d <- data.frame(f = factor(sample(letters[1:5], 60, TRUE)),
                    x = rnorm(60))
    d$y <- rbinom(60, 1, plogis(-1 + 0.5 * d$x)); ff <- y ~ f + x
  } else if (design == "trials") {
    d <- data.frame(x = rnorm(200)); size <- 5
    d$k <- rbinom(200, 5, plogis(1 + 3 * d$x)); d$nk <- 5 - d$k
    ff <- cbind(k, nk) ~ x
  } else if (design == "glmm") {
    d <- data.frame(g = factor(rep(1:20, each = 10)), x = rnorm(200))
    d$y <- rbinom(200, 1, plogis(0.5 * d$x + rnorm(20)[d$g]))
    ff <- y ~ x + (1 | g)
  } else if (design == "complete") {
    d <- data.frame(x = rnorm(50), z = rnorm(50))
    d$y <- as.integer(d$x > 0); ff <- y ~ x + z
  } else if (design == "quasi") {
    d <- data.frame(f = factor(rep(c("a", "b", "c"), each = 40)),
                    x = rnorm(120))
    d$y <- rbinom(120, 1, plogis(0.3 + 0.5 * d$x)); d$y[d$f == "c"] <- 0L
    ff <- y ~ f + x
  } else if (design == "quasi_x") {
    d <- data.frame(x = c(rnorm(40, -2), rep(0, 20), rnorm(40, 2)))
    d$y <- c(rep(0L, 40), rbinom(20, 1, 0.5), rep(1L, 40)); ff <- y ~ x
  } else if (design == "glmm_quasi") {
    d <- data.frame(g = factor(rep(1:20, each = 10)),
                    f = factor(rep(c("a", "b"), 100)), x = rnorm(200))
    d$y <- rbinom(200, 1, plogis(0.5 * d$x + rnorm(20)[d$g]))
    d$y[d$f == "b"] <- 0L; ff <- y ~ f + x + (1 | g)
  }
  list(d = d, ff = ff, size = size)
}
fam_of <- function(g) if (g$size > 1) binomial() else bernoulli()
rows <- list()
designs <- c("plain", "strong", "rare", "cells", "trials", "glmm",
             "complete", "quasi", "quasi_x", "glmm_quasi")
for (design in designs) {
  for (seed in seq_len(S)) {
    g <- gen(design, seed)
    ff_fixed <- reformulas::nobars(g$ff)
    mf <- model.frame(ff_fixed, g$d)
    X <- model.matrix(ff_fixed, mf)
    yv <- if (g$size > 1) g$d$k else g$d$y
    ref <- lp_separated(X, yv, g$size)
    gl <- conds(glm(ff_fixed, family = binomial(), data = g$d))
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(conds(frm(g$ff, family = fam_of(g), data = g$d)),
                  error = function(e) NULL)
    t1 <- proc.time()[["elapsed"]]
    if (is.null(r)) {
      rows[[length(rows) + 1L]] <- data.frame(design = design, seed = seed,
        lp_separated = ref, glm_01 = NA, code = NA, evals = NA,
        sep_warn = NA, other_warn = NA, first = "ERROR", fit_s = t1 - t0)
      next
    }
    sepw <- grepl("separate the outcomes", r$warnings)
    rows[[length(rows) + 1L]] <- data.frame(
      design = design, seed = seed, lp_separated = ref,
      glm_01 = any(grepl("numerically 0 or 1", gl$warnings)),
      code = r$value$opt$convergence, evals = r$value$opt$evals %||% NA,
      sep_warn = sum(sepw), other_warn = sum(!sepw),
      first = substr(c(r$warnings, "")[1L], 1, 80), fit_s = t1 - t0)
  }
}
res <- do.call(rbind, rows)
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
for (d in designs) {
  x <- res[res$design == d, ]
  cat(sprintf(paste0("%-10s n=%2d LP separated %2d | named %2d (of LP ",
                     "separated %2d, of not %2d) | other warning %2d | ",
                     "code!=0 %2d | glm 0/1 %2d | median evals %5.0f ",
                     "fit s %.2f\n"),
              d, nrow(x), sum(x$lp_separated, na.rm = TRUE),
              sum(x$sep_warn > 0, na.rm = TRUE),
              sum(x$sep_warn > 0 & x$lp_separated, na.rm = TRUE),
              sum(x$sep_warn > 0 & !x$lp_separated, na.rm = TRUE),
              sum(x$other_warn > 0, na.rm = TRUE),
              sum(x$code != 0, na.rm = TRUE), sum(x$glm_01, na.rm = TRUE),
              stats::median(x$evals, na.rm = TRUE),
              stats::median(x$fit_s)))
}
