# Lane optima, item 1, prototype measurement: on top of a fit, refit
# from the opposite sign of each mo() coefficient (simplex coordinates
# reset to the barycenter), and for the softmax also the vertex escape
# of dev/nanse-mo-escape.R, and count the seeds reaching the exact
# maximum of dev/optima-mo-study.R.
#   Rscript dev/optima-mo-flip.R soft|sphere seeds out.tsv
args <- commandArgs(trailingOnly = TRUE)
param <- args[1]
seeds <- eval(parse(text = args[2]))
out <- args[3]
Sys.setenv(FRMTMB_MO_PARAM = if (param == "sphere") "sphere" else "")
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("param", param, "frmtmb from", find.package("frmtmb"), "\n")
src <- readLines("dev/optima-mo-study.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^simplex_of", src) - 1L)]))
ms <- utils::getFromNamespace("mo_simplex", "frmtmb")

nl <- function(obj, p) {
  o <- stats::nlminb(p, obj$fn, obj$gr,
                     control = list(eval.max = 1000, iter.max = 1000))
  o$par <- stats::setNames(o$par, names(p))
  o
}
vertex_escape <- function(obj, par, sat = 1e-6, eps = 0.01, rounds = 10) {
  nm <- names(par)
  blocks <- split(seq_along(par), nm)
  blocks <- blocks[grepl("^zeta", names(blocks))]
  best <- par
  fbest <- obj$fn(par)
  n <- 0L
  for (r in seq_len(rounds)) {
    cand <- NULL
    fc <- fbest
    for (b in blocks) {
      w <- ms(best[b])
      if (min(w) > sat) next
      for (j in seq_along(w)) {
        e <- numeric(length(w))
        e[j] <- 1
        w2 <- (1 - eps) * w + eps * e
        p <- best
        p[b] <- log(w2[-1] / w2[1])
        f2 <- obj$fn(p)
        if (is.finite(f2) && f2 < fc - 1e-8) {
          fc <- f2
          cand <- p
        }
      }
    }
    if (is.null(cand)) break
    o <- nl(obj, cand)
    n <- n + o$evaluations[[1]]
    if (o$objective < fbest - 1e-10) {
      best <- o$par
      fbest <- o$objective
    } else break
  }
  list(par = best, objective = fbest, evals = n)
}
flip <- function(f, par, fbest) {
  obj <- f$obj
  nm <- names(par)
  bnames <- names(f$frame$par_template$beta)
  ib <- which(nm == "beta")
  mo <- f$frame$linpreds[[1]]$mo
  n <- 0L
  for (mi in mo) {
    p <- par
    j <- ib[f$frame$linpreds[[1]]$idx[mi$col]]
    p[j] <- -p[j]
    p[nm == mi$zeta] <- 0
    o <- nl(obj, p)
    n <- n + o$evaluations[[1]]
    if (is.finite(o$objective) && o$objective < fbest - 1e-10) {
      par <- o$par
      fbest <- o$objective
    }
  }
  list(par = par, objective = fbest, evals = n)
}
# the opposite sign cone of each saturated term, with the sign held by a
# bound so the run cannot slide back through zero, then a free polish
cone <- function(f, par, fbest, sat = 1e-8) {
  obj <- f$obj
  nm <- names(par)
  ib <- which(nm == "beta")
  mo <- f$frame$linpreds[[1]]$mo
  n <- 0L
  for (mi in mo) {
    if (min(ms(par[nm == mi$zeta])) > sat) next
    p <- par
    j <- ib[f$frame$linpreds[[1]]$idx[mi$col]]
    s <- sign(p[j])
    p[j] <- -p[j]
    p[nm == mi$zeta] <- 0
    if (identical(Sys.getenv("OPTIMA_SCREEN"), "1")) {
      # which vertex direction leads into the opposite cone fastest, read
      # off the gradient in b at b = 0 (others held)
      D <- mi$D
      Q <- utils::getFromNamespace("mo_sphere_basis", "frmtmb")(D)
      P <- rep(-1 / sqrt(D), D)
      tox <- function(w) {
        u <- sqrt(w)
        as.numeric(crossprod(Q, u)) / (1 - sum(u * P))
      }
      sc <- vapply(seq_len(D), function(k) {
        q <- par
        q[j] <- 0
        q[nm == mi$zeta] <- tox(replace(numeric(D), k, 1))
        s * obj$gr(q)[j]
      }, 0)
      k <- which.max(sc)
      if (sc[k] <= 0 && identical(Sys.getenv("OPTIMA_SKIP"), "1")) next
      w0 <- 0.9 * replace(numeric(D), k, 1) + 0.1 / D
      p[nm == mi$zeta] <- tox(w0)
    }
    lo <- rep(-Inf, length(p)); up <- rep(Inf, length(p))
    if (s > 0) up[j] <- 0 else lo[j] <- 0
    o <- stats::nlminb(p, obj$fn, obj$gr, lower = lo, upper = up,
                       control = list(eval.max = 1000, iter.max = 1000))
    n <- n + o$evaluations[[1]]
    o2 <- nl(obj, stats::setNames(o$par, nm))
    n <- n + o2$evaluations[[1]]
    if (is.finite(o2$objective) && o2$objective < fbest - 1e-10) {
      par <- o2$par
      fbest <- o2$objective
    }
  }
  list(par = par, objective = fbest, evals = n)
}
rows <- list()
for (s in seeds) {
  d <- mk(s)
  ex <- exact_max(d)
  f <- suppressWarnings(frm(ls ~ mo(income) * age, data = d))
  p0 <- f$opt$par
  f0 <- f$opt$objective
  fl <- if (identical(Sys.getenv("OPTIMA_FLIP"), "cone")) {
    cone(f, p0, f0)
  } else flip(f, p0, f0)
  ve <- if (param == "soft") vertex_escape(f$obj, p0) else
    list(objective = f0, evals = 0L, par = p0)
  vf <- if (param == "soft") {
    a <- flip(f, ve$par, ve$objective)
    b <- vertex_escape(f$obj, a$par)
    list(objective = b$objective, evals = ve$evals + a$evals + b$evals)
  } else fl
  rows[[length(rows) + 1L]] <- data.frame(
    seed = s, ll_exact = ex$ll, gap_fit = ex$ll + f0,
    gap_flip = ex$ll + fl$objective, gap_vertex = ex$ll + ve$objective,
    gap_both = ex$ll + vf$objective, ev_fit = f$opt$evals,
    ev_flip = fl$evals, ev_vertex = ve$evals, ev_both = vf$evals,
    autoscaled = f$autoscaled)
}
X <- do.call(rbind, rows)
utils::write.table(X, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("seeds", nrow(X), "\n")
