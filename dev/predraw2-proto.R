# Lane predraw2, prototype on an unchanged build: for each fit, predict()'s
# own parameter draws (the 0.82.0 drawer, rebuilt here from its internals),
# the drop in log-likelihood at each draw against the drop the draw law
# implies there, and the importance weights likelihood / proposal over
# 1000 draws with their Pareto k.
#
# e = drop - m / 2 - sum(log g'(x) - log g'(0)), the excess drop of one
# draw: m is the draw's Mahalanobis distance under the normal law the
# draws come from, g the side law's map on a narrowed coordinate. e is 0
# at every draw when the log-likelihood is the quadratic V describes.
#
# Usage: Rscript dev/predraw2-proto.R <lib or "base"> <group> <out.tsv>
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
source("dev/predraw2-cases.R")
group <- args[2]
out <- args[3]

pieces <- function(object, ndraws, seed) {
  ds <- ns$fit_draw_space(object)
  V <- ds$V
  if (!is.null(ds$jc) && !isTRUE(ds$boundary_only) &&
        !isTRUE(ds$mo_const_only)) V <- NULL
  if (any(ds$lost %in% "limit") && !ns$draw_limit_quiet(object, ds$lost)) {
    V <- NULL
  }
  if (is.null(ds$root)) V <- ns$draw_limit_restore(V, ds, object)
  L <- if (!is.null(ds$root) && !is.null(V)) {
    t(ds$root)
  } else if (is.null(V) || !all(is.finite(V))) {
    NULL
  } else {
    tryCatch(chol(V + diag(0, nrow(V))), error = function(e) {
      if (is.null(ds$jc)) return(NULL)
      e <- eigen(V, symmetric = TRUE)
      sqrt(pmax(e$values, 0)) * t(e$vectors)
    })
  }
  if (is.null(L)) return(NULL)
  map <- ds$map
  v0 <- ns$fit_outer_vector(object, map)
  p <- length(v0)
  set.seed(seed)
  Z <- matrix(stats::rnorm(p * ndraws), p, ndraws)
  law <- suppressMessages(ns$draw_edge_law(object, map, V, colSums(L^2),
                                           !is.null(ds$root),
                                           rownames(V) %||% map$names))
  Lf <- law$Lc %||% L
  X <- crossprod(Lf, Z)
  D <- v0 + if (is.null(law)) X else ns$draw_edge_apply(law, L, Z)
  held <- ns$latres_held_pos(object, map)
  if (length(held)) D[held, ] <- v0[held]
  list(D = D, Z = Z, X = X, Lf = Lf, law = law, v0 = v0, held = held,
       off = p - length(object$opt$par))
}

# slope of draw_edge_map() at deviation x of a coordinate with sd s, law v
map_slope <- function(x, s, v) {
  if (length(v) == 1L) return(rep(v, length(x)))
  z <- abs(x) / s
  k <- pmin(floor(z), 3L)
  vv <- c(0, v)
  sl_g <- ifelse(z > 4, max(0, v[4] - v[3]), vv[k + 2L] - vv[k + 1L])
  gz <- stats::approx(0:4, vv, xout = pmin(z, 4))$y
  gz[z > 4] <- v[4] + (z[z > 4] - 4) * max(0, v[4] - v[3])
  piece <- cbind(z, z * v[2] / 2, gz)
  w <- max.col(-piece, ties.method = "first")
  cbind(1, v[2] / 2, sl_g)[cbind(seq_along(z), w)]
}

lw_of <- function(object, pc, cols) {
  obj <- object$obj
  p0 <- object$opt$par
  np <- length(p0)
  pos <- pc$off + seq_len(np)
  saved <- ns$obj_state_save(obj)
  on.exit(ns$obj_state_restore(obj, saved), add = TRUE)
  f0 <- obj$fn(p0)
  drop <- vapply(cols, function(s) {
    v <- tryCatch(suppressWarnings(obj$fn(pc$D[pos, s])),
                  error = function(e) NA_real_)
    if (is.finite(v)) v - f0 else Inf
  }, 0)
  law <- pc$law
  edge <- law$edge %||% integer(0)
  sv <- colSums(pc$Lf^2)
  R <- setdiff(pos[sv[pos] > 0], c(pc$held, edge))
  qf <- qr(pc$Lf[, R, drop = FALSE])
  rk <- qf$rank
  m <- colSums(qr.qty(qf, pc$Z[, cols, drop = FALSE])[seq_len(rk), ,
                                                       drop = FALSE]^2)
  lj <- numeric(length(cols))
  stratum <- rep("", length(cols))
  if (!is.null(law)) {
    mapped <- which(!vapply(law$lo, function(v) length(v) == 1L && v == 1,
                            TRUE) |
                      !vapply(law$hi, function(v) length(v) == 1L && v == 1,
                              TRUE))
    for (j in intersect(mapped, c(R, edge))) {
      x <- pc$X[j, cols]
      neg <- x < 0
      sl <- numeric(length(x))
      sl[neg] <- map_slope(x[neg], law$s[j], law$lo[[j]])
      sl[!neg] <- map_slope(x[!neg], law$s[j], law$hi[[j]])
      s0 <- ifelse(neg, map_slope(-1e-12, law$s[j], law$lo[[j]]),
                   map_slope(1e-12, law$s[j], law$hi[[j]]))
      moving <- sl > 0
      if (j %in% edge) {
        stratum <- paste0(stratum, ifelse(moving, "s", "f"))
        m[moving] <- m[moving] + x[moving]^2 / sv[j]
      }
      lj[moving] <- lj[moving] + log(sl[moving]) - log(s0[moving])
    }
  }
  list(drop = drop, m = m, e = drop - m / 2 - lj, lw = -drop + m / 2 + lj,
       stratum = stratum, rank = rk + length(edge))
}

psis_k <- function(lw, stratum) {
  ks <- tapply(seq_along(lw), stratum, function(i) {
    if (length(i) < 30L) return(NA_real_)
    r <- suppressWarnings(loo::psis(lw[i], r_eff = 1))
    c(k = r$diagnostics$pareto_k, ess = r$diagnostics$n_eff / length(i))
  })
  ks <- do.call(rbind, ks)
  c(k = max(ks[, 1]), ess = min(ks[, 2]))
}

rows <- list()
for (id in case_ids(group)) {
  t0 <- proc.time()[[3]]
  cf <- tryCatch(case_fit(id), error = function(e) NULL)
  fit <- cf$fit
  if (is.null(fit)) {
    cat(id, "fit failed\n")
    next
  }
  tf <- proc.time()[[3]] - t0
  pc <- tryCatch(pieces(fit, 1000L, 20261113L), error = function(e) {
    cat(id, "pieces:", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(pc)) {
    rows[[id]] <- data.frame(id = id, law = "none", p = NA, med20 = NA,
                             max20 = NA, n2_20 = NA, med = NA, q90 = NA,
                             mean_drop = NA, chi_mean = NA, k = NA, ess = NA,
                             t_fit = tf, t20 = NA, t1000 = NA)
    next
  }
  t1 <- proc.time()[[3]]
  a20 <- lw_of(fit, pc, 1:20)
  t2 <- proc.time()[[3]]
  a <- lw_of(fit, pc, 1:1000)
  t3 <- proc.time()[[3]]
  ec <- function(r) {
    # within a stratum, e is anchored at the estimate
    r$e
  }
  e20 <- ec(a20)
  e <- ec(a)
  kk <- tryCatch(psis_k(a$lw, a$stratum), error = function(err) c(NA, NA))
  lawtag <- if (is.null(pc$law)) "normal" else
    paste0("edge", length(pc$law$edge), "nar", length(pc$law$narrowed))
  rows[[id]] <- data.frame(
    id = id, law = lawtag, p = a$rank, med20 = stats::median(e20),
    max20 = max(e20), n2_20 = sum(e20 > 2), med = stats::median(e),
    q90 = unname(stats::quantile(e, 0.9)), mean_drop = mean(a$drop),
    chi_mean = a$rank / 2, k = kk[[1]], ess = kk[[2]], t_fit = tf,
    t20 = t2 - t1, t1000 = t3 - t2)
  r <- rows[[id]]
  cat(sprintf("%-38s %-10s p %2d med20 %8.3f max20 %8.3f n2 %2d | med %8.3f q90 %8.3f k %6.2f ess %5.3f | t %.2f %.2f\n",
              id, r$law, r$p, r$med20, r$max20, r$n2_20, r$med, r$q90, r$k,
              r$ess, r$t20, r$t1000))
}
res <- do.call(rbind, rows)
utils::write.table(res, out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("rows", nrow(res), "\n")
