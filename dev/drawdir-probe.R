# Lane drawdir: the per-direction reading of the draw check, measured
# before it is built. For each fit, through the build's own draw law (the
# normal law of vcov() and the side law), the objective is evaluated at
# +-c along each of a set of directions of the deviates and the excess
# e = drop - m / 2 - log Jacobian is read (draw_logratio(), the summed
# check's own per-draw quantity). Three sets of directions:
#   eig  the eigen-directions of the draw law's covariance (left singular
#        vectors of its factor)
#   cor  the eigen-directions of its correlation matrix (scale free)
#   par  one per parameter: the parameter moved c of its own standard
#        deviations and the others by their regression on it (column j
#        of the factor), so the law's drop is c^2 / 2
# Also the summed check's verdict and the quantiles of its 40 test draws'
# excesses (for the two-sided reading of backlog 5644).
# Usage: Rscript dev/drawdir-probe.R <lib or "base"> <out.tsv> <group>
#        [id filter regex]
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), "\n")
out <- args[2]
source("dev/drawdir-cases.R")
ids <- dc_ids(args[3])
if (length(args) > 3) ids <- grep(args[4], ids, value = TRUE)
cs <- as.numeric(strsplit(Sys.getenv("DRAWDIR_CS", "1.5,2,3"), ",")[[1]])
real <- ns$draw_reweight
cap <- new.env()
grab <- function(object, dw) {
  dw$off <- length(dw$v0) - length(object$opt$par)
  cap$dw <- dw
  dw$D
}
get_dw <- function(fit) {
  cap$dw <- NULL
  assignInNamespace("draw_reweight", grab, ns = "frmtmb")
  on.exit(assignInNamespace("draw_reweight", real, ns = "frmtmb"))
  set.seed(1)
  fq(ns$predict_par_drawer(fit, TRUE, 2L))
  cap$dw
}
mk <- function(dw, Z) {
  D <- dw$v0 + if (is.null(dw$law)) crossprod(dw$L, Z) else
    ns$draw_edge_apply(dw$law, dw$L, Z)
  if (length(dw$held)) D[dw$held, ] <- dw$v0[dw$held]
  D
}
dirs <- function(dw, kind, R) {
  Lf <- dw$law$Lc %||% dw$L
  B <- Lf[, R, drop = FALSE]
  if (kind == "par") {
    return(sweep(B, 2L, sqrt(colSums(B^2)), "/"))
  }
  if (kind == "cor") B <- sweep(B, 2L, sqrt(colSums(B^2)), "/")
  s <- svd(B)
  s$u[, s$d > max(s$d) * 1e-10, drop = FALSE]
}
lab <- function(dw, U, k, R, names) {
  Lf <- dw$law$Lc %||% dw$L
  x <- drop(crossprod(Lf, U[, k]))[R] / sqrt(colSums(Lf[, R, drop = FALSE]^2))
  names[R][which.max(abs(x))]
}
rows <- list()
for (id in ids) {
  t0 <- proc.time()[[3]]
  r <- tryCatch(dc_fit(id), error = function(e) NULL)
  if (is.null(r) || is.null(r$fit)) {
    cat(id, "fit failed\n")
    next
  }
  fit <- r$fit
  t1 <- proc.time()[[3]]
  dw <- tryCatch(get_dw(fit), error = function(e) NULL)
  t2 <- proc.time()[[3]]
  if (is.null(dw)) {
    cat(id, "no draws (plug-in)\n")
    next
  }
  chk <- fq(ns$draw_check(fit, dw))
  t3 <- proc.time()[[3]]
  nm <- rownames(vcov(fit, full = TRUE)) %||% names(dw$v0)
  if (length(nm) != length(dw$v0)) nm <- paste0("v", seq_along(dw$v0))
  Z0 <- ns$draw_check_z(length(dw$v0), 2L * ns$draw_check_n)
  r0 <- ns$draw_logratio(fit, dw, mk(dw, Z0), Z0)
  np <- length(fit$opt$par)
  pos <- dw$off + seq_len(np)
  Lf <- dw$law$Lc %||% dw$L
  sv <- colSums(Lf^2)
  R <- setdiff(pos[sv[pos] > 0], c(dw$held, dw$law$edge))
  row <- data.frame(id = id, np = np, nR = length(R),
                    fire = isTRUE(chk$fire), excess = chk$excess,
                    rel = chk$rel, n = chk$n,
                    law = !is.null(dw$law),
                    t_fit = t1 - t0, t_law = t2 - t1, t_chk = t3 - t2)
  q <- stats::quantile(r0$e, c(0.1, 0.25, 0.5, 0.75, 0.9), names = FALSE)
  row[c("q10", "q25", "q50", "q75", "q90")] <- as.list(q)
  row$med40 <- stats::median(r0$e)
  for (kind in c("eig", "cor", "par")) {
    U <- if (length(R)) dirs(dw, kind, R) else matrix(0, length(dw$v0), 0)
    k <- ncol(U)
    for (cc in cs) {
      tag <- paste0(kind, cc)
      if (!k) {
        row[paste0(tag, c("_max", "_min", "_arg", "_argmin", "_t",
                          "_mrng"))] <- list(NA, NA, "", "", 0, NA)
        next
      }
      Z <- cbind(cc * U, -cc * U)
      ta <- proc.time()[[3]]
      rr <- ns$draw_logratio(fit, dw, mk(dw, Z), Z)
      tb <- proc.time()[[3]]
      e <- rr$e
      e[!is.finite(e)] <- Inf
      jm <- which.max(e)
      jn <- which.min(e)
      row[[paste0(tag, "_max")]] <- e[jm]
      row[[paste0(tag, "_min")]] <- e[jn]
      row[[paste0(tag, "_arg")]] <- paste0(if (jm > k) "-" else "+",
                                           lab(dw, U, (jm - 1) %% k + 1, R, nm))
      row[[paste0(tag, "_argmin")]] <- paste0(if (jn > k) "-" else "+",
                                              lab(dw, U, (jn - 1) %% k + 1, R,
                                                  nm))
      row[[paste0(tag, "_t")]] <- tb - ta
      row[[paste0(tag, "_mrng")]] <- diff(range(rr$m))
    }
  }
  rows[[length(rows) + 1L]] <- row
  cat(sprintf(paste("%-40s np %3d fire %-5s ex %7.3f rel %6.3f | q25 %6.2f",
                    "q75 %6.2f | eig2 %7.2f/%7.2f cor2 %7.2f/%7.2f par2",
                    "%7.2f/%7.2f (%s) | t %.1f/%.1f/%.1f/%.2f\n"),
              id, np, row$fire, row$excess, row$rel, row$q25, row$q75,
              row$eig2_max, row$eig2_min, row$cor2_max, row$cor2_min,
              row$par2_max, row$par2_min, row$par2_arg, row$t_fit,
              row$t_law, row$t_chk, row$par2_t))
}
all <- do.call(rbind, lapply(rows, function(x) x[names(rows[[1]])]))
utils::write.table(all, out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("DONE", nrow(all), "rows\n")
