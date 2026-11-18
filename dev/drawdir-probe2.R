# Lane drawdir: per-parameter detail of the per-direction reading. For
# each parameter j of the objective (not held, not on a flat side):
#   A  its direction of the draw law: j moved c of its own standard
#      deviations and the others by their regression on it (column j of
#      the law's factor), through the side law, read by draw_logratio()
#      (excess eA = drop - c^2 / 2 - log Jacobian);
#   B  its own slice, the others at the estimates: j moved c of its
#      conditional standard deviations (1 / sqrt of the precision's
#      diagonal), excess eB = drop - c^2 / 2, no side law. B is what
#      draw_edge_law() reads (at c = 1 on each side).
# Also each direction's angle to its coordinate: s_j / cs_j (1 where j is
# uncorrelated with the others).
# Usage: Rscript dev/drawdir-probe2.R <lib or "base"> <out.tsv> <group>
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
cs <- as.numeric(strsplit(Sys.getenv("DRAWDIR_CS", "1,1.5,2"), ",")[[1]])
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
rows <- list()
for (id in ids) {
  r <- tryCatch(dc_fit(id), error = function(e) NULL)
  if (is.null(r) || is.null(r$fit)) {
    cat(id, "fit failed\n")
    next
  }
  fit <- r$fit
  dw <- tryCatch(get_dw(fit), error = function(e) NULL)
  if (is.null(dw)) {
    cat(id, "no draws (plug-in)\n")
    next
  }
  chk <- fq(ns$draw_check(fit, dw))
  nm <- rownames(vcov(fit, full = TRUE)) %||% names(dw$v0)
  if (length(nm) != length(dw$v0)) nm <- paste0("v", seq_along(dw$v0))
  np <- length(fit$opt$par)
  pos <- dw$off + seq_len(np)
  Lf <- dw$law$Lc %||% dw$L
  sv <- colSums(Lf^2)
  R <- setdiff(pos[sv[pos] > 0], c(dw$held, dw$law$edge))
  if (!length(R)) next
  B <- Lf[, R, drop = FALSE]
  s <- sqrt(colSums(B^2))
  U <- sweep(B, 2L, s, "/")
  # conditional sd of each R coordinate given the other R coordinates
  VR <- crossprod(B)
  csd <- tryCatch(1 / sqrt(diag(solve(VR))), error = function(e) rep(NA, length(R)))
  obj <- fit$obj
  saved <- ns$obj_state_save(obj)
  f0 <- obj$fn(fit$opt$par)
  at <- function(v) {
    ns$obj_state_restore(obj, saved)
    f1 <- tryCatch(suppressWarnings(obj$fn(v[pos])), error = function(e) NA)
    if (is.finite(f1)) f1 - f0 else Inf
  }
  for (cc in cs) {
    Z <- cbind(cc * U, -cc * U)
    rr <- ns$draw_logratio(fit, dw, mk(dw, Z), Z)
    k <- length(R)
    eB <- vapply(c(1, -1), function(sg) vapply(seq_len(k), function(i) {
      v <- dw$v0
      v[R[i]] <- v[R[i]] + sg * cc * csd[i]
      at(v) - cc^2 / 2
    }, 0), numeric(k))
    rows[[length(rows) + 1L]] <- data.frame(
      id = id, np = np, fire = isTRUE(chk$fire), excess = chk$excess,
      law = !is.null(dw$law), c = cc, par = nm[R], ratio = s / csd,
      eAp = rr$e[seq_len(k)], eAm = rr$e[k + seq_len(k)],
      eBp = eB[, 1], eBm = eB[, 2])
  }
  ns$obj_state_restore(obj, saved)
  x <- rows[[length(rows)]]
  j <- which.max(pmax(x$eAp, x$eAm))
  cat(sprintf("%-40s np %3d fire %-5s | c %.1f max eA %8.3f (%s, s/cs %.2f) eB %8.3f | max eA-eB %8.3f\n",
              id, np, chk$fire, x$c[1], max(x$eAp, x$eAm), x$par[j],
              x$ratio[j], max(x$eBp[j], x$eBm[j]),
              max(c(x$eAp - pmax(x$eBp, 0), x$eAm - pmax(x$eBm, 0)))))
}
utils::write.table(do.call(rbind, rows), out, sep = "\t", row.names = FALSE,
                   quote = FALSE)
cat("DONE\n")
