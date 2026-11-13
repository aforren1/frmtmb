# Shared helpers of the ppdraw review. Sourced after .libPaths() is set.
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
lane_has_scan <- exists("draw_edge_scan", envir = asNamespace("frmtmb"))
cat("lane code present:", lane_has_scan, "\n")

q <- function(expr) suppressMessages(suppressWarnings(expr))

# smaller |change| of the objective at +-2 draw-law sds, per parameter,
# as draw_edge_scan() asks it
flat_sides <- function(fit) {
  ds <- frmtmb:::fit_draw_space(fit)
  if (is.null(ds$V) || !all(is.finite(ds$V))) return(NULL)
  p <- fit$opt$par
  off <- length(ds$map$names) - length(p)
  sd <- sqrt(diag(ds$V))
  f0 <- fit$obj$fn(p)
  out <- t(vapply(seq_along(p), function(j) {
    s <- sd[off + j]
    if (!is.finite(s) || s <= 0) return(c(p[j], s, NA, NA))
    dr <- vapply(c(-1, 1), function(k) {
      qq <- p
      qq[j] <- p[j] + 2 * k * s
      tryCatch(suppressWarnings(fit$obj$fn(qq)) - f0,
               error = function(e) NA_real_)
    }, 0)
    c(p[j], s, dr)
  }, numeric(4)))
  dimnames(out) <- list(rownames(ds$V)[off + seq_along(p)],
                        c("est", "sd", "minus2sd", "plus2sd"))
  fit$obj$fn(p)
  out
}

# predict()'s mean against fitted()'s, in Monte Carlo standard errors of
# the draws' mean; held parameters named by the lane's message
pz <- function(fit, ndraws = 400L, seed = 1, newdata = NULL) {
  held <- character(0)
  warns <- character(0)
  set.seed(seed)
  p <- withCallingHandlers(
    predict(fit, newdata = newdata, ndraws = ndraws, summary = FALSE),
    frmtmb_draw_held = function(m) {
      held <<- c(held, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    message = function(m) invokeRestart("muffleMessage"),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  m <- rowMeans(p, na.rm = TRUE)
  fm <- mean(fitted(fit, newdata = newdata)[, "Estimate"])
  list(z = (mean(m) - fm) / (stats::sd(m) / sqrt(ndraws)),
       pmean = mean(m), fmean = fm, na = mean(is.na(p)),
       held = held, warns = unique(warns), p = p)
}
