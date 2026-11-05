# Fix round 1, item 6: whose account of M2 is right? se_tier3()'s rows read at y * 1000 with the old steps and
# with the shipped steps, seed 14, noise 1e-4. Then se_row_curved() gated on curvature in the unit, on the
# five ordinary fits of review-ncgp-r2-claim.R and on the rows it exists for.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
patch <- function(nm, f) { environment(f) <- ns; unlockBinding(nm, ns); assign(nm, f, envir = ns); lockBinding(nm, ns) }
o_t3 <- ns$se_tier3; o_steps <- ns$outer_fd_steps; o_cur <- ns$se_row_curved
cat("se_empty_row", ns$se_empty_row, "se_noise_mult", ns$se_noise_mult, "\n")
spy <- function(fit, H, E, p, exact = FALSE, held = NULL, edge = logical(0), lazy = FALSE) {
  u <- fit$par_units %||% rep(1, NROW(H))
  Hu <- abs(H * outer(u, u)); Eu <- E * outer(u, u)
  cat(sprintf("  rows units %s\n  row max Hu %s\n  floor (mult * max Eu) %s  exact %s\n", paste(signif(u, 3), collapse = " "), paste(signif(apply(Hu, 1, max), 3), collapse = " "), paste(signif(se_noise_mult * apply(Eu, 1, max), 3), collapse = " "), exact))
  o_t3(fit, H, E, p, exact = exact, held = held, edge = edge, lazy = lazy)
}
old_steps <- o_steps
b <- deparse(body(o_steps)); b2 <- sub("if (any(vapply", "if (FALSE && any(vapply", paste(b, collapse = "\n"), fixed = TRUE); stopifnot(b2 != paste(b, collapse = "\n"))
body(old_steps) <- parse(text = b2)[[1]]
set.seed(14); x <- sort(runif(40, 0, 8)); y0 <- sin(x) + 1e-4 * rnorm(40)
patch("se_tier3", spy)
for (c0 in c(1, 1000)) for (arm in c("shipped", "oldsteps", "oldsteps,norow")) {
  patch("outer_fd_steps", if (arm == "shipped") o_steps else old_steps)
  patch("se_row_curved", if (arm == "oldsteps,norow") function(j, fit, p, u) FALSE else o_cur)
  cat(sprintf("RM2 c %g arm %s\n", c0, arm))
  f <- suppressWarnings(frm(y ~ gp(x), data.frame(y = y0 * c0, x = x)))
  cat(sprintf("  SE/c intercept %.4f; all SE %s\n", sqrt(as.matrix(vcov(f, full = TRUE))[1, 1]) / c0, paste(signif(sqrt(diag(as.matrix(vcov(f, full = TRUE)))), 4), collapse = " ")))
}
patch("se_tier3", o_t3); patch("outer_fd_steps", o_steps); patch("se_row_curved", o_cur)
cat("== a gate on curvature in the unit\n")
gated <- function(j, fit, p, u) {
  obj <- fit$obj
  saved <- obj_state_save(obj); on.exit(obj_state_restore(obj, saved), add = TRUE)
  h <- 1e-3 * u[j]
  d <- function(s) { q <- p; q[j] <- p[j] + s * h; up <- drop(obj$gr(q))[j]; q[j] <- p[j] - s * h; (up - drop(obj$gr(q))[j]) / (2 * s * h) }
  d1 <- tryCatch(d(1), error = function(e) NA_real_); d2 <- tryCatch(d(4), error = function(e) NA_real_)
  isTRUE(d1 > 0 && d2 > 0 && abs(d1 / d2 - 1) <= 0.1 && d1 * u[j]^2 > se_empty_row)
}
for (cs in list(list(seed = 2, sg = 1, c0 = 1), list(seed = 2, sg = 1, c0 = 1000), list(seed = 4, sg = 0.1, c0 = 1), list(seed = 4, sg = 0.1, c0 = 3), list(seed = 4, sg = 0.1, c0 = 1000))) {
  set.seed(60000 + cs$seed); n <- 50L
  x <- sort(round(runif(n, 0, 8), 4)); z <- rnorm(n); g <- factor(sample(letters[1:5], n, TRUE))
  y <- cs$c0 * (sin(2 * x) + 0.3 * z + cs$sg * rnorm(n)); d <- data.frame(x, z, y)
  fml <- if (cs$seed == 2) y ~ z + gp(x) else y ~ gp(x, z)
  w0 <- character(); f0 <- withCallingHandlers(frm(fml, d), warning = function(c) { w0 <<- c(w0, conditionMessage(c)); invokeRestart("muffleWarning") })
  patch("se_row_curved", gated)
  w1 <- character(); f1 <- withCallingHandlers(frm(fml, d), warning = function(c) { w1 <<- c(w1, conditionMessage(c)); invokeRestart("muffleWarning") })
  patch("se_row_curved", o_cur)
  cat(sprintf("RM2G s=%g c=%g shipped SE %s w %d | gated SE %s w %d\n", cs$sg, cs$c0, paste(signif(sqrt(diag(as.matrix(vcov(f0, full = TRUE)))), 3), collapse = " "), length(w0), paste(signif(sqrt(diag(as.matrix(vcov(f1, full = TRUE)))), 3), collapse = " "), length(w1)))
  flush.console()
}
for (c0 in c(1000, 1e6, 1e9)) {
  patch("se_row_curved", gated)
  f <- suppressWarnings(frm(y ~ gp(x), data.frame(y = y0 * c0, x = set_x <- sort({ set.seed(14); runif(40, 0, 8) }))))
  patch("se_row_curved", o_cur)
  cat(sprintf("RM2G units c %g SE/c %.4f\n", c0, sqrt(as.matrix(vcov(f, full = TRUE))[1, 1]) / c0))
}
cat("DONE\n")
