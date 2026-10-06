# Lane fixes, punch round: the s() null-space column against brms
# 2.23.0's standata Xs, smooth by smooth, and the fit against the base
# build (a reparameterization: logLik and fitted values must agree).
#   Rscript dev/fixes-sx.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), " mgcv",
    as.character(packageVersion("mgcv")), "\n")
set.seed(1)
d <- mgcv::gamSim(eg = 6, n = 200, scale = 2, verbose = FALSE)
d$z <- runif(200)
d$g <- factor(sample(c("a", "b", "c"), 200, TRUE))
d$fac <- factor(d$fac)
cases <- list(
  "s(x1) + s(x2)" = list(y ~ s(x1) + s(x2)),
  "s(x1, bs = 'cr', k = 6)" = list(y ~ s(x1, bs = "cr", k = 6)),
  "s(x1, by = g) + g" = list(y ~ s(x1, by = g) + g),
  "s(x1, by = z)" = list(y ~ s(x1, by = z)),
  "s(x1, x2)" = list(y ~ s(x1, x2)),
  "te(x1, x2)" = list(y ~ te(x1, x2)),
  "t2(x1, x2)" = list(y ~ t2(x1, x2)),
  "sigma ~ s(x0)" = list(y ~ s(x1), sigma ~ s(x0)))
nd <- d[c(3, 50, 120), ]
nd$x1 <- nd$x1 + 0.01
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fo <- do.call(bf, cs)
  bfo <- do.call(brms::bf, cs)
  fit <- tryCatch(suppressWarnings(frm(fo, data = d)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-28s frm ERROR %s\n", nm, conditionMessage(fit)))
    next
  }
  sd <- brms::standata(bfo, data = d)
  out <- character()
  for (lp in fit$frame$linpreds) {
    sfx <- if (identical(lp$dpar, "mu")) "" else paste0("_", lp$dpar)
    Xs <- sd[[paste0("Xs", sfx)]]
    fx <- grep("[.]fx[0-9]+$", colnames(lp$X), value = TRUE)
    if (!length(fx) && (is.null(Xs) || !ncol(Xs))) next
    a <- unname(lp$X[, fx, drop = FALSE])
    b <- unname(Xs)
    same <- identical(dim(a), dim(b)) && max(abs(a - b)) == 0
    out <- c(out, sprintf("%s fixed cols %d/%d bitwise %s maxdiff %.2e",
                          lp$dpar, ncol(a), NCOL(b), same,
                          if (identical(dim(a), dim(b))) max(abs(a - b))
                          else NA))
  }
  ll <- as.numeric(logLik(fit))
  p <- tryCatch(as.numeric(fitted(fit, newdata = nd)[, "Estimate"]),
                error = function(e) rep(NA, 3))
  cat(sprintf("%-28s %s | logLik %.10f | newdata fitted %s\n", nm,
              paste(out, collapse = "; "), ll,
              paste(format(p, digits = 10), collapse = " ")))
  fe <- fixef(fit)
  sx <- grep("^(sigma_)?s.*_[0-9]+$", rownames(fe), value = TRUE)
  if (length(sx)) cat("    ", paste(sprintf("%s=%.6f", sx, fe[sx, 1]),
                                    collapse = " "), "\n")
}

cat("\n## the wiggly part and the newdata path, s(x1) + s(x2) and by = g\n")
for (fo in list(y ~ s(x1) + s(x2), y ~ s(x1, by = g) + g)) {
  fit <- frm(bf(fo), data = d)
  sd <- brms::standata(brms::bf(fo), data = d)
  zs <- grep("^Zs_", names(sd), value = TRUE)
  # each smooth basis is a block of its own; its columns of the
  # predictor's Z are the wiggly part
  bks <- Filter(function(bk) identical(bk[["covstruct"]], "smooth"),
                fit$frame$re_blocks)
  Zl <- fit$frame$linpreds[[1]]$Z
  zf <- lapply(bks, function(bk) {
    as.matrix(Zl[, bk[["b_idx"]], drop = FALSE])
  })
  cat("  blocks:", length(bks), " brms Zs:", length(zs), "\n")
  for (k in seq_along(zs)) {
    a <- unname(zf[[k]]); b <- unname(as.matrix(sd[[zs[k]]]))
    cat(sprintf("  %s: %s bitwise %s; Z Z' rel %.2e\n", deparse(fo), zs[k],
                identical(dim(a), dim(b)) && max(abs(a - b)) == 0,
                max(abs(tcrossprod(a) - tcrossprod(b))) /
                  max(abs(tcrossprod(b)))))
  }
  a <- as.numeric(fitted(fit, newdata = d[1:8, ])[, "Estimate"])
  b <- as.numeric(fitted(fit)[1:8, "Estimate"])
  cat(sprintf("  newdata at the fitted rows vs in sample: max rel %.2e\n",
              max(abs(a - b) / abs(b))))
}
