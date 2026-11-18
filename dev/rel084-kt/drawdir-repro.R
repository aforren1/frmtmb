# rel084-kt copy of dev/drawdir-repro.R: rellib-r22 removed from .libPaths() (DRAWDIR_LIB, then user lib). No other change.
# Lane drawdir: a copy of dev/reviews/predraw2-short-a.R (the short check of
# predraw2) with its library from DRAWDIR_LIB (empty: rellib-r22 alone).
# Run from the worktree: DRAWDIR_LIB=... Rscript dev/drawdir-repro.R
# Short check of predraw2's fix round, item 1(a): does a large count of
# well-identified parameters hide one weakly identified direction from
# the two-condition check (draw_check_verdict())?
# Part S (exact for block-separable objectives): the check's own e and m
#   from a fired mixture fit, plus a block of K well-identified
#   coordinates: (i) exactly quadratic (e 0, m chi-square K), (ii) the
#   review's poisson y ~ 0 + f design (its own e and m).
# Part R (real fits): the same mixture data with K noise covariates,
#   y ~ x + z1..zK in both components; verdict, and the predictive means'
#   error against the truth, unfired (lane) and with the check forced on.
# Usage (from the lane tree): Rscript dev/reviews/predraw2-short-a.R
local({
lib <- Sys.getenv("DRAWDIR_LIB")
.libPaths(c(if (nzchar(lib)) lib,
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), " tol", ns$draw_check_tol,
    "rel", ns$draw_check_rel, "\n")
q <- function(expr) suppressWarnings(suppressMessages(expr))
n <- 300L
source("dev/ppdraw-sweep-gen.R", local = TRUE)
real <- ns$draw_reweight
cap <- new.env()
grab <- function(object, dw) {
  dw$off <- length(dw$v0) - length(object$opt$par)
  cap$dw <- dw
  dw$D
}
em <- function(fit) {
  cap$dw <- NULL
  assignInNamespace("draw_reweight", grab, ns = "frmtmb")
  on.exit(assignInNamespace("draw_reweight", real, ns = "frmtmb"))
  set.seed(1)
  q(ns$predict_par_drawer(fit, TRUE, 2L))
  dw <- cap$dw
  if (is.null(dw)) return(NULL)
  Z0 <- ns$draw_check_z(length(dw$v0), 2L * ns$draw_check_n)
  D0 <- dw$v0 + if (is.null(dw$law)) crossprod(dw$L, Z0) else
    ns$draw_edge_apply(dw$law, dw$L, Z0)
  if (length(dw$held)) D0[dw$held, ] <- dw$v0[dw$held]
  r <- ns$draw_logratio(fit, dw, D0, Z0)
  chk <- ns$draw_check(fit, dw)
  list(e = r$e, m = r$m, p = length(fit$opt$par), chk = chk)
}
verdict <- function(e, m) {
  f <- seq_len(ns$draw_check_n)
  v <- ns$draw_check_verdict(e[f], m[f])
  v$n <- length(f)
  if (v$near) {
    v <- ns$draw_check_verdict(e, m)
    v$n <- length(e)
  }
  v
}
fmt <- function(v) sprintf("%s ex %.3g rel %.3g n %d",
                           if (v$fire) "FIRE" else "----", v$excess, v$rel, v$n)

## Part S
cat("\n== Part S: fired mixture fits + a separable block of K coords ==\n")
fired <- list(mixture_gaussian = c(20, 23, 25, 29, 36, 40),
              mixture_poisson = c(5, 11, 13, 18, 20, 24, 25, 26, 31, 33, 36, 39))
# the review's K design, seed 1 (no level without counts)
set.seed(1)
dk <- data.frame(f = factor(rep(1:200, each = 5)))
dk$y <- stats::rpois(1000, exp(log(2) + stats::rnorm(200, 0, 0.3))[dk$f])
fk <- q(frm(bf(y ~ 0 + f), family = poisson(), data = dk))
ek <- em(fk)
cat(sprintf("K200 poisson block alone: p %d | %s | cached fire %s ex %.4g\n",
            ek$p, fmt(verdict(ek$e, ek$m)), ek$chk$fire, ek$chk$excess))
Ks <- c(0, 4, 8, 16, 32, 64, 128, 256)
mixcache <- list()
for (fam in names(fired)) for (sd in fired[[fam]]) {
  d <- gen(fam, "boundary", sd)
  s <- spec(fam)
  fit <- tryCatch(q(frm(s[[1]], family = s[[2]], data = d)),
                  error = function(e) NULL)
  if (is.null(fit)) { cat(fam, sd, "fit failed\n"); next }
  a <- em(fit)
  if (is.null(a)) { cat(fam, sd, "no draws\n"); next }
  mixcache[[paste(fam, sd)]] <- a
  v0 <- verdict(a$e, a$m)
  cat(sprintf("%-16s seed %2d p %d alone: %s (cached %s %.4g) | median m %.2f\n",
              fam, sd, a$p, fmt(v0), a$chk$fire, a$chk$excess,
              stats::median(a$m[seq_len(v0$n)])))
  out <- character(0)
  for (K in Ks) {
    set.seed(9000 + K)
    mk <- if (K == 0) numeric(40) else stats::rchisq(40, K)
    v <- verdict(a$e, a$m + mk)
    out <- c(out, sprintf("K%d %s %.3g/%.3g", K, if (v$fire) "F" else "-",
                          v$excess, v$rel))
  }
  cat("    quadratic block:", paste(out, collapse = " | "), "\n")
  v2 <- verdict(a$e + ek$e, a$m + ek$m)
  cat("    + poisson K200 block (p", a$p + ek$p, "):", fmt(v2), "\n")
}

## Part R
cat("\n== Part R: the same data with K noise covariates in both components ==\n")
tm_of <- function(fam, d) if (fam == "mixture_gaussian") 0.5 * d$x else
  exp(1 + 0.3 * d$x)
run <- function(fit, force) {
  if (force) {
    chk <- fit$cache[["draw_check"]]
    chk$fire <- TRUE
    assign("draw_check", chk, envir = fit$cache)
  }
  said <- character(0)
  set.seed(1)
  p <- withCallingHandlers(
    predict(fit, ndraws = 1000, summary = FALSE),
    warning = function(w) {
      said <<- c(said, class(w)[1]); invokeRestart("muffleWarning")
    }, message = function(m) {
      said <<- c(said, class(m)[1]); invokeRestart("muffleMessage")
    })
  list(p = p, said = said)
}
for (fam in names(fired)) for (sd in fired[[fam]]) for (K in c(0, 5, 10, 20)) {
  d <- gen(fam, "boundary", sd)
  set.seed(5000 + sd)
  zs <- character(0)
  for (j in seq_len(K)) {
    zs <- c(zs, paste0("z", j))
    d[[paste0("z", j)]] <- stats::rnorm(n)
  }
  fo <- stats::as.formula(paste("y ~", paste(c("x", zs), collapse = " + ")))
  fit <- tryCatch(q(frm(bf(fo), family = spec(fam)[[2]], data = d)),
                  error = function(e) NULL)
  if (is.null(fit)) { cat(fam, sd, "K", K, "fit failed\n"); next }
  r1 <- tryCatch(run(fit, FALSE), error = function(e) NULL)
  chk <- fit$cache[["draw_check"]]
  if (is.null(r1) || is.null(chk) || is.na(chk$excess)) {
    cat(sprintf("%-16s seed %2d K %2d no verdict\n", fam, sd, K)); next
  }
  tm <- tm_of(fam, d)
  fe <- q(fitted(fit))[, "Estimate"]
  mae1 <- mean(abs(colMeans(r1$p) - tm))
  mu1 <- rowMeans(r1$p)
  z1 <- (mean(mu1) - mean(fe)) / (stats::sd(mu1) / sqrt(1000))
  line <- sprintf("%-16s seed %2d K %2d p %2d %s ex %7.3f rel %6.3f n %d | mae lane %.4f z %6.2f",
                  fam, sd, K, length(fit$opt$par),
                  if (chk$fire) "FIRE" else "----", chk$excess, chk$rel, chk$n,
                  mae1, z1)
  if (!chk$fire) {
    if (exists("draw_check_said", envir = fit$cache, inherits = FALSE)) {
      rm("draw_check_said", envir = fit$cache)
    }
    r2 <- tryCatch(run(fit, TRUE), error = function(e) NULL)
    if (!is.null(r2)) {
      mu2 <- rowMeans(r2$p)
      line <- paste0(line, sprintf(" | forced mae %.4f z %6.2f said %s",
                                   mean(abs(colMeans(r2$p) - tm)),
                                   (mean(mu2) - mean(fe)) /
                                     (stats::sd(mu2) / sqrt(1000)),
                                   paste(unique(r2$said), collapse = ",")))
    }
  }
  cat(line, "\n")
}
})
