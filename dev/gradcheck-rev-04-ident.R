## Reviewer, claim 2: bitwise identity of the FIT between the two
## builds, on model shapes the worker's 720-fit grid does not contain.
## usage: Rscript gradcheck-rev-04-ident.R <core-lib> <out.rds>
## Writes one row per model: coef, logLik, convergence, max|grad|, the
## warning strings, and the objective at the optimum.
a <- commandArgs(TRUE)
LIB <- a[1]; OUT <- a[2]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), "\n")
`%||%` <- function(x, y) if (is.null(x)) y else x

gw <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                             warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
  })
  list(fit = val, w = w)
}
row_of <- function(tag, r) {
  f <- r$fit
  if (inherits(f, "condition")) {
    return(list(tag = tag, err = conditionMessage(f)))
  }
  g <- try(drop(f$obj$gr(f$opt$par)), silent = TRUE)
  list(tag = tag, err = NA_character_,
       coef = coef(f), par = f$opt$par,
       logLik = as.numeric(logLik(f)),
       objective = as.numeric(f$obj$fn(f$opt$par)),
       opt_obj = f$opt$objective,
       convergence = f$opt$convergence,
       gmax = if (inherits(g, "try-error")) NA_real_ else max(abs(g)),
       grad = if (inherits(g, "try-error")) NA_real_ else g,
       warns = sort(r$w))
}

models <- list()

models$mixture <- function() {
  set.seed(7001)
  n <- 900
  z <- rbinom(n, 1, 0.4)
  dd <- data.frame(x = rnorm(n))
  dd$y <- ifelse(z == 1, rnorm(n, 5, 1), rnorm(n, 0, 1)) + 0.5 * dd$x
  frm(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = dd)
}
models$mixture_collapsed <- function() {
  # both components chase the same mode, so one collapses
  set.seed(7002)
  n <- 800
  dd <- data.frame(x = rnorm(n))
  dd$y <- rnorm(n, 1 + 0.5 * dd$x, 1)
  frm(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = dd)
}
models$nonlinear <- function() {
  set.seed(7003)
  n <- 700
  dd <- data.frame(t = runif(n, 0, 5))
  dd$y <- rnorm(n, 3 * exp(-0.8 * dd$t), 0.2)
  frm(bf(y ~ a * exp(-b * t), a ~ 1, b ~ 1, nl = TRUE),
      family = gaussian(), data = dd,
      start = list(a = 2, b = 1))
}
models$rescor <- function() {
  set.seed(7004)
  n <- 900
  dd <- data.frame(x = rnorm(n))
  e <- matrix(rnorm(2 * n), n, 2) %*% chol(matrix(c(1, 0.6, 0.6, 1), 2))
  dd$y1 <- 1 + 2 * dd$x + e[, 1]
  dd$y2 <- -1 + 0.5 * dd$x + e[, 2]
  frm(bf(y1 ~ x) + gaussian() + bf(y2 ~ x) + gaussian() + set_rescor(TRUE),
      data = dd)
}
models$ordinal_thres_gr <- function() {
  set.seed(7005)
  n <- 1200
  dd <- data.frame(x = rnorm(n), grp = factor(rep(c("a", "b"), n / 2)))
  eta <- 0.9 * dd$x + ifelse(dd$grp == "b", 0.4, 0)
  dd$yo <- cut(eta + rlogis(n), breaks = c(-Inf, -0.8, 0.6, Inf),
               labels = FALSE)
  frm(bf(yo | thres(gr = grp) ~ x), family = cumulative(), data = dd)
}
models$mi <- function() {
  set.seed(7006)
  n <- 600
  dd <- data.frame(x = rnorm(n))
  dd$y <- rnorm(n, 1 + 2 * dd$x, 1)
  dd$xm <- dd$x
  dd$xm[sample(n, 60)] <- NA
  frm(bf(y ~ mi(xm)) + gaussian() + bf(xm | mi() ~ 1) + gaussian() +
        set_rescor(FALSE), data = dd)
}
models$importance <- function() {
  set.seed(7007)
  n <- 900
  dd <- data.frame(g = factor(rep(1:45, 20)))
  dd$x <- rnorm(nrow(dd))
  re <- rnorm(45, 0, 0.6)
  dd$y <- rpois(nrow(dd), exp(0.3 + 0.5 * dd$x + re[dd$g]))
  frm(bf(y ~ x + (1 | g)), family = poisson(), data = dd,
      control = frmtmb_control(importance = TRUE, importance_draws = 64,
                               importance_seed = 11L))
}
models$gp <- function() {
  set.seed(7008)
  n <- 250
  dd <- data.frame(t = sort(runif(n, 0, 10)))
  dd$y <- rnorm(n, sin(dd$t) + 0.3 * dd$t, 0.3)
  frm(bf(y ~ gp(t)), family = gaussian(), data = dd)
}
models$smooth <- function() {
  set.seed(7009)
  n <- 400
  dd <- data.frame(x = sort(runif(n, 0, 6)))
  dd$y <- rnorm(n, sin(dd$x), 0.4)
  frm(bf(y ~ s(x)), family = gaussian(), data = dd)
}
models$reml <- function() {
  set.seed(7010)
  ng <- 30
  dd <- data.frame(g = factor(rep(seq_len(ng), 25)))
  dd$x <- rnorm(nrow(dd))
  re <- rnorm(ng, 0, 0.8)
  dd$y <- rnorm(nrow(dd), 1 + 0.6 * dd$x + re[dd$g], 1)
  frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd, REML = TRUE)
}
models$zi <- function() {
  set.seed(7011)
  n <- 1500
  dd <- data.frame(x = rnorm(n))
  lam <- exp(0.5 + 0.4 * dd$x)
  dd$y <- ifelse(rbinom(n, 1, 0.3) == 1, 0L, rpois(n, lam))
  frm(bf(y ~ x, zi ~ 1), family = zero_inflated_poisson(), data = dd)
}
models$separation <- function() {
  set.seed(7012)
  n <- 400
  dd <- data.frame(x = rnorm(n))
  dd$y <- as.integer(dd$x > 0)             # perfectly separated
  frm(bf(y ~ x), family = bernoulli(), data = dd)
}
models$rare_top <- function() {
  set.seed(7013)
  n <- 4000
  dd <- data.frame(x = rnorm(n))
  eta <- 0.8 * dd$x
  dd$yo <- cut(eta + rlogis(n), breaks = c(-Inf, 0, 2, 9), labels = FALSE)
  cat("rare_top cell counts:", paste(table(dd$yo), collapse = " "), "\n")
  frm(bf(yo ~ x), family = cumulative(), data = dd)
}
models$student_nu <- function() {
  set.seed(7014)
  n <- 6000
  dd <- data.frame(x = rnorm(n))
  dd$y <- rnorm(n, 1 + 2 * dd$x, 1)        # nu runs to its upper edge
  frm(bf(y ~ x), family = student(), data = dd)
}
models$glmm_big <- function() {
  set.seed(7015)
  ng <- 2000
  dd <- data.frame(g = factor(rep(seq_len(ng), 5)))
  dd$x <- rnorm(nrow(dd))
  re <- rnorm(ng, 0, 0.7)
  dd$y <- rpois(nrow(dd), exp(0.3 + 0.5 * dd$x + re[dd$g]))
  frm(bf(y ~ x + (1 | g)), family = poisson(), data = dd)
}

out <- list()
for (nm in names(models)) {
  t0 <- proc.time()[["elapsed"]]
  r <- gw(models[[nm]]())
  el <- proc.time()[["elapsed"]] - t0
  out[[nm]] <- row_of(nm, r)
  out[[nm]]$elapsed <- el
  o <- out[[nm]]
  cat(sprintf("%-18s %7.2fs  err=%s conv=%s logLik=%s gmax=%s nwarn=%d\n",
              nm, el, o$err %||% NA,
              format(o$convergence %||% NA),
              format(o$logLik %||% NA, digits = 15),
              format(o$gmax %||% NA, digits = 6),
              length(o$warns %||% character(0))))
  for (w in o$warns %||% character(0)) cat("     WARN:", w, "\n")
}
saveRDS(out, OUT)
cat("DONE ident\n")
