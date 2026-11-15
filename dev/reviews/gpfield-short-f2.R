# Short check of lane gpfield's fix round, F2: the Newton route of
# gp_field_draw(). Per replicate: does gp_field_newton() settle (and in how
# many steps), and how far is its z from the tape's inner mode, in
# posterior sds of z (the inner Hessian's inverse at the replicate)?
# Also: does the factor cache leak into the fit, and is predict()
# reproducible across calls with other ndraws and newdata?
# Usage: Rscript gpfield-short-f2.R <case> <ndraws>
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
a <- commandArgs(TRUE)
cs <- a[1]; nd_ <- as.integer(a[2])
gf_libs("lane")
ns <- asNamespace("frmtmb")
set.seed(9950 + match(cs, c("pois300", "bern200", "poislow200", "bern60",
                             "binom80", "gauss200")))
mkx <- function(n) sort(round(runif(n, 0, 10), 3))
row <- switch(cs,
  pois300 = { x <- sort(round(runif(300, 0, 50), 4))
    list(d = data.frame(x = x, y = rpois(300, exp(1 + sin(x / 3)))),
         fam = poisson()) },
  bern200 = { x <- mkx(200)
    list(d = data.frame(x = x, y = rbinom(200, 1, plogis(2 * sin(3 * x)))),
         fam = bernoulli()) },
  poislow200 = { x <- mkx(200)
    list(d = data.frame(x = x, y = rpois(200, exp(-1 + 1.5 * sin(3 * x)))),
         fam = poisson()) },
  bern60 = { x <- mkx(60)
    list(d = data.frame(x = x, y = rbinom(60, 1, plogis(2 * sin(3 * x)))),
         fam = bernoulli()) },
  binom80 = { x <- mkx(80)
    list(d = data.frame(x = x, tr = 20,
                        y = rbinom(80, 20, plogis(3 * sin(x)))),
         fam = binomial()) },
  gauss200 = { x <- mkx(200)
    list(d = data.frame(x = x, y = sin(x) + rnorm(200, 0, 0.3)),
         fam = gaussian()) })
fml <- if (cs == "binom80") bf(y | trials(tr) ~ gp(x)) else bf(y ~ gp(x))
q <- gf_fit(fml, data = row$d, family = row$fam)
fit <- q$fit
cat(sprintf("F2FIT %s conv %d ll %.6f warn [%s]\n", cs, fit$opt$convergence,
            as.numeric(logLik(fit)),
            paste(unique(sub(":.*", "", q$warn)), collapse = ";")))
rec <- new.env(); rec$r <- list()
orig_draw <- ns$gp_field_draw
orig_newton <- ns$gp_field_newton
set_steps <- function(k) {
  unlockBinding("gp_field_newton_steps", ns)
  assign("gp_field_newton_steps", k, envir = ns)
  lockBinding("gp_field_newton_steps", ns)
}
gf_put("gp_field_draw", function(fs, plan, xn) {
  nt <- plan$newton
  steps <- NA_integer_; zn <- NULL
  if (!is.null(nt)) {
    for (k in 1:8) {
      set_steps(k)
      zn <- tryCatch(orig_newton(fs, nt), error = function(e) NULL)
      if (!is.null(zn)) { steps <- k; break }
    }
    set_steps(8L)
  }
  # the tape's mode and the posterior sds of z there
  obj <- plan$obj; env <- obj$env
  par <- ns$outer_from_template(fs$estimates, obj, fs$frame)
  sv <- ns$obj_state_save(obj)
  v <- tryCatch(suppressWarnings(obj$fn(par)), error = function(e) NA)
  g <- plan$blocks[[1L]]
  zt <- env$last.par[g$pos]
  H <- as.matrix(env$spHess(env$last.par, random = TRUE))
  psd <- sqrt(diag(solve(0.5 * (H + t(H)))))
  ns$obj_state_restore(obj, sv)
  err <- if (!is.null(zn)) max(abs(zn$z - zt) / psd) else NA
  rec$r[[length(rec$r) + 1L]] <- c(newton = !is.null(nt), steps = steps,
                                    err = err, finite = is.finite(v))
  orig_draw(fs, plan, xn)
})
newd <- data.frame(x = c(row$d$x[c(5, 50)], max(row$d$x) + 1, 99))
if (!is.null(row$d$tr)) newd$tr <- 20
size0 <- as.numeric(object.size(fit))
set.seed(3)
p1 <- predict(fit, newdata = newd, summary = FALSE, ndraws = nd_)
gf_put("gp_field_draw", orig_draw)
m <- do.call(rbind, rec$r)
st <- table(factor(m[, "steps"], levels = 1:8), useNA = "always")
cat(sprintf("F2NEWTON %s replicates %d | newton plan %s | settled %d, fell back to tape %d | steps %s | max err %.2e median %.2e (posterior sds of z)\n",
            cs, nrow(m), format(all(m[, "newton"] == 1)),
            sum(!is.na(m[, "steps"])), sum(is.na(m[, "steps"])),
            paste(names(st), st, sep = ":", collapse = " "),
            max(m[, "err"], na.rm = TRUE), stats::median(m[, "err"], na.rm = TRUE)))
# cache: does anything stay on the fit?
cached <- vapply(fit$frame$re_blocks, function(b) !is.null(b$gp_factor_cache), NA)
cat(sprintf("F2CACHE %s cache on fit blocks %s | object.size before %.0f after %.0f\n",
            cs, format(any(cached)), size0, as.numeric(object.size(fit))))
# reproducible across an intervening predict with other ndraws and newdata
set.seed(3)
p2 <- predict(fit, newdata = newd, summary = FALSE, ndraws = nd_)
set.seed(4)
pz <- predict(fit, newdata = data.frame(x = c(1, 2)), summary = FALSE,
              ndraws = 7)
set.seed(3)
p3 <- predict(fit, newdata = newd, summary = FALSE, ndraws = nd_)
cat(sprintf("F2REPRO %s first==second %s second==after-other %s | sd %s\n", cs,
            identical(p1, p2), identical(p2, p3),
            paste(round(apply(p1, 2, sd), 3), collapse = " ")))
