# Short check of lane wcheck2's fix round: check_ald_laplace().
#   Rscript dev/reviews/wcheck2-short-ald.R <base|lane> <part>
# part zial:    zero_inflated_asym_laplace (1 | g) with zero-dominated
#               groups: frm()'s verdict, TMB's inner Hessian against
#               1 / sd^2, TMB's Laplace objective profiled along log sd
#               (rest re-optimized), and the exact-mode Laplace by hand
#               (inner max by optimize(), rest re-optimized by nlminb()).
# part weights: asym_laplace (1 | g) with weights in (0.2, 3), the same.
# part student: asym_laplace with (1 | g) beside a Student-t block
#               (1 | gr(h, dist = "student")): verdict and TMB profile
#               along g's log sd.
# part excl:    the check's exclusions, spot checks.
a <- commandArgs(TRUE)
lib_arg <- a[1]
part <- a[2]
source("dev/rel083-kt/wcheck2/wcheck2-common.R")
ns <- asNamespace("frmtmb")
grid <- c(-2, -1, 0, 1, 2, 4, 8)

get_obj <- function(fo, fam, d, ...) {
  ob <- q(frm(fo, family = fam, data = d, dry_run = "objective", ...))
  if (!is.null(ob$obj)) ob$obj else ob
}
safe_fn <- function(obj) function(v) {
  r <- tryCatch(obj$fn(v), error = function(e) NA_real_)
  if (!is.finite(r)) 1e10 else r
}
# TMB's profile along theta[ith] (log sd), every other outer parameter
# re-optimized from the previous step, from the optimum outward
tmb_profile <- function(obj, ith) {
  f0 <- safe_fn(obj)
  p <- obj$par
  op <- nlminb(p, f0, control = list(iter.max = 500, eval.max = 800))
  p <- op$par
  it <- which(names(p) == "theta")[ith]
  base <- p[it]
  free <- setdiff(seq_along(p), it)
  run <- function(gs) {
    pp <- p
    vapply(gs, function(t) {
      pp[it] <- base + t
      f <- function(r) { q2 <- pp; q2[free] <- r; f0(q2) }
      o <- nlminb(pp[free], f, control = list(iter.max = 500, eval.max = 800))
      if (o$objective < 1e9) pp[free] <<- o$par
      -o$objective
    }, 0)
  }
  up <- run(grid[grid >= 0])
  down <- rev(run(rev(grid[grid <= 0])))
  v <- c(down[-length(down)], up)
  list(v = v, base = base, opt = -op$objective, par = p)
}
show <- function(label, v) {
  cat(sprintf("%s\n  log sd + %s\n  value     %s\n  never falls: %s (largest fall %.3g)\n",
              label, paste(grid, collapse = " "),
              paste(format(v, digits = 8), collapse = " "),
              all(diff(v) >= -1e-6 * pmax(1, abs(v[-1]))), max(0, -diff(v))))
}
# TMB's inner Hessian at the inner mode for the outer parameters p,
# against the prior's precision 1 / sd^2 of a scalar (1 | g) block
hess_check <- function(obj, p, ith) {
  obj$fn(p)
  H <- as.matrix(obj$env$spHess(obj$env$last.par, random = TRUE))
  sd <- exp(p[which(names(p) == "theta")[ith]])
  cat(sprintf("  inner Hessian: %d x %d, max |H - I / sd^2| = %.3g (1 / sd^2 = %.6g)\n",
              nrow(H), ncol(H), max(abs(H - diag(1 / sd^2, nrow(H)))),
              1 / sd^2))
}
# exact-mode Laplace of y ~ x + (1 | g), ALD (and zi), at log sd ls:
# sum_g max_b [sum_i w_i log f(y_i | mu_i + b) - b^2 / (2 sd^2)], the
# prior's constant and -1/2 log H cancelling; rest maximized by nlminb()
ald_lpdf <- function(y, mu, sigma, p) {
  u <- (y - mu) / sigma
  log(p) + log(1 - p) - log(sigma) - 0.5 * (abs(u) + (2 * p - 1) * u)
}
exact_laplace <- function(d, w, ls, zi_fam, start) {
  sd <- exp(ls)
  nz <- if (zi_fam) d$y != 0 else rep(TRUE, nrow(d))
  gi <- split(seq_len(nrow(d)), d$g)
  fn <- function(th) {
    sigma <- exp(th[3]); p <- plogis(th[4])
    mu <- th[1] + th[2] * d$x
    tot <- 0
    for (rows in gi) {
      rr <- rows[nz[rows]]
      if (!length(rr)) next
      h <- function(b) sum(w[rr] * ald_lpdf(d$y[rr], mu[rr] + b, sigma, p)) -
        b^2 / (2 * sd^2)
      tot <- tot + optimize(h, c(-50, 50), maximum = TRUE, tol = 1e-10)$objective
    }
    if (zi_fam) {
      k0 <- sum(w[!nz]); k1 <- sum(w[nz])
      zi <- k0 / (k0 + k1)
      tot <- tot + k0 * log(zi) + k1 * log(1 - zi)
    }
    -tot
  }
  o <- nlminb(start, fn)
  list(v = -o$objective, par = o$par)
}
exact_profile <- function(d, w, zi_fam, ls0) {
  st <- c(median(d$y[d$y != 0]), 0, log(sd(d$y[d$y != 0]) / 2), 0)
  vapply(ls0 + grid, function(ls) {
    r <- exact_laplace(d, w, ls, zi_fam, st)
    st <<- r$par
    r$v
  }, 0)
}

if (part == "zial") {
  for (cfg in list(c(seed = 20261113, all0 = 10, per = 3),
                   c(seed = 20261114, all0 = 4, per = 2),
                   c(seed = 20261115, all0 = 0, per = 3))) {
    set.seed(cfg[["seed"]])
    G <- 30; n <- 4 * G
    d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(G), each = 4)))
    u <- rnorm(G, 0, 0.8)[as.integer(d$g)]
    d$y <- 1 + 0.5 * d$x + u + rnorm(n)
    # groups 1..all0 all zero; in the rest, `per` of 4 rows zero
    zero <- as.integer(d$g) <= cfg[["all0"]] |
      ave(seq_len(n), d$g, FUN = seq_along) <= cfg[["per"]]
    d$y[zero] <- 0
    cat(sprintf("--- zial seed %d: %d of %d rows zero, %d all-zero groups\n",
                cfg[["seed"]], sum(d$y == 0), n, cfg[["all0"]]))
    fo <- bf(y ~ x + (1 | g))
    fam <- zero_inflated_asym_laplace()
    cat("  frm():", outcome(frm(fo, family = fam, data = d)), "\n")
    obj <- get_obj(fo, fam, d)
    tp <- tmb_profile(obj, 1L)
    hess_check(obj, tp$par, 1L)
    show(sprintf("  TMB Laplace, optimum at log sd %.4g", tp$base), tp$v)
    show("  exact-mode Laplace by hand, same log sd grid",
         exact_profile(d, rep(1, n), TRUE, tp$base))
  }
}

if (part == "weights") {
  for (seed in c(20261113, 20261114)) {
    set.seed(seed)
    G <- 20; n <- 4 * G
    d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(G), each = 4)))
    d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 0.7)[as.integer(d$g)] + rnorm(n)
    d$w <- runif(n, 0.2, 3)
    cat(sprintf("--- weights seed %d, w in (%.3g, %.3g)\n", seed, min(d$w),
                max(d$w)))
    fo <- bf(y | weights(w) ~ x + (1 | g))
    cat("  frm():", outcome(frm(fo, family = asym_laplace(), data = d)), "\n")
    obj <- get_obj(fo, asym_laplace(), d)
    tp <- tmb_profile(obj, 1L)
    hess_check(obj, tp$par, 1L)
    show(sprintf("  TMB Laplace, optimum at log sd %.4g", tp$base), tp$v)
    show("  exact-mode Laplace by hand", exact_profile(d, d$w, FALSE, tp$base))
  }
}

if (part == "student") {
  for (cfg in list(list(nu = 2.5, h = "same"), list(nu = 5, h = "same"),
                   list(nu = 2.5, h = "crossed"))) {
    set.seed(20261116)
    G <- 20; n <- 8 * G
    d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(G), each = 8)))
    d$h <- if (cfg$h == "same") factor(paste0("h", d$g)) else
      factor(rep(seq_len(8), G))
    d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 1.5)[as.integer(d$g)] +
      rt(nlevels(d$h), 2.5)[as.integer(d$h)] + rnorm(n, 0, 0.5)
    fo <- bf(as.formula(sprintf(
      "y ~ x + (1 | g) + (1 | gr(h, dist = 'student', dist_nu = %g))", cfg$nu)))
    cat(sprintf("--- student nu %g, h %s\n", cfg$nu, cfg$h))
    fr <- q(frm(fo, family = asym_laplace(), data = d, dry_run = "frame"))
    cat("  check_ald_laplace():", tryCatch({
      ns$check_ald_laplace(fr$spec, fr, NULL); "NULL (open)"
    }, error = function(e) paste("REFUSED:", substr(conditionMessage(e), 1, 90))),
    "\n")
    ig <- NA_integer_
    for (bk in fr$re_blocks) {
      cat(sprintf("  block %s covstruct %s student %s theta_idx %s\n",
                  bk$term_label, bk$covstruct, !is.null(bk$dist_nu),
                  paste(bk$theta_idx, collapse = ",")))
      if (is.null(bk$dist_nu)) ig <- bk$theta_idx[1]
    }
    obj <- get_obj(fo, asym_laplace(), d)
    tp <- tmb_profile(obj, ig)
    show(sprintf("  TMB Laplace along g's log sd, optimum at %.4g", tp$base),
         tp$v)
  }
}

if (part == "excl") {
  set.seed(20261114)
  n <- 80
  d <- data.frame(x = rnorm(n), g = factor(rep(1:20, each = 4)),
                  tt = rep(1:4, 20))
  d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.7)[as.integer(d$g)] + rnorm(n)
  d$sdx <- 0.2
  AL <- asym_laplace()
  chk <- function(lab, fo, ...) {
    fr <- q(frm(fo, family = AL, data = d, dry_run = "frame", ...))
    r <- tryCatch({ ns$check_ald_laplace(fr$spec, fr, list(...)$prior,
                                         exact = isTRUE(list(...)$quadrature))
                    "NULL (open)" },
                  error = function(e) paste("REFUSED", substr(conditionMessage(e), 1, 60)))
    cat(sprintf("%-38s check: %s\n%-38s frm(): %s\n", lab, r, "",
                outcome(frm(fo, family = AL, data = d, ...))))
  }
  chk("trunc(lb = -20) + (1 | g)", bf(y | trunc(lb = -20) ~ x + (1 | g)))
  chk("me(x, sdx) + (1 | g)", bf(y ~ me(x, sdx) + (1 | g)))
  chk("ar(tt, gr = g) + (1 | g)", bf(y ~ x + ar(time = tt, gr = g) + (1 | g)))
  chk("quadrature = TRUE, (1 | g)", bf(y ~ x + (1 | g)), quadrature = TRUE)
  chk("nonlinear body + (1 | g)",
      bf(y ~ a + b * x, a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
      prior = c(set_prior("normal(0, 5)", nlpar = "a"),
                set_prior("normal(0, 5)", nlpar = "b")))
  chk("sd prior + (1 | g)", bf(y ~ x + (1 | g)),
      prior = set_prior("exponential(1)", class = "sd"))
  chk("(1 | g) control: no exclusion", bf(y ~ x + (1 | g)))
}
cat("DONE\n")
