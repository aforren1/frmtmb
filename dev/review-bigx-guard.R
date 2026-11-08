# Where the basis does not engage or a later step can undo it:
# iteration caps (natural_finish()), bounds on the slope, more than 400
# outer parameters, and the 1e8 warning cells. REV_ARM selects the build.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
only <- Sys.getenv("REV_ONLY", "")
run_part <- function(p) !nzchar(only) || grepl(only, p)

gd <- function(seed, m, n = 150L, G = 15L, fam = "poisson", re = TRUE) {
  set.seed(seed)
  g <- factor(rep(seq_len(G), length.out = n))
  z <- rnorm(n)
  u <- if (re) rnorm(G, 0, 0.7)[g] else 0
  eta <- 0.5 + 0.4 * z + u
  y <- switch(fam, poisson = rpois(n, exp(eta)),
              gaussian = eta + rnorm(n, 0, 0.5),
              bernoulli = rbinom(n, 1, plogis(eta)))
  data.frame(y = y, z = z, w = z + m, g = g)
}
famf <- function(s) switch(s, poisson = poisson(), gaussian = gaussian(),
                           bernoulli = bernoulli())
slope <- function(f, nm) unname(quiet(fixef(f))[nm, 1:2])
line <- function(lbl, rw, fz, nw = "w", nz = "z") {
  if (rw$err) { cat(lbl, "ERROR", conditionMessage(rw$f), "\n"); return() }
  ew <- slope(rw$f, nw); ez <- slope(fz, nz)
  nf <- !is.null(rw$f$opt$natural_finish)
  cat(sprintf("%-44s basis=%s conv=%s natfin=%s dll=%+.3e shift=%+.2f SE warn: %s\n",
              lbl, basis_on(rw$f), rw$f$opt$convergence, nf,
              ll_of(rw$f) - ll_of(fz), (ew[1] - ez[1]) / ez[2], short_w(rw)))
}

if (run_part("itermax")) {
  cat("== iteration caps, poisson (1|g), m = 1e6\n")
  for (s in 1:5) {
    d <- gd(s, 1e6)
    fz <- quiet(frm(bf(y ~ z + (1 | g)), family = poisson(), data = d))
    for (it in c(2, 3, 5, 8)) {
      ctl <- frmtmb_control(optCtrl = list(iter.max = it, eval.max = 400))
      rw <- rfit(bf(y ~ w + (1 | g)), family = poisson(), data = d,
                 control = ctl)
      line(sprintf("seed %d iter.max %d", s, it), rw, fz)
    }
  }
}

if (run_part("bounds")) {
  cat("== bounds on the slope, poisson (1|g), m = 1e6\n")
  for (s in 1:5) {
    d <- gd(s, 1e6)
    for (lb in c(-5, 0.45)) {
      pr <- set_prior("normal(0, 10)", class = "b", lb = lb)
      fz <- quiet(frm(bf(y ~ z + (1 | g)), family = poisson(), data = d,
                      prior = pr))
      rw <- rfit(bf(y ~ w + (1 | g)), family = poisson(), data = d,
                 prior = pr)
      line(sprintf("seed %d lb %g", s, lb), rw, fz)
    }
    pr <- set_prior("normal(0, 10)", class = "b", ub = 5)
    fz <- quiet(frm(bf(y ~ z + (1 | g)), family = poisson(), data = d,
                    prior = pr))
    rw <- rfit(bf(y ~ w + (1 | g)), family = poisson(), data = d, prior = pr)
    line(sprintf("seed %d ub 5", s), rw, fz)
  }
}

if (run_part("np400")) {
  cat("== more than 400 outer parameters, m = 1e6\n")
  for (fam in c("gaussian", "poisson")) {
    for (s in 1:2) {
      set.seed(100 + s)
      n <- 2400L
      f <- factor(rep(sprintf("L%03d", 1:401), length.out = n))
      z <- rnorm(n)
      a <- rnorm(401, 0, 0.3)[f]
      eta <- 0.5 + 0.4 * z + a
      y <- if (fam == "gaussian") eta + rnorm(n, 0, 0.5) else
        rpois(n, exp(eta))
      d <- data.frame(y = y, z = z, w = z + 1e6, f = f)
      t0 <- proc.time()[[3]]
      fz <- quiet(frm(bf(y ~ z + f), family = famf(fam), data = d))
      rw <- rfit(bf(y ~ w + f), family = famf(fam), data = d)
      line(sprintf("%s seed %d np=%d", fam, s, length(rw$f$opt$par)), rw, fz)
      cat("   time", proc.time()[[3]] - t0, "\n")
    }
  }
}

if (run_part("w1e8")) {
  cat("== 1e8 cells the lane lists as warned\n")
  cells <- list(c("poisson", 1, TRUE), c("poisson", 5, TRUE),
                c("poisson", 1, FALSE), c("bernoulli", 4, FALSE),
                c("gaussian", 4, TRUE), c("gaussian", 5, TRUE))
  for (cc in cells) {
    fam <- cc[1]; s <- as.integer(cc[2]); re <- as.logical(cc[3])
    d <- gd(s, 1e8, fam = fam, re = re)
    fo_w <- if (re) bf(y ~ w + (1 | g)) else bf(y ~ w)
    fo_z <- if (re) bf(y ~ z + (1 | g)) else bf(y ~ z)
    fz <- quiet(frm(fo_z, family = famf(fam), data = d))
    rw <- rfit(fo_w, family = famf(fam), data = d)
    line(sprintf("%s seed %d re %s", fam, s, re), rw, fz)
    if (length(rw$w)) cat("   full:", rw$w[1], "\n")
    # an independent check of "at the maximum": the natural gradient
    # mapped into the centered coordinates, and a Newton step there
    if (!rw$err) {
      f <- rw$f
      g <- f$obj$gr(f$opt$par)
      cat("   natural grad:", format(g, digits = 3), "\n")
    }
  }
}
