# Shapes beyond the lane's grid, each against the same model on the
# centered column. REV_ARM selects the build. Seeds in each block.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
only <- Sys.getenv("REV_ONLY", "")

fx <- function(f, nm) {
  fe <- tryCatch(quiet(fixef(f)), error = function(e) NULL)
  if (is.null(fe) || !nm %in% rownames(fe)) return(c(NA, NA))
  unname(fe[nm, 1:2])
}

cmp <- function(lbl, aw, az, nw, nz, extra = NULL) {
  if (nzchar(only) && !grepl(only, lbl)) return(invisible())
  t0 <- proc.time()[[3]]
  rw <- do.call(rfit, aw)
  rz <- do.call(rfit, az)
  if (rw$err || rz$err) {
    cat(sprintf("%-46s ERROR w: %s | z: %s\n", lbl,
                if (rw$err) conditionMessage(rw$f) else "ok",
                if (rz$err) conditionMessage(rz$f) else "ok"))
    return(invisible())
  }
  ew <- fx(rw$f, nw); ez <- fx(rz$f, nz)
  dll <- ll_of(rw$f) - ll_of(rz$f)
  cat(sprintf(paste0("%-46s basis=%s conv=%s dll=%+.3e shift=%+.3e SE ",
                     "seratio=%.5f t=%.1f\n    warn w: %s\n    msg w: %s\n"),
              lbl, basis_on(rw$f), rw$f$opt$convergence, dll,
              (ew[1] - ez[1]) / ez[2], ew[2] / ez[2],
              proc.time()[[3]] - t0, short_w(rw),
              if (length(rw$msg)) substr(paste(unique(rw$msg), collapse = " | "), 1, 200) else "none"))
  if (length(rz$w)) cat("    warn z:", short_w(rz), "\n")
  if (!is.null(extra)) extra(rw$f, rz$f)
  invisible(list(w = rw$f, z = rz$f))
}

mk <- function(seed, n = 150, G = 15, m = 1e6) {
  set.seed(seed)
  g <- factor(rep(seq_len(G), length.out = n))
  z <- rnorm(n); z2 <- rnorm(n)
  f <- factor(sample(letters[1:3], n, TRUE))
  o <- sample(1:4, n, TRUE)
  t <- runif(n); wt <- runif(n, 0.5, 2); expo <- runif(n, 1, 3)
  u <- rnorm(G, 0, 0.7)[g]
  eta <- 0.5 + 0.4 * z + u + 0.2 * (o - 2)
  data.frame(g = g, z = z, w = z + m, z2 = z2, f = f, o = o, t = t,
             wt = wt, expo = expo, u = u,
             yp = rpois(n, exp(eta + log(expo))),
             yb = rbinom(n, 1, plogis(eta)),
             yg = eta + rnorm(n, 0, 0.5),
             yo = cut(eta + rlogis(n), c(-Inf, 0, 0.8, 1.6, Inf),
                      ordered_result = TRUE, labels = FALSE))
}

for (s in 1:2) {
  d <- mk(s)
  P <- poisson(); B <- bernoulli(); Gs <- gaussian()
  cat("== data seed", s, "\n")
  cmp(paste("mo(): yp ~ w + mo(o) + (1|g)", s),
      list(bf(yp ~ w + mo(o) + (1 | g)), family = P, data = d),
      list(bf(yp ~ z + mo(o) + (1 | g)), family = P, data = d), "w", "z")
  cmp(paste("offset: yp ~ w + offset(log(expo)) + (1|g)", s),
      list(bf(yp ~ w + offset(log(expo)) + (1 | g)), family = P, data = d),
      list(bf(yp ~ z + offset(log(expo)) + (1 | g)), family = P, data = d),
      "w", "z")
  cmp(paste("weights: yp|weights(wt) ~ w + (1|g)", s),
      list(bf(yp | weights(wt) ~ w + (1 | g)), family = P, data = d),
      list(bf(yp | weights(wt) ~ z + (1 | g)), family = P, data = d),
      "w", "z")
  cmp(paste("s(): yp ~ w + s(t) + (1|g)", s),
      list(bf(yp ~ w + s(t) + (1 | g)), family = P, data = d),
      list(bf(yp ~ z + s(t) + (1 | g)), family = P, data = d), "w", "z")
  cmp(paste("factor: yp ~ w + f + (1|g)", s),
      list(bf(yp ~ w + f + (1 | g)), family = P, data = d),
      list(bf(yp ~ z + f + (1 | g)), family = P, data = d), "w", "z")
  cmp(paste("random slope: yp ~ w + (1 + w|g)", s),
      list(bf(yp ~ w + (1 + w | g)), family = P, data = d),
      list(bf(yp ~ z + (1 + z | g)), family = P, data = d), "w", "z")
  cmp(paste("cs(): acat yo ~ cs(w) + (1|g)", s),
      list(bf(yo ~ cs(w) + (1 | g)), family = acat(), data = d),
      list(bf(yo ~ cs(z) + (1 | g)), family = acat(), data = d),
      "w[1]", "z[1]")
  cmp(paste("cumulative yo ~ w + (1|g)", s),
      list(bf(yo ~ w + (1 | g)), family = cumulative(), data = d),
      list(bf(yo ~ z + (1 | g)), family = cumulative(), data = d), "w", "z")
  cmp(paste("nlpar: yg ~ a + b*t, a ~ 1 + w + (1|g)", s),
      list(bf(yg ~ a + b * t, a ~ 1 + w + (1 | g), b ~ 1, nl = TRUE),
           family = Gs, data = d,
           prior = prior(normal(0, 5), nlpar = "b")),
      list(bf(yg ~ a + b * t, a ~ 1 + z + (1 | g), b ~ 1, nl = TRUE),
           family = Gs, data = d,
           prior = prior(normal(0, 5), nlpar = "b")),
      "a_w", "a_z")
  cmp(paste("nl body: yg ~ b0 + b1 * w (1 + w scale), nl", s),
      list(bf(yg ~ b0 + b1 * w, b0 ~ 1 + (1 | g), b1 ~ 1, nl = TRUE),
           family = Gs, data = d),
      list(bf(yg ~ b0 + b1 * z, b0 ~ 1 + (1 | g), b1 ~ 1, nl = TRUE),
           family = Gs, data = d), "b1_Intercept", "b1_Intercept")
  cmp(paste("zi: zero_inflated_poisson zi ~ w", s),
      list(bf(yp ~ w + (1 | g), zi ~ w), family = zero_inflated_poisson(),
           data = d),
      list(bf(yp ~ z + (1 | g), zi ~ z), family = zero_inflated_poisson(),
           data = d), "w", "z")
  cmp(paste("bernoulli quadrature: yb ~ w + (1|g)", s),
      list(bf(yb ~ w + (1 | g)), family = B, data = d, quadrature = TRUE),
      list(bf(yb ~ z + (1 | g)), family = B, data = d, quadrature = TRUE),
      "w", "z")
  cmp(paste("REML gaussian yg ~ w + (1|g)", s),
      list(bf(yg ~ w + (1 | g)), family = Gs, data = d, REML = TRUE),
      list(bf(yg ~ z + (1 | g)), family = Gs, data = d, REML = TRUE),
      "w", "z")
  cmp(paste("optimizer optim: yp ~ w + (1|g)", s),
      list(bf(yp ~ w + (1 | g)), family = P, data = d,
           control = frmtmb_control(optimizer = "optim")),
      list(bf(yp ~ z + (1 | g)), family = P, data = d,
           control = frmtmb_control(optimizer = "optim")), "w", "z")
}
