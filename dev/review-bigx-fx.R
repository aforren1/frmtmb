# Short check of the fix round. REV_ARM selects the build, REV_ONLY a part.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
only <- Sys.getenv("REV_ONLY", "")
run_part <- function(p) !nzchar(only) || grepl(only, p)
wcls <- function(r) {
  if (!length(r$w)) return("none")
  paste(paste0("[", r$wc, "] ", substr(r$w, 1, 150)), collapse = " || ")
}
fe <- function(f, nm) {
  x <- tryCatch(quiet(fixef(f)), error = function(e) NULL)
  if (is.null(x) || !nm %in% rownames(x)) c(NA, NA) else unname(x[nm, 1:2])
}
cmpl <- function(lbl, rw, fz, nw, nz) {
  if (rw$err) { cat(lbl, "ERROR", conditionMessage(rw$f), "\n"); return() }
  a <- fe(rw$f, nw); b <- fe(fz, nz)
  cat(sprintf("%-48s conv=%s dll=%+.3e shift=%+.2e SE se_ratio=%.5f\n   warn(%d): %s\n",
              lbl, rw$f$opt$convergence, ll_of(rw$f) - ll_of(fz),
              (a[1] - b[1]) / b[2], a[2] / b[2], length(rw$w), wcls(rw)))
  if (length(rw$msg)) cat("   msg:", substr(paste(unique(rw$msg), collapse = " | "), 1, 200), "\n")
}
gd <- function(seed, m, n = 150L, G = 15L) {
  set.seed(seed)
  g <- factor(rep(seq_len(G), length.out = n))
  z <- rnorm(n); z2 <- rnorm(n)
  f <- factor(sample(c("a", "b", "c"), n, TRUE))
  u <- rnorm(G, 0, 0.7)[g]
  eta <- 0.5 + 0.4 * z - 0.3 * z2 + 0.3 * (f == "b") + u
  data.frame(y = rpois(n, exp(eta)), z = z, w = z + m, z2 = z2,
             w2 = z2 - m / 2, f = f, g = g,
             yg = 1 + 0.4 * z + u + rnorm(n, 0, 0.5))
}
P <- poisson()

if (run_part("bounds")) {
  cat("== B1 bounds, poisson (1|g), 1e6\n")
  for (s in 1:3) {
    d <- gd(s, 1e6)
    pr2 <- set_prior("normal(0, 10)", class = "b", lb = -5)
    fz <- quiet(frm(bf(y ~ z + z2 + (1 | g)), family = P, data = d, prior = pr2))
    rw <- rfit(bf(y ~ w + w2 + (1 | g)), family = P, data = d, prior = pr2)
    cmpl(sprintf("two shifted slopes, lb=-5 on b, seed %d", s), rw, fz, "w", "z")
    if (!rw$err) cat("   w2 shift/SE:", (fe(rw$f, "w2")[1] - fe(fz, "z2")[1]) / fe(fz, "z2")[2], "\n")
    pr3 <- set_prior("normal(0, 10)", class = "b", lb = 0.35)
    fz <- quiet(frm(bf(y ~ z + z2 + (1 | g)), family = P, data = d, prior = pr3))
    rw <- rfit(bf(y ~ w + w2 + (1 | g)), family = P, data = d, prior = pr3)
    cmpl(sprintf("two shifted slopes, lb=0.35 (binds w2), seed %d", s), rw, fz, "w", "z")
    # a bound on a factor level's coefficient only
    prf <- tryCatch(set_prior("normal(0, 5)", class = "b", coef = "fb", lb = -1),
                    error = function(e) e)
    if (inherits(prf, "error")) cat("coef-level bound refused:", conditionMessage(prf), "\n") else {
      fz <- rfit(bf(y ~ z + f + (1 | g)), family = P, data = d, prior = prf)
      rw <- rfit(bf(y ~ w + f + (1 | g)), family = P, data = d, prior = prf)
      if (fz$err) cat("factor-level bound, z fit ERROR", conditionMessage(fz$f), "\n") else
        cmpl(sprintf("lb on factor level fb only, seed %d", s), rw, fz$f, "w", "z")
    }
    # a bound on the intercept
    pri <- set_prior("normal(0, 10)", class = "Intercept", lb = -50)
    fz <- quiet(frm(bf(y ~ z + (1 | g)), family = P, data = d, prior = pri))
    rw <- rfit(bf(y ~ w + (1 | g)), family = P, data = d, prior = pri)
    cmpl(sprintf("lb on Intercept, seed %d", s), rw, fz, "w", "z")
    if (!rw$err) cat("   ill_design warnings:", sum(grepl("frmtmb_ill_design", rw$wc)), "\n")
    # a bound on the intercept of a two-slope design
    rw <- rfit(bf(y ~ w + w2 + (1 | g)), family = P, data = d, prior = pri)
    cat("   two slopes + lb Intercept: ill_design warnings:",
        if (rw$err) "ERROR" else sum(grepl("frmtmb_ill_design", rw$wc)), "|",
        if (!rw$err) wcls(rw), "\n")
  }
}

if (run_part("nl")) {
  cat("== M1 nl bodies, 1e6\n")
  for (s in 1:3) {
    d <- gd(s, 1e6)
    rw <- rfit(bf(yg ~ b0 + b1 * w, b0 ~ 1 + (1 | g), b1 ~ 1, nl = TRUE),
               family = gaussian(), data = d)
    fz <- quiet(frm(bf(yg ~ b0 + b1 * z, b0 ~ 1 + (1 | g), b1 ~ 1, nl = TRUE),
                    family = gaussian(), data = d))
    cmpl(sprintf("b0 + b1*w, seed %d", s), rw, fz, "b1_Intercept", "b1_Intercept")
    # a curved body in a shifted covariate: exact reparametrization
    # b0' = b0 * exp(b1 * m); poisson mean b0 * exp(b1 * w)
    for (m in c(10, 1e2, 1e3)) {
      d$wm <- d$z + m
      rw <- rfit(bf(y ~ b0 * exp(b1 * wm), b0 ~ 1, b1 ~ 1, nl = TRUE),
                 family = poisson(link = "identity"), data = d,
                 prior = c(set_prior("", nlpar = "b0", lb = 0)))
      rz <- rfit(bf(y ~ b0 * exp(b1 * z), b0 ~ 1, b1 ~ 1, nl = TRUE),
                 family = poisson(link = "identity"), data = d,
                 prior = c(set_prior("", nlpar = "b0", lb = 0)))
      if (rz$err) cat("centered curved body ERROR", conditionMessage(rz$f), "\n") else
        cmpl(sprintf("b0*exp(b1*(z+%g)) identity poisson, seed %d", m, s), rw,
             rz$f, "b1_Intercept", "b1_Intercept")
    }
    # a warning from the kept point must reach the user: b2 * b3 is
    # identified only as a product, so the fit's checks name a flat ridge
    d$t <- runif(nrow(d))
    pr23 <- c(set_prior("normal(0.5, 10)", nlpar = "b2"),
              set_prior("normal(0.5, 10)", nlpar = "b3"))
    rw <- rfit(bf(yg ~ b0 + b1 * w + b2 * b3 * t, b0 ~ 1 + (1 | g), b1 ~ 1,
                  b2 ~ 1, b3 ~ 1, nl = TRUE), family = gaussian(), data = d,
               prior = pr23)
    rz <- rfit(bf(yg ~ b0 + b1 * z + b2 * b3 * t, b0 ~ 1 + (1 | g), b1 ~ 1,
                  b2 ~ 1, b3 ~ 1, nl = TRUE), family = gaussian(), data = d,
               prior = pr23)
    if (!rw$err && !rz$err) {
      cmpl(sprintf("b0 + b1*w + b2*b3*t (ridge), seed %d", s), rw, rz$f,
           "b1_Intercept", "b1_Intercept")
      cat("   centered fit's warnings(", length(rz$w), "): ", wcls(rz), "\n", sep = "")
    }
  }
  cat("== mo(), 1e6\n")
  for (s in 1:2) {
    set.seed(s)
    n <- 150; g <- factor(rep(1:15, length.out = n)); z <- rnorm(n)
    o <- sample(1:4, n, TRUE)
    eta <- 0.5 + 0.4 * z + 0.2 * (o - 2) + rnorm(15, 0, 0.7)[g]
    d <- data.frame(y = rpois(n, exp(eta)), z = z, w = z + 1e6, o = o, g = g)
    fz <- quiet(frm(bf(y ~ z + mo(o) + (1 | g)), family = P, data = d))
    rw <- rfit(bf(y ~ w + mo(o) + (1 | g)), family = P, data = d)
    cmpl(sprintf("w + mo(o) + (1|g), seed %d", s), rw, fz, "w", "z")
  }
}

if (run_part("rank")) {
  cat("== M2 rank rule\n")
  # column order: x3 = x1 + x2 + tiny noise; any two span the model
  for (s in 1:2) {
    set.seed(600 + s); n <- 200
    x1 <- rnorm(n); x2 <- rnorm(n)
    d <- data.frame(x1 = x1, x2 = x2, x3 = x1 + x2 + 1e-9 * rnorm(n),
                    junk = 1e-3 * rnorm(n))
    d$y <- 1 + 2 * x1 - x2 + rnorm(n)
    d$xd <- d$junk + 1e-12 * rnorm(n) + 0 * x1
    for (fo in list(y ~ x1 + x2 + x3, y ~ x3 + x1 + x2, y ~ x3 + x2 + x1)) {
      r <- rfit(bf(fo), family = gaussian(), data = d)
      g <- glm(fo, data = d)
      cat(sprintf("%-22s kept %s  ll %.6f glm %.6f\n   msg: %s\n   msgclass/warn: %s\n",
                  deparse(fo), paste(names(beta_of(r$f)), collapse = ","),
                  ll_of(r$f), as.numeric(logLik(g)),
                  paste(unique(r$msg), collapse = " | "), wcls(r)))
    }
    # the important column after a near-alias of the intercept:
    # junk2 = 1 + 1e-9 noise (aliased with the intercept), then x1
    d$k <- 1 + 1e-9 * rnorm(n)
    for (fo in list(y ~ k + x1 + x2, y ~ x1 + x2 + k)) {
      r <- rfit(bf(fo), family = gaussian(), data = d)
      g <- glm(fo, data = d)
      cat(sprintf("%-22s kept %s  ll %.6f glm %.6f\n   msg: %s\n",
                  deparse(fo), paste(names(beta_of(r$f)), collapse = ","),
                  ll_of(r$f), as.numeric(logLik(g)),
                  paste(unique(r$msg), collapse = " | ")))
    }
  }
  # message class
  set.seed(1); n <- 100
  d <- data.frame(x = rnorm(n)); d$k <- 5; d$y <- rnorm(n)
  cls <- NULL
  withCallingHandlers(quiet0 <- frm(bf(y ~ x + k), data = d),
                      message = function(m) { cls <<- class(m); invokeRestart("muffleMessage") })
  cat("drop message class:", paste(cls, collapse = "/"), "\n")
  # a rare factor level next to a shifted column
  for (s in 1:2) {
    set.seed(700 + s); n <- 2000
    f <- factor(c("rare", sample(c("a", "b", "c"), n - 1, TRUE)))
    z <- rnorm(n)
    eta <- 0.3 + 0.4 * z + 0.3 * (f == "b") + 0.5 * (f == "rare")
    d <- data.frame(y = rpois(n, exp(eta)), z = z, w = z + 1e6, f = f)
    fz <- quiet(frm(bf(y ~ z + f), family = P, data = d))
    rw <- rfit(bf(y ~ w + f), family = P, data = d)
    cmpl(sprintf("w + f (1 rare row of 2000), seed %d", 700 + s), rw, fz, "w", "z")
    cat("   kept:", paste(names(beta_of(rw$f)), collapse = ","), "\n")
    a <- fe(rw$f, "frare"); b <- fe(fz, "frare")
    cat("   frare shift/SE:", (a[1] - b[1]) / b[2], "\n")
  }
  # interactions and powers of shifted columns, against the centered model
  for (m in c(1e2, 1e3, 1e4, 1e5)) {
    set.seed(42); n <- 300
    z1 <- rnorm(n); z2 <- rnorm(n)
    d <- data.frame(z1 = z1, z2 = z2, w1 = z1 + m, w2 = z2 + m)
    d$y <- 0.5 + 0.4 * z1 - 0.3 * z2 + 0.3 * z1 * z2 + rnorm(n, 0, 0.5)
    d$yp <- rpois(n, exp(0.2 + 0.3 * z1 + 0.2 * z1^2 / 4))
    fz <- quiet(frm(bf(y ~ z1 * z2), data = d))
    rw <- rfit(bf(y ~ w1 * w2), data = d)
    g <- glm(y ~ w1 * w2, data = d)
    cat(sprintf("gaussian w1*w2, m=%g: kept %s dll(vs centered) %+.3e glm dll %+.3e\n   inter shift/SE %+.2e  msg: %s warn: %s\n", m,
                paste(names(beta_of(rw$f)), collapse = ","), ll_of(rw$f) - ll_of(fz),
                as.numeric(logLik(g)) - ll_of(fz),
                (fe(rw$f, "w1:w2")[1] - fe(fz, "z1:z2")[1]) / fe(fz, "z1:z2")[2],
                paste(unique(rw$msg), collapse = " | "), wcls(rw)))
    fz <- quiet(frm(bf(yp ~ z1 + I(z1^2)), family = P, data = d))
    rw <- rfit(bf(yp ~ w1 + I(w1^2)), family = P, data = d)
    cat(sprintf("poisson w1 + I(w1^2), m=%g: kept %s dll %+.3e  quad shift/SE %+.2e msg: %s warn: %s\n", m,
                paste(names(beta_of(rw$f)), collapse = ","), ll_of(rw$f) - ll_of(fz),
                (fe(rw$f, "Iw1E2")[1] - fe(fz, "Iz1E2")[1]) / fe(fz, "Iz1E2")[2],
                paste(unique(rw$msg), collapse = " | "), wcls(rw)))
  }
  # the drop's limit on a plain shift
  for (m in c(1e9, 3e9, 1e10)) {
    set.seed(3); n <- 150; z <- rnorm(n)
    d <- data.frame(z = z, w = z + m, y = 0.5 + 0.4 * z + rnorm(n, 0, 0.5))
    fz <- quiet(frm(bf(y ~ z), data = d))
    rw <- rfit(bf(y ~ w), data = d)
    cat(sprintf("gaussian y ~ w, m=%g: kept %s dll %+.3e msg: %s warn: %s\n", m,
                paste(names(beta_of(rw$f)), collapse = ","), ll_of(rw$f) - ll_of(fz),
                paste(unique(rw$msg), collapse = " | "), wcls(rw)))
  }
}

if (run_part("profile")) {
  cat("== M3 profile and uniroot\n")
  for (s in 1:3) {
    d <- gd(s, 1e6)
    fz <- quiet(frm(bf(y ~ z + (1 | g)), family = P, data = d))
    fw <- quiet(frm(bf(y ~ w + (1 | g)), family = P, data = d))
    V0 <- vcov(fw); ll0 <- ll_of(fw); lp0 <- fw$obj$env$last.par
    r <- rfit_ci <- withCallingHandlers(
      tryCatch(list(p = confint(fw, parm = "w", method = "profile"),
                    u = confint(fw, parm = "w", method = "uniroot")),
               error = function(e) e),
      warning = function(x) { cat("   warning:", conditionMessage(x), "\n"); invokeRestart("muffleWarning") })
    pz <- quiet(confint(fz, parm = "z", method = "profile"))
    uz <- quiet(confint(fz, parm = "z", method = "uniroot"))
    if (inherits(r, "error")) cat("seed", s, "ERROR", conditionMessage(r), "\n") else {
      cat(sprintf("seed %d profile w %s | z %s\n        uniroot w %s | z %s\n", s,
                  paste(format(unlist(r$p[1, 1:2]), digits = 7), collapse = " "),
                  paste(format(unlist(pz[1, 1:2]), digits = 7), collapse = " "),
                  paste(format(unlist(r$u[1, 1:2]), digits = 7), collapse = " "),
                  paste(format(unlist(uz[1, 1:2]), digits = 7), collapse = " ")))
    }
    V1 <- vcov(fw)
    cat("   vcov identical after:", identical(V0, V1), " logLik identical:",
        identical(ll0, ll_of(fw)), " last.par identical:",
        identical(lp0, fw$obj$env$last.par), "\n")
    # predictions read the tape's state
    cat("   fitted after == fitted centered (rel):",
        max(abs(fitted(fw)[, 1] / fitted(fz)[, 1] - 1)), "\n")
  }
}
