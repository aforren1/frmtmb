# Review of lane wcheck2, item 4977: the pair rule against an independent
# linear program. For a dispersion design X (sigma, log link), weights w
# and the set of columns whose prior is free, the Laplace objective of
# y | weights(w) ~ x + (1 | obs) along a ray d in sigma's coefficients
# grows at the rate
#   r(d) = sum_{s_i < 0} (w_i - 1) |s_i| - sum_{s_i > 0} w_i s_i,  s = X d
# (a shrinking row gains its excess, a growing row loses its weight),
# with the centered intercept held (colMeans(X)'d = 0) under a bounding
# Intercept prior, and every column with a bounding prior held at 0.
# Unbounded when max r > 0. The LP is solved two ways: quadprog with a
# 1e-7 ridge (variables d, a = s^-, c = s^+, sum(a + c) <= 1), and a
# random search over unit directions. The check's verdict comes from
# check_weights_unbounded() on the model's frame.
#   Rscript dev/reviews/wcheck2-review-pair.R <lane|lib> <part>
# part: grid (the verdict against the LP over designs), profile (the
# objective along chosen directions)
a <- commandArgs(TRUE)
lib_arg <- a[1]
part <- a[2]
source("dev/wcheck2-common.R")
suppressMessages(library(quadprog))
P <- function(...) set_prior(..., dpar = "sigma")
pic <- P("normal(0, 1)", class = "Intercept")
olre <- "gives every observation its own random effect"
ns <- asNamespace("frmtmb")

lp_max <- function(X, w, fixed = integer(), hold_ic = TRUE) {
  n <- nrow(X); p <- ncol(X)
  e <- ifelse(w > 0, w - 1, 0)
  l <- w
  nv <- p + 2 * n
  Aeq <- cbind(X, diag(n), -diag(n))           # X d + a - c = 0
  if (hold_ic) Aeq <- rbind(Aeq, c(colMeans(X), rep(0, 2 * n)))
  for (j in fixed) {
    r <- numeric(nv); r[j] <- 1; Aeq <- rbind(Aeq, r)
  }
  Ain <- rbind(cbind(matrix(0, 2 * n, p), diag(2 * n)),
               c(rep(0, p), rep(-1, 2 * n)))
  A <- rbind(Aeq, Ain)
  b <- c(rep(0, nrow(Aeq)), rep(0, 2 * n), -1)
  dv <- c(rep(0, p), e, -l)
  sol <- solve.QP(diag(1e-7, nv), dv, t(A), b, meq = nrow(Aeq))
  x <- sol$solution
  list(val = sum(dv * x), d = x[seq_len(p)])
}

# random search: d in the null space of the held constraints
rand_max <- function(X, w, fixed = integer(), hold_ic = TRUE, K = 20000) {
  p <- ncol(X)
  C <- NULL
  if (hold_ic) C <- rbind(C, colMeans(X))
  for (j in fixed) { r <- numeric(p); r[j] <- 1; C <- rbind(C, r) }
  N <- if (is.null(C)) diag(p) else MASS::Null(t(C))
  if (!ncol(N)) return(list(val = 0))
  e <- ifelse(w > 0, w - 1, 0)
  rate <- function(z) {
    d <- N %*% z
    s <- as.numeric(X %*% d)
    (sum(e * pmax(-s, 0)) - sum(w * pmax(s, 0))) / sum(abs(s))
  }
  Z <- matrix(rnorm(K * ncol(N)), K)
  v <- apply(Z, 1, rate)
  list(val = max(v, na.rm = TRUE))
}

verdict <- function(d, fo, prior) {
  fr <- allow_quiet(frm(fo, family = gaussian(), data = d, dry_run = "frame"))
  r <- tryCatch(ns$check_weights_unbounded(fr$spec, fr, prior),
                error = function(e) conditionMessage(e))
  if (is.null(r)) "pass" else if (grepl("no maximum", r)) "REFUSE" else
    paste("ERR", substr(r, 1, 80))
}
allow_quiet <- function(e) suppressWarnings(suppressMessages(e))

sigma_X <- function(d, fo) {
  fr <- allow_quiet(frm(fo, family = gaussian(), data = d, dry_run = "frame"))
  as.matrix(fr$linpreds[[grep("sigma", names(fr$linpreds))[1]]]$X)
}

if (part == "grid") {
  set.seed(1)
  designs <- list(
    three = list(lev = 3, sizes = c(30, 30, 30), fo = ~ 1 + fac),
    four_unbal = list(lev = 4, sizes = c(10, 20, 30, 40), fo = ~ 1 + fac),
    inter = list(lev = 0, sizes = NULL, fo = ~ 1 + f1 * f2),
    fac_x = list(lev = 3, sizes = c(30, 30, 30), fo = ~ 1 + fac + z),
    sum3 = list(lev = 3, sizes = c(30, 30, 30), fo = ~ 1 + fac,
                contr = TRUE),
    three_bnorm = list(lev = 3, sizes = c(30, 30, 30), fo = ~ 1 + fac,
                       pb = "facb"),
    three_w0 = list(lev = 3, sizes = c(30, 30, 30), fo = ~ 1 + fac,
                    w0 = TRUE))
  for (dn in names(designs)) {
    ds <- designs[[dn]]
    if (dn == "inter") {
      n <- 96
      d <- data.frame(f1 = factor(rep(c("A", "B"), each = 48)),
                      f2 = factor(rep(rep(c("C", "D"), each = 24), 2)))
    } else {
      n <- sum(ds$sizes)
      d <- data.frame(fac = factor(rep(letters[seq_len(ds$lev)], ds$sizes)))
    }
    d$obs <- factor(seq_len(n)); d$x <- rnorm(n); d$z <- rnorm(n)
    d$y <- 0.5 + 0.4 * d$x + rnorm(n)
    if (isTRUE(ds$contr)) contrasts(d$fac) <- contr.sum(3)
    fo <- bf(y | weights(w) ~ x + (1 | obs))
    fo <- bf(as.formula(paste("y | weights(w) ~ x + (1 | obs)")),
             as.formula(paste("sigma", deparse(ds$fo))))
    prior <- pic
    if (!is.null(ds$pb)) prior <- c(pic, P("normal(0, 1)", class = "b",
                                           coef = ds$pb))
    d$w <- 1
    X <- sigma_X(d, fo)
    fixed <- if (!is.null(ds$pb)) which(colnames(X) == ds$pb) else integer()
    tab <- c(agree_ref = 0, agree_pass = 0, false_refusal = 0, miss = 0,
             lp_disagree = 0)
    worst_miss <- NULL
    for (k in 1:150) {
      # cell weights in (0.2, 1.8), half the configurations mixed within
      # a cell (row weights spread +-0.3 around the cell's)
      if (dn == "inter") {
        cellid <- interaction(d$f1, d$f2)
      } else cellid <- d$fac
      cw <- runif(nlevels(cellid), 0.2, 1.8)
      w <- cw[as.integer(cellid)]
      if (k %% 2 == 0) w <- pmax(0.05, w + runif(n, -0.3, 0.3))
      if (isTRUE(ds$w0)) w[sample(n, 9)] <- 0
      d$w <- w
      v <- verdict(d, fo, prior)
      lm <- lp_max(X, w, fixed)
      rm <- rand_max(X, w, fixed, K = 4000)
      lp_pos <- lm$val > 1e-5
      if (xor(lp_pos, rm$val > 1e-6)) tab[["lp_disagree"]] <- tab[["lp_disagree"]] + 1
      if (v == "REFUSE" && lp_pos) tab[["agree_ref"]] <- tab[["agree_ref"]] + 1
      if (v == "pass" && !lp_pos) tab[["agree_pass"]] <- tab[["agree_pass"]] + 1
      if (v == "REFUSE" && !lp_pos) {
        tab[["false_refusal"]] <- tab[["false_refusal"]] + 1
        cat(sprintf("  FALSE REFUSAL %s k %d lp %.3g rand %.3g w %s\n", dn,
                    k, lm$val, rm$val, paste(format(cw, digits = 3),
                                             collapse = " ")))
      }
      if (v == "pass" && lp_pos) {
        tab[["miss"]] <- tab[["miss"]] + 1
        if (is.null(worst_miss) || lm$val > worst_miss$val) {
          worst_miss <- list(val = lm$val, k = k, cw = cw,
                             d = lm$d)
        }
      }
      if (startsWith(v, "ERR")) cat("  ", dn, k, v, "\n")
    }
    cat(sprintf("%-12s cols %s | %s\n", dn,
                paste(colnames(X), collapse = ","),
                paste(names(tab), tab, sep = " ", collapse = ", ")))
    if (!is.null(worst_miss)) {
      cat(sprintf("  largest miss: k %d LP rate %.4g, cell weights %s, d %s\n",
                  worst_miss$k, worst_miss$val,
                  paste(format(worst_miss$cw, digits = 3), collapse = " "),
                  paste(format(worst_miss$d, digits = 3), collapse = " ")))
    }
  }
}

# the objective along a direction d of sigma's coefficients (b0 first),
# the rest re-optimized, from the start
profile_dir <- function(ob, v, g) {
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p <- ob$par
  isd <- which(names(p) == "betad")[seq_along(v)]
  base <- p[isd]
  free <- setdiff(seq_along(p), isd)
  vapply(g, function(t) {
    p[isd] <- base + t * v
    fr <- function(r) {
      q2 <- p; q2[free] <- r
      val <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
      if (!is.finite(val)) 1e10 else val
    }
    op <- tryCatch(nlminb(p[free], fr, control = list(iter.max = 500,
                                                       eval.max = 800)),
                   error = function(e) NULL)
    if (is.null(op)) return(NA_real_)
    if (op$objective < 1e9) p[free] <<- op$par
    -op$objective
  }, 0)
}

if (part == "profile") {
  # (A) the interaction design, cell weights chosen so that the LP is
  # positive and no pair is: the check passes. (B) 4 unbalanced cells
  # near the pair boundary, refused. (C) fac + z with a free z: the LP
  # moves z. For each: the verdict, the LP rate and direction, the
  # objective's slope along it, and frm()'s outcome.
  run <- function(lab, d, fo, prior, fixed = integer()) {
    X <- sigma_X(d, fo)
    lm <- lp_max(X, d$w, fixed)
    v <- verdict(d, fo, prior)
    # scale d so that sum |X d| = n (one unit of log sigma per row on
    # average); the predicted slope is then n * LP value
    s <- as.numeric(X %*% lm$d)
    dd <- lm$d * nrow(X) / sum(abs(s))
    pred <- lm$val * nrow(X)
    ob <- q(frm(fo, family = gaussian(), data = d, prior = prior,
                dry_run = "objective"))
    g <- c(0, 2, 4, 8, 12, 16, 20)
    vals <- profile_dir(ob, dd, g)
    cat(sprintf("--- %s: verdict %s, LP rate %.4g, d %s\n", lab, v, lm$val,
                paste(format(dd, digits = 3), collapse = " ")))
    show_profile(sprintf("  profile (predicted slope %.4g)", pred), g, vals)
    cat("  frm():", outcome(frm(fo, family = gaussian(), data = d,
                                 prior = prior)), "\n")
  }
  set.seed(7)
  n <- 96
  d <- data.frame(f1 = factor(rep(c("A", "B"), each = 48)),
                  f2 = factor(rep(rep(c("C", "D"), each = 24), 2)))
  d$obs <- factor(seq_len(n)); d$x <- rnorm(n)
  d$y <- 0.5 + 0.4 * d$x + rnorm(n)
  fo <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + f1 * f2)
  cell <- interaction(d$f1, d$f2)
  # A.C, B.C, A.D, B.D
  for (cw in list(c(1.0, 1.0, 1.0, 1.25), c(0.4, 1.0, 1.0, 1.25),
                  c(1.25, 0.4, 1.0, 1.0), c(1.25, 1.0, 0.4, 1.0))) {
    d$w <- cw[as.integer(cell)]
    run(paste("inter, A.C B.C A.D B.D =", paste(cw, collapse = " ")), d, fo,
        pic)
  }
  set.seed(8)
  sz <- c(10, 20, 30, 40)
  n <- sum(sz)
  d <- data.frame(fac = factor(rep(letters[1:4], sz)))
  d$obs <- factor(seq_len(n)); d$x <- rnorm(n); d$z <- rnorm(n)
  d$y <- 0.5 + 0.4 * d$x + rnorm(n)
  fo <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac)
  # a falls (10 rows, E/n 0.12), d rises (40 rows, W/n 0.1): G +0.02
  for (cw in list(c(1.12, 1, 1, 0.1), c(1.08, 1, 1, 0.1))) {
    d$w <- cw[as.integer(d$fac)]
    run(paste("four unbalanced, a b c d =", paste(cw, collapse = " ")), d,
        fo, pic)
  }
  # fac + z with z free: rows high in z fall, low in z rise
  fo <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac + z)
  d$w <- ifelse(d$z > 0, 1.3, 0.9)
  run("four + z, w 1.3 where z > 0 else 0.9", d, fo, pic)
}

if (part == "miss") {
  # the grid's largest misses: the LP is positive and the check passes
  run <- function(lab, d, fo, prior, fixed = integer()) {
    X <- sigma_X(d, fo)
    lm <- lp_max(X, d$w, fixed)
    v <- verdict(d, fo, prior)
    s <- as.numeric(X %*% lm$d)
    dd <- lm$d * nrow(X) / sum(abs(s))
    pred <- lm$val * nrow(X)
    ob <- q(frm(fo, family = gaussian(), data = d, prior = prior,
                dry_run = "objective"))
    g <- c(0, 2, 4, 8, 12, 16, 20)
    vals <- profile_dir(ob, dd, g)
    cat(sprintf("--- %s: verdict %s, LP rate %.4g, d %s\n", lab, v, lm$val,
                paste(format(dd, digits = 3), collapse = " ")))
    show_profile(sprintf("  profile (predicted slope %.4g)", pred), g, vals)
    cat("  frm():", outcome(frm(fo, family = gaussian(), data = d,
                                 prior = prior)), "\n")
  }
  set.seed(7)
  n <- 96
  d <- data.frame(f1 = factor(rep(c("A", "B"), each = 48)),
                  f2 = factor(rep(rep(c("C", "D"), each = 24), 2)))
  d$obs <- factor(seq_len(n)); d$x <- rnorm(n)
  d$y <- 0.5 + 0.4 * d$x + rnorm(n)
  cell <- interaction(d$f1, d$f2)
  fo <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + f1 * f2)
  d$w <- c(0.254, 0.881, 1.779, 0.654)[as.integer(cell)]
  run("inter, A.C B.C A.D B.D = 0.254 0.881 1.779 0.654", d, fo, pic)
  set.seed(9)
  n <- 90
  d <- data.frame(fac = factor(rep(c("a", "b", "c"), each = 30)))
  d$obs <- factor(seq_len(n)); d$x <- rnorm(n)
  d$y <- 0.5 + 0.4 * d$x + rnorm(n)
  fo <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac)
  d$w <- c(0.226, 0.914, 1.783)[as.integer(d$fac)]
  ds <- d
  contrasts(ds$fac) <- contr.sum(3)
  run("sum contrasts, a b c = 0.226 0.914 1.783", ds, fo, pic)
  run("treatment contrasts, same weights", d, fo, pic)
  d$w <- c(0.524, 0.416, 1.774)[as.integer(d$fac)]
  pb <- c(pic, P("normal(0, 1)", class = "b", coef = "facb"))
  run("normal(0, 1) on facb, a b c = 0.524 0.416 1.774", d, fo, pb,
      fixed = 2L)
}
cat("DONE\n")
