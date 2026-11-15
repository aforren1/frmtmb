# Every item of lane refits, before and after. ARM=before reads rellib-r21.
# Run from the worktree root: Rscript dev/refits-repro.R
source("dev/refits-lib.R")
sec <- Sys.getenv("SEC", "ABCDEFGHIJKL")
want <- function(s) grepl(s, sec, fixed = TRUE)

## A. 1259: character ordinal response with non-numeric labels
if (want("A")) {
  set.seed(2501)
  n <- 60
  d <- data.frame(x = rnorm(n))
  lab <- c("none", "mild", "moderate", "severe")
  idx <- pmin(4, pmax(1, round(2.5 + d$x + rnorm(n, 0, .7))))
  d$y <- lab[idx]
  d$yn <- as.character(idx)
  d$yo <- factor(d$y, levels = lab, ordered = TRUE)
  go("A char labels", class(frm(bf(y ~ x), family = cumulative(), data = d)))
  go("A char digits", {
    f <- frm(bf(yn ~ x), family = cumulative(), data = d); fixef(f)[, 1] })
  go("A ordered factor", {
    f <- frm(bf(yo ~ x), family = cumulative(), data = d); fixef(f)[, 1] })
}

## B. 1285: fitted(mv, newdata) on a multivariate ordinal fit
if (want("B")) {
  set.seed(77)
  n <- 300
  d2 <- data.frame(x = rnorm(n), z = rnorm(n))
  d2$fc <- factor(sample(c("a", "b", "c"), n, TRUE))
  eff <- c(a = -1, b = 0, c = 1.5)[as.character(d2$fc)]
  p1 <- plogis(-0.3 + eff); p2 <- (1 - p1) * plogis(0.5 - eff)
  u <- runif(n); u2 <- runif(n)
  d2$yo <- ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L))
  d2$yo2 <- ifelse(u2 < p1, 1L, ifelse(u2 < p1 + p2, 2L, 3L))
  mv <- suppressWarnings(
    frm(bf(yo ~ x + cs(fc)) + bf(yo2 ~ z + cs(fc)), family = sratio(),
        data = d2))
  mv2 <- suppressWarnings(
    frm(bf(yo ~ x) + bf(yo2 ~ z), family = sratio(), data = d2))
  nd <- data.frame(x = 0, z = 0, fc = factor("a", levels = c("a", "b", "c")))
  nd2 <- data.frame(x = c(0, 1), z = c(0, -1),
                    fc = factor(c("a", "c"), levels = c("a", "b", "c")))
  go("B cs 1 row", dim(fitted(mv, newdata = nd)))
  go("B no cs 1 row", dim(fitted(mv2, newdata = nd)))
  go("B no cs 2 rows", dim(fitted(mv2, newdata = nd2)))
  go("B 1 row equals row 1 of 2 rows", {
    a <- fitted(mv2, newdata = nd)
    b <- fitted(mv2, newdata = nd2)
    c(identical(dimnames(a)[[3]], dimnames(b)[[3]]),
      max(abs(a[1, , ] - b[1, , ])))
  })
  go("B one response 1 row", dim(fitted(mv2, newdata = nd, resp = "yo")))
}

## C. 1567: y ~ . - w refit on data without w
if (want("C")) {
  set.seed(3)
  d4 <- data.frame(x = rnorm(80), z = rnorm(80), w = rnorm(80))
  d4$y <- d4$x + rnorm(80, 0, exp(0.3 * d4$w))
  go("C y ~ . - w refit without w", {
    f <- frm(y ~ . - w, data = d4)
    print(formula(f))
    d5 <- d4; d5$w <- NULL
    f2 <- update(f, newdata = d5)
    c(identical(fixef(f2), fixef(update(f, newdata = d4))),
      rownames(fixef(f2)))
  })
  go("C y ~ x + w - w refit without w", {
    f <- frm(y ~ x + w - w, data = d4)
    d5 <- d4; d5$w <- NULL
    rownames(fixef(update(f, newdata = d5)))
  })
}

## D. 1687: refit(newresp =) of a recoded bernoulli fit
if (want("D")) {
  set.seed(9)
  d7 <- data.frame(x = rnorm(100))
  d7$yb <- ifelse(runif(100) < plogis(d7$x), -1, -2)
  fb <- frm(yb ~ x, family = bernoulli(), data = d7)
  go("D refit original values", {
    r <- refit(fb, newresp = d7$yb)
    c(coef_equal = max(abs(fixef(r)[, 1] - fixef(fb)[, 1])))
  })
  go("D refit simulated codes", {
    s <- simulate(fb, nsim = 1, seed = 1)[[1]]
    r <- refit(fb, newresp = s); fixef(r)[, 1] })
  go("D refit factor original", {
    d7$yf <- factor(ifelse(d7$yb == -1, "yes", "no"))
    ff <- frm(yf ~ x, family = bernoulli(), data = d7)
    r <- refit(ff, newresp = d7$yf)
    max(abs(fixef(r)[, 1] - fixef(ff)[, 1]))
  })
}

## E. 2939: (1 + f || gr(g, cov = A))
if (want("E")) {
  set.seed(11)
  ng <- 12
  A <- 0.5 ^ abs(outer(1:ng, 1:ng, "-"))
  dimnames(A) <- list(paste0("g", 1:ng), paste0("g", 1:ng))
  d8 <- data.frame(g = factor(rep(paste0("g", 1:ng), each = 10),
                              levels = paste0("g", 1:ng)),
                   f = factor(rep(c("a", "b"), 60)))
  d8$y <- rnorm(120)
  go("E (0 + f || gr(g, cov = A))", {
    f <- frm(y ~ 1 + (0 + f || gr(g, cov = A)), data = d8,
             data2 = list(A = A))
    as.numeric(logLik(f))
  })
}

## F. 3225: group-level cs()
if (want("F")) {
  set.seed(12)
  d9 <- data.frame(g = factor(rep(1:15, each = 20)), x = rnorm(300))
  d9$y <- sample(1:4, 300, TRUE)
  go("F (cs(1) | g) sratio", {
    f <- frm(y ~ x + (cs(1) | g), family = sratio(), data = d9)
    as.numeric(logLik(f))
  })
}

## G. 3430: vcov_cluster outside the model; update() evaluates at the caller
if (want("G")) {
  set.seed(1)
  G <- 20
  dd <- data.frame(g = factor(rep(seq_len(G), each = 6)),
                   cl = factor(rep(seq_len(G / 2), each = 12)),
                   x = rnorm(G * 6))
  dd$y <- 1 + 0.5 * dd$x + rnorm(G * 6)
  fit <- frm(y ~ x, data = dd)
  v1 <- go("G cluster outside, data present",
           sqrt(diag(vcov_cluster(fit, ~ cl))))
  dd_saved <- dd
  dd$cl <- factor(rep(seq_len(G / 4), each = 24))
  go("G cluster outside, data changed",
     sqrt(diag(vcov_cluster(fit, ~ cl))))
  rm(dd)
  go("G cluster outside, data gone", sqrt(diag(vcov_cluster(fit, ~ cl))))
  # update(): na.action of the stored call, made inside a function
  mk <- function(d) {
    my_na <- stats::na.exclude
    frm(y ~ x, data = d, na.action = my_na)
  }
  dna <- dd_saved; dna$x[3] <- NA
  fn <- mk(dna)
  go("G update() of a fit whose na.action was a local", {
    u <- update(fn)
    c(length(residuals(u)), length(residuals(fn)))
  })
  go("G update(newdata) of a fit whose na.action was a local", {
    u <- update(fn, newdata = dna)
    c(length(residuals(u)), length(residuals(fn)))
  })
  mk2 <- function(d) {
    s <- TRUE
    frm(y ~ x, data = d, se = s)
  }
  fs <- mk2(dd_saved)
  go("G update() of a fit whose se was a local", class(update(fs)))
  s <- "not a flag"
  go("G update() with a global of the same name", class(update(fs)))
  rm(s)
}

## H. 4375: confint(parm = "x") on a saved ar(cov = TRUE) fit
if (want("H")) {
  set.seed(4)
  dh <- data.frame(t = rep(1:12, 5), h = factor(rep(1:5, each = 12)),
                   x = rnorm(60))
  dh$y <- 1 + 0.5 * dh$x + rnorm(60)
  fh <- frm(y ~ x + ar(time = t, gr = h, cov = TRUE), data = dh)
  old <- fh
  old$frame$autocor[[1L]]$sigma_is <- NULL
  go("H confint(parm = 'x') old", confint(old, parm = "x"))
  go("H confint(parm = 'x') new", confint(fh, parm = "x"))
  go("H confint() old", confint(old))
}

## I. 3188: marginaleffects predictions() on a matrix response
if (want("I")) {
  set.seed(5)
  dm <- data.frame(x = rnorm(200), N = 10L)
  pr <- cbind(1, exp(0.3 + 0.5 * dm$x), exp(-0.2 * dm$x))
  pr <- pr / rowSums(pr)
  dm$Y <- t(sapply(1:200, function(i) rmultinom(1, 10, pr[i, ])))
  colnames(dm$Y) <- c("a", "b", "c")
  fm <- frm(Y | trials(N) ~ x, family = multinomial(), data = dm)
  go("I predictions()", {
    p <- marginaleffects::predictions(fm)
    print(head(as.data.frame(p)[, c("group", "estimate", "std.error")]))
    dim(p)
  })
  go("I avg_predictions()", {
    p <- marginaleffects::avg_predictions(fm)
    as.data.frame(p)[, c("group", "estimate", "std.error")]
  })
  go("I fitted() mean", {
    ft <- fitted(fm)
    colMeans(ft[, 1, ])
  })
}

## J. 4018: recover_data() keeps the class of the error
if (want("J")) {
  set.seed(6)
  dj <- data.frame(x = rnorm(50), f = factor(rep(c("a", "b"), 25)))
  dj$y <- 1 + dj$x + rnorm(50)
  fj <- frm(y ~ x + f, data = dj)
  go("J emmeans dpar unknown", {
    r <- tryCatch(emmeans::emmeans(fj, ~ f, dpar = "nope"),
                  error = function(e) e)
    list(class = class(r), msg = conditionMessage(r))
  })
}

## K. 3427: getCall() of a delta update
if (want("K")) {
  set.seed(7)
  dk <- data.frame(x = rnorm(40), x3 = rnorm(40), g = factor(rep(1:8, 5)))
  dk$z <- dk$x + rnorm(8)[dk$g] + rnorm(40)
  fk <- frm(z ~ x + x3 + (1 | g), data = dk)
  go("K getCall(update(fit, ~ . - x3))", {
    u <- update(fk, ~ . - x3)
    cl <- getCall(u)
    print(cl)
    c(nchar = nchar(paste(deparse(cl), collapse = "")),
      same = identical(fixef(update(u)), fixef(u)))
  })
  go("K getCall(update(fit, . ~ . + x3)) of the reduced", {
    u <- update(fk, ~ . - x3)
    u2 <- update(u, . ~ . + x3)
    print(getCall(u2))
    max(abs(fixef(u2)[, 1] - fixef(fk)[, 1]))
  })
  go("K delta on a sigma model", {
    fs <- frm(bf(z ~ x + x3, sigma ~ x), data = dk)
    u <- update(fs, ~ . - x3)
    print(getCall(u))
    rownames(fixef(u))
  })
}

## L. 3434: influence() with s(x): the deletion refits keep the basis
if (want("L")) {
  set.seed(8)
  dl <- data.frame(x = runif(60, 0, 3), z = rnorm(60))
  dl$y <- sin(dl$x) + 0.3 * dl$z + rnorm(60, 0, 0.3)
  fl <- frm(y ~ z + s(x, k = 5), data = dl)
  go("L influence obs 1..60: max |dfbetas| per column", {
    inf <- influence(fl)
    # a deletion moves a coefficient by about its se / sqrt(n); a moved
    # basis moves the null-space column by far more
    round(apply(abs(dfbetas(inf)), 2, max), 3)
  })
  go("L deletion of row 1 against the fit's basis written out", {
    inf <- influence(fl)
    # the fit's own basis at the 59 rows: the wiggly columns as one iid
    # block with one sd (homdiag), the null-space column as a fixed one;
    # that is the model a deletion refit that keeps the basis fits
    si <- fl$frame$linpreds[[1]]$smooths[[1]]
    M <- frmtmb:::smooth_basis_at(si, dl[-1, ])
    nr <- sum(si$nr)
    dd <- dl[-1, ]
    for (k in seq_len(nr)) dd[[paste0("r", k)]] <- M[, k]
    dd$fx <- M[, nr + 1L]
    dd$id <- factor(1)
    fo <- as.formula(paste0("y ~ z + fx + homdiag(0 + ",
                            paste0("r", seq_len(nr), collapse = " + "),
                            " | id)"))
    ref <- frm(fo, data = dd)
    rb <- fixef(ref)[, 1]
    c(table_fx = unname(inf$fixed[1, "s(x).fx1"]),
      written_fx = unname(rb["fx"]),
      table_z = unname(inf$fixed[1, "z"]), written_z = unname(rb["z"]),
      full_fx = unname(inf$fixed_full["s(x).fx1"]))
  })
}
cat("\nDONE\n")
