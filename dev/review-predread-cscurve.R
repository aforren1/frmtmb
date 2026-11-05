# Item 3: frm_curve() SEs and derivatives under layered families, against a Monte
# Carlo draw from the fit's coefficient covariance through the reader itself, finite
# differences of the fit, and brms's sampled draws on the same data. Seed 7951.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7951)
n <- 600
d <- data.frame(x = rnorm(n), z = rnorm(n), f = factor(sample(c("a", "b", "c"), n, TRUE)))
eta <- 0.8 * d$x
d$ord <- as.integer(cut(eta + ifelse(d$z > 0, 1.3, -1.3) * rlogis(n), c(-Inf, -0.5, 0.8, 2, Inf), labels = FALSE))
d$ord2 <- as.integer(cut(eta + 0.8 * as.integer(d$f) + rlogis(n), c(-Inf, 0.5, 2, 3.5, Inf), labels = FALSE))
comp <- rbinom(n, 1, 0.5)
d$cat <- as.integer(cut(0.7 * d$x + rlogis(n), c(-Inf, 0, 1.5, Inf), labels = FALSE))
d$mixc <- rpois(n, ifelse(comp == 1, exp(1.5 + 0.4 * d$x), exp(0.2 - 0.3 * d$x)))
grid <- data.frame(x = c(-1, 0.3, 1.2), z = c(0.7, -1, 0.4), f = factor(c("b", "c", "b"), levels = levels(d$f)))
mc_se <- function(fit, nd, ...) {
  # draws of the whole outer parameter vector, each read through frm_linpred()
  V <- as.matrix(vcov(fit, full = TRUE))
  est <- unlist(fit$estimates[intersect(names(fit$estimates), unique(sub("[0-9]+$", "", rownames(V))))])
  L <- chol((V + t(V)) / 2)
  pars <- fit$obj$par
  par_names <- names(pars)
  sims <- replicate(800, {
    p <- pars + as.numeric(crossprod(L, rnorm(length(pars))))
    f <- fit
    f$estimates <- utils::relist(unlist(fit$estimates), fit$estimates)
    idx <- split(seq_along(p), par_names)
    for (nm in names(idx)) f$estimates[[nm]][] <- p[idx[[nm]]]
    as.numeric(frm_linpred(f, newdata = nd, type = "link", ...))
  })
  apply(sims, 1, sd)
}
rel <- function(a, b) signif(max(abs(a - b) / abs(b)), 3)
for (sp in list(
  list(name = "sratio_cs", f = "bf(ord ~ x + cs(z))", fam = "sratio()", var = "z"),
  list(name = "acat_cs_f", f = "bf(ord2 ~ x + cs(f))", fam = "acat()", var = "x"),
  list(name = "cumulative_cs", f = "bf(ord ~ x + cs(z))", fam = "cumulative()", var = "z"),
  list(name = "categorical", f = "bf(cat ~ x)", fam = "categorical()", var = "x"),
  list(name = "mixture", f = "bf(mixc ~ x)", fam = "mixture(poisson(), poisson())", var = "x"))) {
  cat("\n===", sp$name, sp$f, sp$fam, "\n")
  fit <- quiet(frm(rev_ev(sp$f, "frmtmb"), family = rev_ev(sp$fam, "frmtmb"), data = d))
  print(round(unlist(fit$estimates[c("beta", "thres")]), 4))
  cu <- frm_curve(fit, newdata = grid, seed = 1)
  li <- frm_linpred(fit, newdata = grid, se.fit = TRUE)
  cat(".se vs se.fit rel:", rel(cu$.se, as.numeric(li$se.fit)), " check:", attr(cu, "check")$cov_rel_error, "\n")
  mc <- tryCatch(mc_se(fit, grid), error = function(e) e)
  if (inherits(mc, "error")) {
    cat("mc ERR", conditionMessage(mc), "\n")
  } else {
    cat(".se  ", signif(cu$.se, 4), "\nMC se", signif(mc, 4), "\nmax rel", rel(cu$.se, mc), " (MC error about", signif(1 / sqrt(2 * 800), 2), ")\n")
  }
  # the derivative against a central difference of the reader itself
  de <- frm_curve_deriv(fit, sp$var, newdata = grid, seed = 1)
  gp <- grid
  gm <- grid
  h <- 1e-5
  if (is.numeric(grid[[sp$var]])) {
    gp[[sp$var]] <- gp[[sp$var]] + h
    gm[[sp$var]] <- gm[[sp$var]] - h
    fd <- (as.numeric(frm_linpred(fit, newdata = gp, type = "link")) -
             as.numeric(frm_linpred(fit, newdata = gm, type = "link"))) / (2 * h)
    cat("deriv ", signif(de$.estimate, 6), "\nFD    ", signif(fd, 6), "\nmax abs diff", signif(max(abs(de$.estimate - fd)), 3), "\n")
    # the deriv's SE against a MC of the central difference
    dmc <- tryCatch({
      V <- as.matrix(vcov(fit, full = TRUE))
      L <- chol((V + t(V)) / 2)
      pars <- fit$obj$par
      idx <- split(seq_along(pars), names(pars))
      sims <- replicate(800, {
        p <- pars + as.numeric(crossprod(L, rnorm(length(pars))))
        f <- fit
        for (nm in names(idx)) f$estimates[[nm]][] <- p[idx[[nm]]]
        (as.numeric(frm_linpred(f, newdata = gp, type = "link")) -
            as.numeric(frm_linpred(f, newdata = gm, type = "link"))) / (2 * h)
      })
      apply(sims, 1, sd)
    }, error = function(e) e)
    if (!inherits(dmc, "error")) {
      cat("deriv .se", signif(de$.se, 4), "\nMC se    ", signif(dmc, 4), "\n")
    } else cat("deriv mc ERR", conditionMessage(dmc), "\n")
  }
  # brms's sampled draws, flat priors, the same data
  bfit <- tryCatch(quiet(brms::brm(rev_ev(sp$f, "brms"), family = rev_ev(sp$fam, "brms"), data = d,
                                   chains = 2, iter = 1500, refresh = 0, seed = 3, silent = 2, cores = 1,
                                   init = 0)), error = function(e) e)
  if (inherits(bfit, "error")) {
    cat("brms ERR", conditionMessage(bfit), "\n")
    next
  }
  pl <- posterior_linpred(bfit, newdata = grid, re_formula = NA)
  sv <- apply(pl, setdiff(seq_along(dim(pl)), 1L), sd)
  mv <- apply(pl, setdiff(seq_along(dim(pl)), 1L), mean)
  cat("brms mean", signif(as.numeric(mv), 4), "\nfrm_curve", signif(cu$.estimate, 4), "\n")
  cat("brms sd  ", signif(as.numeric(sv), 4), "\n.se      ", signif(cu$.se, 4), "\nratio .se / sd", signif(cu$.se / as.numeric(sv), 3), "\n")
  Sb <- cov(matrix(pl, nrow = dim(pl)[1L]))
  Sf <- attr(cu, "Sigma")
  cat("cov between rows: max rel of off-diagonals in units of the sds:",
      signif(max(abs(cov2cor(as.matrix(Sf)) - cov2cor(Sb))), 3), "\n")
  crit_brms <- {
    M <- matrix(pl, nrow = dim(pl)[1L])
    zz <- abs(sweep(M, 2, colMeans(M))) / matrix(apply(M, 2, sd), nrow(M), ncol(M), byrow = TRUE)
    quantile(apply(zz, 1, max), 0.95)
  }
  cat("simultaneous crit: frm_curve", signif(cu$.crit_sim[1], 4), " from brms draws", signif(crit_brms, 4), "\n")
}
