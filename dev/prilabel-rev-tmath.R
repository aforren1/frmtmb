# Review of lane prilabel: (a) per-coefficient sd prior math against a
# hand-computed penalty at fixed points; (b) the shared Student-t scale
# against my own closed-form multivariate t and Monte Carlo, on the
# TAPE at several points; brms's Stan code for the sharing rule; fits
# with three t terms on one factor. Data as test-student-share.R
# st_data() (seed 21); perturbation seed 31; MC seed 41.
.libPaths(c("C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
})
lg <- file.path(wt, "dev/prilabel-rev-log")
sink(file.path(lg, "tmath.txt"), split = TRUE)
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "; brms",
    format(packageVersion("brms")), "\n")
st_data <- function(n = 360, seed = 21) {
  set.seed(seed)
  g <- factor(sample(sprintf("g%02d", 1:24), n, TRUE))
  h <- factor(sample(sprintf("h%d", 1:8), n, TRUE))
  x <- rnorm(n)
  z <- rnorm(n)
  w <- sqrt(5 / rchisq(24, 5))
  u0 <- (w * rnorm(24, 0, 0.8))[as.integer(g)]
  u1 <- (w * rnorm(24, 0, 0.4))[as.integer(g)]
  y <- 1 + 0.5 * x + 0.3 * z + u0 + u1 * x + rnorm(8, 0, 0.5)[as.integer(h)] +
    rnorm(n, 0, 0.6)
  data.frame(y, x, z, g, h)
}
dd <- st_data()
dd$hg <- factor(ifelse(as.integer(dd$g) <= 12L, "u", "v"))
dd$yz <- dd$y + rnorm(nrow(dd))
try_ <- function(expr) tryCatch(expr, error = function(e) {
  paste("ERROR:", conditionMessage(e))
})
th_pos <- function(fit) which(names(fit$obj$par) == "theta")

cat("\n######## (a) per-coefficient sd prior: penalty at fixed points\n")
# pen(p) = fn_prior(p) - fn_flat(p); the hand value is minus the log
# density of exp(theta_k) plus theta_k (Jacobian of sd = exp(theta))
pen_check <- function(form, pr, targets, logdens, label) {
  f0 <- frm(form, data = dd)
  f1 <- frm(form, data = dd, prior = pr)
  pos <- th_pos(f1)
  set.seed(31)
  cat("\n--", label, "\n")
  for (r in 1:3) {
    p <- f1$obj$par + rnorm(length(f1$obj$par), 0, 0.3)
    pen <- as.numeric(f1$obj$fn(p) - f0$obj$fn(p))
    t <- p[pos[targets]]
    hand <- -sum(logdens(exp(t)) + t)
    cat(sprintf("  point %d: pen %.12f hand %.12f diff %.3e\n", r, pen,
                hand, pen - hand))
  }
  invisible(f1)
}
pr_n <- set_prior("normal(0, 0.5)", class = "sd", group = "g", coef = "x")
ln <- function(s) dnorm(s, 0, 0.5, log = TRUE)
# us block (1 + x | g): theta 1, 2 log sds, 3 cor; (1 | h) theta 4
pen_check(bf(y ~ x + (1 + x | g) + (1 | h)), pr_n, 2L, ln,
          "us (1 + x | g) + (1 | h), coef x, normal(0, 0.5)")
pen_check(bf(y ~ x + (1 + x | g) + (1 | h)),
          set_prior("exponential(2)", class = "sd", group = "g",
                    coef = "Intercept"), 1L,
          function(s) dexp(s, 2, log = TRUE),
          "us, coef Intercept, exponential(2)")
pen_check(bf(y ~ x + (1 + x | gr(g, dist = "student")) + (1 | h)), pr_n, 2L,
          ln, "Student-t us_t block, coef x (the scale)")
f_by <- pen_check(bf(y ~ x + (1 + x | gr(g, by = hg))), pr_n, c(2L, 5L), ln,
                  "gr(g, by = hg), coef x: both by-levels")
cat("  gr(by) theta layout:",
    paste(variables(f_by)[grepl("^sd_|^cor_", variables(f_by))],
          collapse = ", "), "\n")
# the class-wide sd prior on a one-sd model, for the constant
pen_check(bf(y ~ x + (1 | g)), set_prior("normal(0, 0.5)", class = "sd"),
          1L, ln, "class-wide sd prior on (1 | g), the convention")
pen_check(bf(y ~ x + (1 | g)),
          set_prior("normal(0, 0.5)", class = "sd", group = "g",
                    coef = "Intercept"), 1L, ln,
          "per-coef spelling of the same sd")

cat("\n######## (b1) brms's Stan code: which terms share dfm\n")
show_dfm <- function(f, ...) {
  sc <- brms::stancode(f, data = dd, ...)
  l <- strsplit(sc, "\n")[[1]]
  print(grep("dfm_|df_[0-9]", l, value = TRUE))
}
st <- function(v) paste0("gr(", v, ", dist = 'student')")
cat("\n-- mu and sigma terms on g\n")
show_dfm(brms::bf(y ~ x + (1 | gr(g, dist = "student")),
                  sigma ~ (1 | gr(g, dist = "student"))))
cat("\n-- three t terms on g: two in mu, one in sigma\n")
f3b <- brms::bf(y ~ x + (1 | gr(g, dist = "student")) +
                  (0 + x | gr(g, dist = "student")),
                sigma ~ (1 | gr(g, dist = "student")))
show_dfm(f3b)
cat("\n-- two responses, a t term on g in each\n")
show_dfm(brms::bf(y ~ x + (1 | gr(g, dist = "student"))) +
           brms::bf(yz ~ x + (1 | gr(g, dist = "student"))) +
           brms::set_rescor(FALSE))

cat("\n######## (b2) shared density on the tape, F1\n")
mvt_ld <- function(B, S, nu) {
  # my closed form: lgamma((nu+D)/2) - lgamma(nu/2) - D/2 log(nu pi)
  # - 1/2 log|S| - (nu+D)/2 log(1 + b'S^-1 b / nu), per column of B
  D <- nrow(B)
  Si <- solve(S)
  q <- colSums(B * (Si %*% B))
  sum(lgamma((nu + D) / 2) - lgamma(nu / 2) - D / 2 * log(nu * pi) -
        0.5 * as.numeric(determinant(S)$modulus) -
        (nu + D) / 2 * log1p(q / nu))
}
nu <- 4
F1 <- frm(bf(y ~ x + (1 + x | gr(g, dist = "student", dist_nu = 4)) +
               (0 + z | gr(g, dist = "student", dist_nu = 4))), data = dd)
bks <- F1$frame$re_blocks
cat("t_share:", deparse(lapply(bks, `[[`, "t_share")), "\n")
vc <- function(bk, th) {
  frmtmb:::covstruct_registry[[bk$covstruct]]$vcov(th[bk$theta_idx], bk)
}
lp_eta <- function(fit, pars, dpar) {
  lp <- Filter(function(l) identical(l$dpar, dpar), fit$frame$linpreds)
  if (!length(lp)) return(rep(pars$betad[1], nrow(dd)))
  lp <- lp[[1]]
  e <-as.vector(lp$X %*% pars[[lp$par]][lp$idx])
  if (!is.null(lp$Z)) e <- e + as.vector(lp$Z %*% pars$b)
  e
}
split_par <- function(fit, p) {
  out <- list()
  for (nm in unique(names(p))) out[[nm]] <- unname(p[names(p) == nm])
  out
}
ref_joint <- function(fit, p, Sfun, Bfun, nu) {
  pars <- split_par(fit, p)
  eta <- lp_eta(fit, pars, "mu")
  ls <- lp_eta(fit, pars, "sigma")
  -sum(dnorm(dd$y, eta, exp(ls), log = TRUE)) -
    mvt_ld(Bfun(pars), Sfun(pars), nu)
}
S_F1 <- function(pars) {
  S <- matrix(0, 3, 3)
  S[1:2, 1:2] <- vc(bks[[1]], pars$theta)
  S[3, 3] <- vc(bks[[2]], pars$theta)
  S
}
B_F1 <- function(pars) {
  rbind(matrix(pars$b[bks[[1]]$b_idx], 2), pars$b[bks[[2]]$b_idx])
}
set.seed(31)
p0 <- F1$obj$env$last.par.best
for (r in 1:4) {
  p <- p0 + if (r == 1) 0 else rnorm(length(p0), 0, 0.4)
  tape <- as.numeric(F1$obj$env$f(p))
  ref <- ref_joint(F1, p, S_F1, B_F1, nu)
  cat(sprintf("  point %d: tape %.10f ref %.10f rel %.3e\n", r, tape, ref,
              (tape - ref) / abs(ref)))
}
# Monte Carlo of the closed form at one level: mean over w of the
# gaussian density with covariance w^2 S, w^2 = nu / chi2(nu)
pars <- split_par(F1, p0)
S <- S_F1(pars)
bj <- B_F1(pars)[, 1]
set.seed(41)
w2 <- nu / rchisq(2e6, nu)
Si <- solve(S)
q <- drop(crossprod(bj, Si %*% bj))
ldS <- as.numeric(determinant(S)$modulus)
dens <- exp(-1.5 * log(2 * pi) - 1.5 * log(w2) - ldS / 2 - q / (2 * w2))
cat(sprintf("  MC level 1: log mean %.6f (se %.1e) closed form %.6f\n",
            log(mean(dens)), sd(dens) / sqrt(length(dens)) / mean(dens),
            mvt_ld(matrix(bj), S, nu)))

cat("\n######## (b3) three t terms on g, one in sigma (F3)\n")
F3 <- try_(frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 4)) +
                    (0 + x | gr(g, dist = "student", dist_nu = 4)),
                  sigma ~ (1 | gr(g, dist = "student", dist_nu = 4))),
               data = dd))
if (!is.character(F3)) {
  b3 <- F3$frame$re_blocks
  cat("t_share:", deparse(lapply(b3, `[[`, "t_share")), " conv",
      F3$opt$convergence, " logLik", format(as.numeric(logLik(F3)),
                                            digits = 12), "\n")
  S3 <- function(pars) {
    diag(vapply(b3, function(bk) as.numeric(vc(bk, pars$theta)), 0))
  }
  B3 <- function(pars) {
    do.call(rbind, lapply(b3, function(bk) pars$b[bk$b_idx]))
  }
  set.seed(31)
  q0 <- F3$obj$env$last.par.best
  for (r in 1:4) {
    p <- q0 + if (r == 1) 0 else rnorm(length(q0), 0, 0.4)
    tape <- as.numeric(F3$obj$env$f(p))
    ref <- ref_joint(F3, p, S3, B3, nu)
    cat(sprintf("  point %d: tape %.10f ref %.10f rel %.3e\n", r, tape, ref,
                (tape - ref) / abs(ref)))
  }
  # fitted against fixef + ranef by hand, mu and sigma
  re <- ranef(F3)
  cat("ranef groups:", names(re), "; dims", deparse(lapply(re, dim)), "\n")
  fx <- fitted(F3, dpar = "mu")[, 1]
  pars <- split_par(F3, q0)
  cat(sprintf("  fitted(mu) against the tape's eta: %.3e\n",
              max(abs(fx - lp_eta(F3, pars, "mu")))))
  fs <- try_(fitted(F3, dpar = "sigma")[, 1])
  if (!is.character(fs)) {
    cat(sprintf("  fitted(sigma) against the tape's exp(eta): %.3e\n",
                max(abs(fs - exp(lp_eta(F3, pars, "sigma"))))))
  }
  sm <- try_(simulate(F3, nsim = 3, seed = 2))
  cat("  simulate():", if (is.character(sm)) sm else dim(as.data.frame(sm)),
      "\n")
  sm2 <- try_(simulate(F3, nsim = 2, seed = 2, re_formula = NA))
  cat("  simulate(re_formula = NA):", if (is.character(sm2)) sm2 else "runs",
      "\n")
  sm3 <- try_(simulate(F3, nsim = 2, seed = 2,
                       re_formula = ~ (1 | gr(g, dist = "student"))))
  cat("  simulate(re_formula = one mu term):", substr(sm3, 1, 90), "\n")
  nd <- dd[1:3, ]
  nd$g <- factor("gNEW")
  pn <- try_(predict(F3, newdata = nd, allow_new_levels = TRUE)[, 1])
  cat("  predict new level:", if (is.character(pn)) pn else round(pn, 3),
      "\n")
  # the same model with sigma's term gaussian is a different model
  F3g <- frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 4)) +
                  (0 + x | gr(g, dist = "student", dist_nu = 4)),
                sigma ~ (1 | g)), data = dd)
  cat(sprintf("  logLik with sigma's term gaussian %.6f (t_share %s)\n",
              as.numeric(logLik(F3g)),
              deparse(lapply(F3g$frame$re_blocks, `[[`, "t_share"))))
}

cat("\n######## (b4) battery: shared terms against one diag block\n")
two <- frm(bf(y ~ x + (1 | gr(g, dist = "student")) +
                (0 + x | gr(g, dist = "student"))), data = dd)
one <- frm(bf(y ~ x + (1 + x || gr(g, dist = "student"))), data = dd)
cat(sprintf("two terms %.10f  one || block %.10f\n", logLik(two), logLik(one)))
three <- frm(bf(y ~ x + (1 | gr(g, dist = "student")) +
                  (0 + x | gr(g, dist = "student")) +
                  (0 + z | gr(g, dist = "student"))), data = dd)
onez <- frm(bf(y ~ x + (1 + x + z || gr(g, dist = "student"))), data = dd)
cat(sprintf("three terms %.10f (conv %d)  one || block %.10f (conv %d)\n",
            logLik(three), three$opt$convergence, logLik(onez),
            onez$opt$convergence))
cat("theta three:", round(three$estimates$theta, 6), "\n")
cat("theta one:  ", round(onez$estimates$theta, 6), "\n")
cat("\n-- two responses, a t term on g in each (frmtmb)\n")
mv <- try_(frm(bf(y ~ x + (1 | gr(g, dist = "student"))) +
                 bf(yz ~ x + (1 | gr(g, dist = "student"))) +
                 set_rescor(FALSE), data = dd, dry_run = "frame"))
if (!is.character(mv)) {
  cat("t_share:", deparse(lapply(mv$re_blocks, `[[`, "t_share")), "\n")
} else cat(mv, "\n")
sink()
