source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# M1 structures. Reference: autoscale = FALSE at sx = 1 (well scaled).
# Test: the same data with x multiplied by sx = 1e-3. logLik and fitted()
# are invariant to the rescale, so the default at 1e-3 must match the
# reference; FALSE at 1e-3 shows whether the structure stalls at all.
seed <- as.integer(Sys.getenv("R2_SEED", "21"))
set.seed(seed)
ng <- 20; per <- 15; n <- ng * per
g <- factor(rep(seq_len(ng), each = per))
h <- factor(sample(rep(1:10, length.out = n)))
fb <- factor(rep(c("a", "b"), length.out = ng))[g]
g2 <- factor(sample(rep(seq_len(ng), length.out = n)))
xs <- rnorm(n); z <- rnorm(n)
u0 <- rnorm(ng, 0, 0.8); u1 <- rnorm(ng, 0, 0.5); v1 <- rnorm(10, 0, 0.4)
eta <- 1 + 0.7 * xs + u0[g] + u1[g] * xs + v1[h] * xs + 0.3 * z +
  0.2 * xs * z
y <- eta + rnorm(n) * exp(0.2 * xs)
yp <- rpois(n, exp(0.5 + 0.3 * xs + 0.3 * u1[g] * xs))
base <- data.frame(y, yp, xs, z, g, h, fb, g2)
ctlF <- frmtmb_control(autoscale = FALSE)

cases <- list(
  S1_us = list(f = y ~ x + (1 + x | g)),
  S2_twoblocks = list(f = y ~ x + (1 + x | g) + (0 + x | h)),
  S3_interact = list(f = y ~ x * z + (1 + x | g)),
  S3b_int_slope = list(f = y ~ x * z + (1 + x:z | g)),
  S4_by = list(f = y ~ x + (1 + x | gr(g, by = fb))),
  S5_dblbar = list(f = y ~ x + (1 + x || g)),
  S6_poly = list(f = y ~ poly(x, 2) + (1 + x | g)),
  S6b_Ionly = list(f = y ~ I(x * 1000) + (1 + x | g)),
  S7_idlink = list(f = bf(y ~ x + (1 + x | p | g), sigma ~ x + (1 + x | p | g)),
                   fam = gaussian()),
  S8_reml = list(f = y ~ x + (1 + x | g), REML = TRUE),
  S9_student = list(f = y ~ x + (1 + x | gr(g, dist = "student"))),
  S10_diag = list(f = y ~ x + diag(1 + x | g)),
  S11_homdiag = list(f = y ~ x + homdiag(1 + x | g)),
  S12_cs = list(f = y ~ x + cs(1 + x | g)),
  S13_mm = list(f = y ~ x + (1 + x | mm(g, g2))),
  S14_sigma_slope = list(f = bf(y ~ x, sigma ~ x + (1 + x | g)),
                         fam = gaussian()),
  S15_poisson = list(f = yp ~ x + (1 + x | g), fam = poisson()),
  S16_profile = list(f = y ~ x + (1 + x | g), prof = TRUE)
)
only <- Sys.getenv("R2_CASES", "")
if (nzchar(only)) cases <- cases[strsplit(only, ",")[[1]]]

run <- function(cs, sx, ctl) {
  d <- base; d$x <- d$xs * sx
  ctl <- if (is.null(ctl)) frmtmb_control(profile = isTRUE(cs$prof))
         else frmtmb_control(autoscale = FALSE, profile = isTRUE(cs$prof))
  args <- list(cs$f, data = d, control = ctl, REML = isTRUE(cs$REML))
  if (!is.null(cs$fam)) args$family <- cs$fam
  do.call(fitw, args)
}
plan_z <- function(r) {
  if (inherits(r$fit, "error")) return("error")
  pl <- frmtmb:::autoscale_plan(r$fit$frame)
  if (is.null(pl)) return("no plan")
  paste(vapply(names(pl), function(k) {
    zz <- pl[[k]]$z
    sprintf("%s:cols=%s z=[%s]", k, paste(pl[[k]]$cols, collapse = ","),
            paste(vapply(zz, function(e) sprintf("theta%d/%dcols/s=%.3g",
                  e$theta, length(e$zcols), e$scale), ""), collapse = " "))
  }, ""), collapse = " ; ")
}
llx <- function(r) if (inherits(r$fit, "error"))
  paste("ERR", conditionMessage(r$fit)) else
  sprintf("%.6f c%d", as.numeric(logLik(r$fit)), r$fit$opt$convergence)

for (nm in names(cases)) {
  cs <- cases[[nm]]
  ref <- run(cs, 1, "F")
  dF <- run(cs, 1e-3, "F")
  dD <- run(cs, 1e-3, NULL)
  fitdiff <- if (!inherits(ref$fit, "error") && !inherits(dD$fit, "error"))
    sprintf("%.2e", rel(fitted(dD$fit)[, 1], fitted(ref$fit)[, 1])) else "NA"
  cat(sprintf("%-16s ref=%s | F@1e-3=%s | default@1e-3=%s engaged=%s fitted_rel=%s\n",
              nm, llx(ref), llx(dF), llx(dD), !is.null(dD$tpl), fitdiff))
  cat("   plan@1e-3:", plan_z(dD), "\n")
  if (!is.null(dD$tpl)) cat("   template move:", tpl_move(dD), "\n")
  if (length(dD$warn)) cat("   warn default:", unique(substr(dD$warn, 1, 110)),
                           sep = "\n      ")
}
