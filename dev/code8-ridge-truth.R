# Part C: which parameters does the ridge of c0 + exp(a)^k reach? The
# invariance (a0, k, b) -> (s a0, k / s, s b) leaves a0 * k, k * sd(b) and
# everything else unchanged, so fixing k at the fit's own value removes
# the ridge and nothing else. The fixed-k fit's standard errors of c0,
# c0_x and sigma are then the full model's, and a0 and theta_1 are the
# ridge's. Also: the profile of the objective along each parameter
# (re-optimizing all others), at steps of 0.25, 0.5 and 1 times the
# fixed-k standard error (or 0.5 on the log scale).
#   CODE8_ARM=before Rscript dev/code8-ridge-truth.R
source("dev/code8-common.R")
options(width = 200)
set.seed(77)
dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(8, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
form <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
f <- quiet(frm(form, data = dn))
p <- f$opt$par
nm <- fns$outer_par_names(f)
cat("full fit objective", format(f$opt$objective, digits = 12), "\n")
print(setNames(p, nm))
k0 <- p[match("k_(Intercept)", nm)]
dn$kk <- k0
fk <- quiet(frm(bf(yn ~ c0 + exp(a)^kk, c0 ~ 1 + x, a ~ 1 + (1 | g),
                   nl = TRUE), data = dn,
                start = list(beta = unname(p[nm %in% c("c0_(Intercept)", "c0_x",
                                                "a_(Intercept)")]))))
cat("fixed-k objective", format(fk$opt$objective, digits = 12), "\n")
sk <- quiet(fns$sdr_of(fk))
print(data.frame(name = fns$outer_par_names(fk), est = fk$opt$par,
                 se = sqrt(diag(as.matrix(sk$cov.fixed)))), digits = 6)
cat("fixed-k se_lost:", names(sk$se_lost), "\n")
# profiles along each parameter of the full model
obj <- f$obj
prof1 <- function(j, t) {
  q0 <- p; q0[j] <- p[j] + t
  fn <- function(q) { z <- q0; z[-j] <- q; as.numeric(obj$fn(z)) }
  gr <- function(q) { z <- q0; z[-j] <- q; drop(obj$gr(z))[-j] }
  o <- tryCatch(quiet(nlminb(p[-j], fn, gr, control = list(
    iter.max = 1000, eval.max = 1000, rel.tol = 1e-14))),
    error = function(e) list(objective = NA))
  o$objective - f$opt$objective
}
for (j in seq_along(p)) {
  st <- c(0.1, 0.2, 0.4)
  loss <- vapply(st, function(t) prof1(j, t), 0)
  lossm <- vapply(st, function(t) prof1(j, -t), 0)
  cat(sprintf("%-16s step %s  loss+ %s  loss- %s\n", nm[j],
              paste(st, collapse = ","), paste(format(loss, digits = 3),
                                               collapse = ","),
              paste(format(lossm, digits = 3), collapse = ",")))
}
cat("DONE\n")
