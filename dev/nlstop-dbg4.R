# Lane nlstop: item 4705, the data-only Hessian at the ridge stop and the
# prior's gain along its flattest direction, for kt = 5, 20, 50, 80.
source("dev/nlstop-common.R")
ns <- asNamespace("frmtmb")
set.seed(77)
dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(8, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
form2 <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
psd <- as.numeric(Sys.getenv("PSD", "1000"))
pr <- set_prior(sprintf("normal(0.5, %g)", psd), nlpar = "k")
f <- suppressWarnings(suppressMessages(frm(form2, data = dn, prior = pr)))
nm <- ns$outer_par_names(f)
ia <- match("a_(Intercept)", nm); ik <- match("k_(Intercept)", nm)
it <- match("theta_1", nm)
cat("MAP:", format(f$opt$par, digits = 5), "\n")
nlp <- ns$neg_log_prior_fn(ns$fit_prior_entries(f))
for (kt in c(5, 20, 50, 80, 0)) {
  if (kt == 0) { g <- f } else {
  opt <- function(par, fn, gr, lower, upper, control) {
    o <- stats::nlminb(par, fn, gr, control = list(iter.max = 1000,
                                                   eval.max = 2000))
    p <- o$par
    s <- p[ik] / kt
    p[ia] <- p[ia] * s; p[ik] <- kt; p[it] <- p[it] + log(abs(s))
    list(par = p, objective = fn(p), convergence = 0L,
         message = "relative convergence (4)")
  }
  g <- suppressWarnings(suppressMessages(frm(form2, data = dn, prior = pr,
    control = frmtmb_control(optimizer = opt, restarts = 0))))
  }
  p <- g$opt$par
  gr <- drop(g$obj$gr(p))
  H <- ns$fit_outer_hessian(g)$H
  est <- g$estimates
  phi <- function(x) { e <- est; e$beta[] <- x[1:length(e$beta)]
    e$theta[] <- x[-(1:length(e$beta))]; as.numeric(nlp(e)) }
  np <- length(p); Hp <- matrix(0, np, np); h <- 1e-3; f0 <- phi(p)
  for (i in 1:np) for (j in 1:np) {
    ei <- replace(numeric(np), i, h); ej <- replace(numeric(np), j, h)
    Hp[i, j] <- (phi(p + ei + ej) - phi(p + ei - ej) - phi(p - ei + ej) +
                   phi(p - ei - ej)) / (4 * h^2)
  }
  Hd <- H - Hp
  D <- sqrt(abs(diag(H)))
  ev <- eigen(Hd / outer(D, D), symmetric = TRUE)
  evf <- eigen(H / outer(D, D), symmetric = TRUE)
  v <- ev$vectors[, np] / D
  gain <- 0.5 * sum(v * gr)^2 / drop(t(v) %*% Hp %*% v)
  newton <- 0.5 * sum(gr * solve(H, gr))
  cat(sprintf("kt %g short %.4g newton %.4g | Hd eig min/max %.3g  H eig min/max %.3g | prior gain along Hd-flat %.4g\n",
              kt, g$opt$objective - f$opt$objective, newton,
              ev$values[np] / max(abs(ev$values)),
              evf$values[np] / max(abs(evf$values)), gain))
}
