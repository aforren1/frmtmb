# Per-link sanity: brms agreement, round trip, mu_eta vs numDeriv, and
# tapeability of every field.
options(digits = 17)
L <- frmtmb:::frmtmb_links
new <- c("probit", "probit_approx", "cauchit", "softit", "softplus",
         "squareplus", "sqrt", "log1p", "1/mu^2")

# eta ranges each link is meaningful on
rng <- list(probit = c(-2, -0.5, 0.3, 1.7), probit_approx = c(-2, -0.5, 0.3, 1.7),
            cauchit = c(-2, -0.5, 0.3, 1.7), softit = c(-2, -0.5, 0.3, 1.7),
            softplus = c(-2, -0.5, 0.3, 1.7), squareplus = c(-2, -0.5, 0.3, 1.7),
            sqrt = c(0.4, 1.1, 2.3), log1p = c(-0.7, 0.3, 1.2),
            `1/mu^2` = c(0.3, 1.1, 4.0))

cat("== linkinv against brms:::inv_link, linkfun against brms:::link ==\n")
for (nm in new) {
  e <- rng[[nm]]
  bi <- brms:::inv_link(e, nm)
  fi <- L[[nm]]$linkinv(e)
  if (nm == "probit_approx") {
    # brms's R inv_link is pnorm here; its STAN program is Phi_approx.
    bi <- 1 / (1 + exp(-(0.07056 * e^3 + 1.5976 * e)))
  }
  d1 <- max(abs(fi - bi) / pmax(1, abs(bi)))
  mu <- if (nm == "probit_approx") stats::pnorm(e) else bi
  bl <- brms:::link(mu, nm)
  fl <- L[[nm]]$linkfun(mu)
  d2 <- max(abs(fl - bl) / pmax(1, abs(bl)))
  cat(sprintf("%-14s inv %8.2g  fun %8.2g\n", nm, d1, d2))
}

cat("\n== round trip linkfun(linkinv(eta)) == eta ==\n")
for (nm in setdiff(new, "probit_approx")) {
  e <- rng[[nm]]
  if (nm == "sqrt") e <- abs(e)
  rt <- L[[nm]]$linkfun(L[[nm]]$linkinv(e))
  cat(sprintf("%-14s %8.2g\n", nm, max(abs(rt - e))))
}
cat("probit_approx      (not an inverse pair in brms either; skipped)\n")

cat("\n== mu_eta against numDeriv on linkinv ==\n")
for (nm in new) {
  e <- rng[[nm]]
  an <- L[[nm]]$mu_eta(e)
  nd <- vapply(e, function(z) numDeriv::grad(L[[nm]]$linkinv, z), 0)
  cat(sprintf("%-14s %8.2g\n", nm, max(abs(an - nd) / pmax(1e-8, abs(nd)))))
}

cat("\n== every field tapes and differentiates ==\n")
for (nm in names(L)) {
  e0 <- if (nm %in% names(rng)) rng[[nm]][1] else 0.4
  if (nm %in% c("inverse", "1/mu^2", "sqrt", "log")) e0 <- 0.7
  for (f in c("linkfun", "linkinv", "mu_eta", "logit_eta", "log_eta")) {
    if (is.null(L[[nm]][[f]])) next
    x0 <- if (f == "linkfun") {
      if (nm %in% c("logit", "probit", "probit_approx", "cauchit",
                    "cloglog", "softit")) 0.63
      else if (nm == "logm1") 2.4 else if (nm == "power12") 1.4
      else if (nm == "log1p") 0.4 else if (nm == "identity") 0.4 else 0.9
    } else e0
    r <- tryCatch({
      tp <- RTMB::MakeTape(function(p) sum(L[[nm]][[f]](p)), x0)
      g <- as.numeric(tp$jacobian(x0))
      if (is.finite(tp(x0)) && is.finite(g)) "ok" else "NONFINITE"
    }, error = function(err) paste("ERR", conditionMessage(err)))
    if (!identical(r, "ok")) cat(sprintf("%-14s %-10s %s\n", nm, f, r))
  }
}
cat("(nothing above = every field taped with a finite value and gradient)\n")

cat("\n== robust fields reproduce the round trip where it is trusted ==\n")
for (nm in c("probit", "probit_approx", "softit")) {
  e <- c(-3, -0.5, 0, 0.5, 3)
  p <- L[[nm]]$linkinv(e)
  cat(sprintf("%-14s %8.2g\n", nm,
              max(abs(L[[nm]]$logit_eta(e) - log(p / (1 - p))))))
}
e <- c(-3, -0.5, 0, 0.5, 3)
cat(sprintf("%-14s %8.2g\n", "squareplus",
            max(abs(L$squareplus$log_eta(e) - log(L$squareplus$linkinv(e))))))
