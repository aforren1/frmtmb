# frmtmb.latent's test-prior-dpar.R read which parameters a prior
# reaches through three core internals (prior_design(),
# resolve_priorlist(), par_template_names()). Its new pd_reached() uses
# exported functions only. Here both routes answer the same priors on
# the three latent families, and must agree on each.
#   Rscript dev/smallfix-pd-reached.R <lib;lib>
# Log: dev/smallfix-out/pd-reached.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(a[1], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.latent)})
cat("frmtmb", format(packageVersion("frmtmb")), dirname(find.package("frmtmb")),
    "\n")
old_route <- function(formula, family, data, pl) {
  des <- frmtmb:::prior_design(formula, data, family, list())
  r <- frmtmb:::resolve_priorlist(des, pl)
  pt <- des$frame[["par_template"]]
  sort(unlist(lapply(r$entries, function(e) {
    frmtmb:::par_template_names(pt[[e$comp]], e$comp)[e$idx]
  })))
}
new_route <- function(formula, family, data, pl) {
  u <- frm(formula, family = family, data = data, dry_run = "objective")
  r <- resolve_prior_input(u, pl)
  pt <- par_template(u)
  sort(unlist(lapply(r$entries, function(e) names(pt[[e$comp]])[e$idx])))
}
set.seed(3)
n <- 200
cl <- stats::rbinom(n, 1, 0.4) + 1
pr <- rbind(c(0.85, 0.80, 0.75, 0.90), c(0.15, 0.20, 0.25, 0.10))
Y <- matrix(0L, n, 4)
for (j in 1:4) Y[, j] <- 1L + stats::rbinom(n, 1, pr[cl, j])
dl <- data.frame(x = stats::rnorm(n))
dl$Y <- Y
set.seed(11)
dh <- do.call(rbind, lapply(1:15, function(id) {
  s <- integer(20); s[1] <- 1L
  for (t in 2:20) s[t] <- sample.int(2, 1, prob = if (s[t - 1] == 1)
    c(.9, .1) else c(.2, .8))
  data.frame(id = id, t = 1:20, x = stats::rnorm(20),
             y = stats::rnorm(20, c(0, 3)[s]))
}))
set.seed(42)
cm <- stats::rbinom(n, 1, 0.4)
dm <- data.frame(x = stats::rnorm(n))
dm$Y <- cbind(stats::rnorm(n, ifelse(cm == 1, 0, 3)),
              stats::rnorm(n, ifelse(cm == 1, 0, 4)))
cases <- list(
  list("lca b theta2", bf(Y ~ x), lca(K = 3), dl,
       set_prior("normal(0, 1)", class = "b", dpar = "theta2")),
  list("lca Intercept theta2", bf(Y ~ x), lca(K = 3), dl,
       set_prior("normal(0, 1)", class = "Intercept", dpar = "theta2")),
  list("lca b coef x theta2 + b theta1", bf(Y ~ x), lca(K = 3), dl,
       set_prior("normal(0, 1)", class = "b", coef = "x", dpar = "theta2") +
         set_prior("normal(0, 2)", class = "b", dpar = "theta1")),
  list("hmm b mu1", bf(y ~ x), hmm(K = 2, gaussian(), time = t, group = id),
       dh, set_prior("normal(0, 1)", class = "b", dpar = "mu1")),
  list("hmm Intercept mu2", bf(y ~ x),
       hmm(K = 2, gaussian(), time = t, group = id), dh,
       set_prior("normal(0, 1)", class = "Intercept", dpar = "mu2")),
  list("mvn b mu1d1", bf(Y ~ x), mixture_mvn(K = 2, D = 2), dm,
       set_prior("normal(0, 1)", class = "b", dpar = "mu1d1")),
  list("mvn b mu2d2 + Intercept mu1d2", bf(Y ~ x), mixture_mvn(K = 2, D = 2),
       dm, set_prior("normal(0, 1)", class = "b", dpar = "mu2d2") +
         set_prior("normal(0, 3)", class = "Intercept", dpar = "mu1d2")))
agree <- 0L
for (cs in cases) {
  o <- old_route(cs[[2]], cs[[3]], cs[[4]], cs[[5]])
  w <- new_route(cs[[2]], cs[[3]], cs[[4]], cs[[5]])
  same <- identical(o, w)
  agree <- agree + same
  cat(sprintf("%-32s old: %-28s new: %-28s %s\n", cs[[1]],
              paste(o, collapse = ","), paste(w, collapse = ","),
              if (same) "SAME" else "DIFFER"))
}
cat("agree", agree, "of", length(cases), "\n")
