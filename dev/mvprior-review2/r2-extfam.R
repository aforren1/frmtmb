# Families with several primary dpars outside categorical/mixture:
# lca (mix slot), hmm (no mix slot), lba (eam), multinomial,
# mixture_mvn. Does class b without dpar refuse or broadcast, per arm?
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-review2/r2-prelude.R")
suppressPackageStartupMessages({library(frmtmb.latent); library(frmtmb.eam)})
ns <- asNamespace("frmtmb")
flat <- function(x) substr(gsub("[[:space:]]+", " ", x), 1, 300)
reach <- function(design, pl) {
  r <- ns$resolve_priorlist(design, pl)
  pt <- design$frame[["par_template"]]
  paste(vapply(r$entries, function(e) paste(ns$par_template_names(pt[[e$comp]], e$comp)[e$idx], collapse = "+"), ""), collapse = " ")
}
run <- function(label, f, fam, d, specs) {
  cat("\n==", label, "==\n")
  des <- tryCatch(ns$prior_design(f, d, fam, list()), error = function(e) {cat("MODEL:", flat(conditionMessage(e)), "\n"); NULL})
  if (is.null(des)) return(invisible())
  cat("primary_dpars:", des$spec$responses[[1]]$primary_dpars, " mix slot:",
      !is.null(des$spec$responses[[1]]$family[["mix"]]), "\n")
  dp <- as.data.frame(default_prior(f, data = d, family = fam))
  dp <- dp[dp$class %in% c("b", "Intercept", "sd"), c("class", "coef", "group", "dpar", "nlpar")]
  cat("default rows b/Intercept/sd:", paste(apply(dp, 1, paste, collapse = "|"), collapse = " ; "), "\n")
  for (s in specs) {
    r <- tryCatch(paste("ACCEPT ->", reach(des, s[[2]])),
                  error = function(e) paste("REFUSE:", flat(conditionMessage(e))))
    cat(sprintf("  %-34s %s\n", s[[1]], r))
  }
}
set.seed(3)
n <- 200
cl <- rbinom(n, 1, 0.4) + 1
pr <- rbind(c(0.85, 0.80, 0.75, 0.90), c(0.15, 0.20, 0.25, 0.10))
Y <- matrix(0L, n, 4); for (j in 1:4) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
dd <- data.frame(x = rnorm(n)); dd$Y <- Y
run("lca K=3, Y ~ x", bf(Y ~ x), lca(K = 3), dd, list(
  list("b", set_prior("normal(0, 1)", class = "b")),
  list("Intercept", set_prior("normal(0, 1)", class = "Intercept")),
  list("b dpar=theta1", set_prior("normal(0, 1)", class = "b", dpar = "theta1")),
  list("b dpar=theta2 coef=x", set_prior("normal(0, 1)", class = "b", dpar = "theta2", coef = "x")),
  list("b lb=0 dpar=theta1", set_prior("", class = "b", dpar = "theta1", lb = 0))))

set.seed(11)
hd <- do.call(rbind, lapply(1:15, function(id) {
  s <- integer(20); s[1] <- 1L
  for (t in 2:20) s[t] <- sample.int(2, 1, prob = if (s[t-1] == 1) c(.9, .1) else c(.2, .8))
  data.frame(id = id, t = 1:20, x = rnorm(20), y = rnorm(20, c(0, 3)[s]))
}))
run("hmm K=2 gaussian, y ~ x", bf(y ~ x), hmm(K = 2, gaussian(), time = t, group = id), hd, list(
  list("b", set_prior("normal(0, 1)", class = "b")),
  list("Intercept", set_prior("normal(0, 1)", class = "Intercept")),
  list("b dpar=mu1", set_prior("normal(0, 1)", class = "b", dpar = "mu1"))))

set.seed(1)
ld <- lba_simulate(300, v = c(2.4, 1.6, 1.0), A = 0.5, k = 0.4, ndt = 0.2)
ld$x <- rnorm(nrow(ld))
run("lba(3), rt | vint(choice) ~ x", bf(rt | vint(choice) ~ x), lba(3), ld, list(
  list("b", set_prior("normal(0, 1)", class = "b")),
  list("Intercept", set_prior("normal(0, 1)", class = "Intercept"))))

set.seed(5)
md <- data.frame(x = rnorm(60))
P <- t(sapply(md$x, function(v) {e <- exp(c(0, 0.5 * v, -0.3 * v)); e / sum(e)}))
md$Ym <- t(sapply(1:60, function(i) rmultinom(1, 20, P[i, ])))
md$N <- 20
run("multinomial, Ym | trials(N) ~ x", bf(Ym | trials(N) ~ x), multinomial(K = 3), md, list(
  list("b", set_prior("normal(0, 1)", class = "b")),
  list("b dpar=mu2", set_prior("normal(0, 1)", class = "b", dpar = "mu2")),
  list("Intercept dpar=mu3", set_prior("normal(0, 1)", class = "Intercept", dpar = "mu3"))))
