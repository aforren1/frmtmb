# Lane wt-mvprior, punch round 2: every family with several location
# dpars, core and extensions. Per arm: the location dpars, the b,
# Intercept and sd rows default_prior() lists, and what class b and
# Intercept reach with and without dpar. Extends the reviewer's
# dev/mvprior-review2/r2-extfam.R with rdm(), mixture_mvn() and a
# random effect on hmm().
#   MVPRIOR_ARM=base|lane Rscript dev/mvprior-extfam.R
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-prelude.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
  library(frmtmb.eam)
})
cat("arm", mvprior_arm, "frmtmb from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
flat <- function(x) substr(gsub("[[:space:]]+", " ", x), 1, 260)
reach <- function(design, pl) {
  r <- ns$resolve_priorlist(design, pl)
  pt <- design$frame[["par_template"]]
  nm <- vapply(r$entries, function(e) {
    paste(ns$par_template_names(pt[[e$comp]], e$comp)[e$idx], collapse = "+")
  }, "")
  paste(c(nm, if (length(r$lower)) paste0("lb:", names(r$lower))),
        collapse = " ")
}
run <- function(label, f, fam, d, locs) {
  cat("\n==", label, "==\n")
  des <- tryCatch(ns$prior_design(f, d, fam, list()), error = function(e) {
    cat("MODEL:", flat(conditionMessage(e)), "\n")
    NULL
  })
  if (is.null(des)) return(invisible())
  rs <- des$spec$responses[[1]]
  cat("location dpars:", setdiff(rs$primary_dpars, rs$nlpars), "\n")
  dp <- as.data.frame(default_prior(f, data = d, family = fam))
  dp <- dp[dp$class %in% c("b", "Intercept", "sd"),
           c("class", "coef", "group", "dpar")]
  cat("rows:", paste(apply(dp, 1, paste, collapse = "|"), collapse = " ; "),
      "\n")
  specs <- list(
    b = set_prior("normal(0, 1)", class = "b"),
    Intercept = set_prior("normal(0, 1)", class = "Intercept"),
    sd = set_prior("normal(0, 1)", class = "sd"))
  for (dp1 in locs) {
    specs[[paste("b dpar", dp1)]] <-
      set_prior("normal(0, 1)", class = "b", dpar = dp1)
    specs[[paste("Intercept dpar", dp1)]] <-
      set_prior("normal(0, 1)", class = "Intercept", dpar = dp1)
  }
  for (nm in names(specs)) {
    r <- tryCatch(paste("ACCEPT ->", reach(des, specs[[nm]])),
                  error = function(e) paste("REFUSE:",
                                            flat(conditionMessage(e))))
    cat(sprintf("  %-22s %s\n", nm, r))
  }
}

set.seed(3)
n <- 200
cl <- rbinom(n, 1, 0.4) + 1
pr <- rbind(c(0.85, 0.80, 0.75, 0.90), c(0.15, 0.20, 0.25, 0.10))
Y <- matrix(0L, n, 4)
for (j in 1:4) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
dd <- data.frame(x = rnorm(n))
dd$Y <- Y
run("lca(K = 3), Y ~ x", bf(Y ~ x), lca(K = 3), dd, c("theta1", "theta2"))

set.seed(11)
hd <- do.call(rbind, lapply(1:15, function(id) {
  s <- integer(20)
  s[1] <- 1L
  for (t in 2:20) {
    s[t] <- sample.int(2, 1, prob = if (s[t - 1] == 1) c(.9, .1) else c(.2, .8))
  }
  data.frame(id = id, t = 1:20, x = rnorm(20), y = rnorm(20, c(0, 3)[s]))
}))
run("hmm(K = 2), y ~ x", bf(y ~ x), hmm(K = 2, gaussian(), time = t,
                                        group = id), hd, c("mu1", "mu2"))

set.seed(1)
ld <- lba_simulate(300, v = c(2.4, 1.6, 1.0), A = 0.5, k = 0.4, ndt = 0.2)
ld$x <- rnorm(nrow(ld))
run("lba(3), rt | vint(choice) ~ x", bf(rt | vint(choice) ~ x), lba(3), ld,
    c("v1", "v3"))

set.seed(2)
rd <- rdm_simulate(300, v = c(2.4, 1.6, 1.0), A = 0.5, k = 0.5, ndt = 0.2)
rd$x <- rnorm(nrow(rd))
run("rdm(3), rt | vint(choice) ~ x", bf(rt | vint(choice) ~ x), rdm(3), rd,
    c("v1", "v3"))

set.seed(42)
nm <- 300
clm <- rbinom(nm, 1, 0.4)
Ym <- cbind(rnorm(nm, ifelse(clm == 1, 0, 3)), rnorm(nm, ifelse(clm == 1, 0, 4)))
md <- data.frame(x = rnorm(nm))
md$Y <- Ym
run("mixture_mvn(K = 2, D = 2), Y ~ x", bf(Y ~ x), mixture_mvn(K = 2, D = 2),
    md, character(0))
