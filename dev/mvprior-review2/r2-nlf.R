# Item 5: nlf(sigma ~ a + b * z) rows and targeting, on reviewer data.
who <- Sys.getenv("R2_WHO", "lane")
Sys.setenv(R2_ARM = if (who == "base") "base" else "lane")
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-review2/r2-prelude.R")
if (who == "brms") suppressPackageStartupMessages(library(brms))
ns <- asNamespace(if (who == "brms") "brms" else "frmtmb")
set.seed(71)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n))
d$y <- 1 + 0.5 * d$x + rnorm(n, 0, exp(0.2 + 0.3 * d$z))
f <- ns$bf(y ~ x) + ns$nlf(sigma ~ a + b * z) + ns$lf(a ~ 1, b ~ 1 + w)
dp <- as.data.frame(ns$default_prior(f, data = d, family = gaussian()))
dp[is.na(dp)] <- ""
print(dp[, c("prior", "class", "coef", "group", "resp", "dpar", "nlpar")], row.names = FALSE)
specs <- list(
  "b nlpar=a" = function() ns$set_prior("normal(0, 1)", class = "b", nlpar = "a"),
  "b nlpar=b coef=w" = function() ns$set_prior("normal(0, 1)", class = "b", nlpar = "b", coef = "w"),
  "b nlpar=b" = function() ns$set_prior("normal(0, 1)", class = "b", nlpar = "b"),
  "b dpar=sigma nlpar=a" = function() ns$set_prior("normal(0, 1)", class = "b", dpar = "sigma", nlpar = "a"),
  "b dpar=a" = function() ns$set_prior("normal(0, 1)", class = "b", dpar = "a"),
  "b (location)" = function() ns$set_prior("normal(0, 1)", class = "b"),
  "Intercept nlpar=a" = function() ns$set_prior("normal(0, 1)", class = "Intercept", nlpar = "a"))
for (nm in names(specs)) {
  r <- tryCatch({
    v <- as.data.frame(ns$validate_prior(specs[[nm]](), f, data = d, family = gaussian()))
    v[is.na(v)] <- ""
    u <- v[v$source == "user", ]
    paste("ACCEPT:", paste(u$class, u$coef, u$dpar, u$nlpar, sep = "|", collapse = " ; "))
  }, error = function(e) paste("REFUSE:", substr(gsub("[[:space:]]+", " ", conditionMessage(e)), 1, 200)))
  cat(sprintf("%-22s %s\n", nm, r))
}
if (who != "brms") {
  des <- ns$prior_design(f, d, gaussian(), list())
  pt <- des$frame[["par_template"]]
  for (nm in c("b nlpar=a", "b nlpar=b coef=w", "b nlpar=b")) {
    r <- ns$resolve_priorlist(des, specs[[nm]]())
    cat(nm, "reaches:", vapply(r$entries, function(e) paste(ns$par_template_names(pt[[e$comp]], e$comp)[e$idx], collapse = "+"), ""), "\n")
  }
  # a fit: a tight prior on b_w must pull that one coefficient to 0
  f0 <- suppressWarnings(suppressMessages(frm(f, data = d, family = gaussian())))
  f1 <- suppressWarnings(suppressMessages(frm(f, data = d, family = gaussian(),
          prior = set_prior("normal(0, 0.001)", class = "b", nlpar = "b", coef = "w"))))
  cat("no prior:\n"); print(fixef(f0)[, 1])
  cat("tight prior on nlpar b coef w:\n"); print(fixef(f1)[, 1])
}
