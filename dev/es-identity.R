# esicar lane: the row 19c identity, run standalone before it goes in
# the tier file.
lib_es <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/es-lib"
.libPaths(c(lib_es, "C:/Users/adf44/AppData/Local/R/win-library/4.6", .libPaths()))
library(frmtmb)
wt <- "C:/Users/adf44/source/r/frmtmb-wt-esicar"
Sys.setenv(FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"),
           FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
env <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(wt, "tests/testthat/helper-brms.R"), envir = env)
for (nm in ls(env)) if (is.function(env[[nm]])) environment(env[[nm]]) <- env

w <- matrix(0, 16, 16)
g <- expand.grid(r = 1:4, c = 1:4)
for (i in 1:16) for (j in 1:16) {
  if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) w[i, j] <- 1
}
dimnames(w) <- list(paste0("L", 1:16), paste0("L", 1:16))
set.seed(42)
kmat <- diag(rowSums(w)) - w + matrix(1 / (1e-3 * 16)^2, 16, 16)
phi <- 1.2 * drop(crossprod(chol(solve(kmat)), rnorm(16)))
loc <- factor(rep(rownames(w), each = 6), levels = rownames(w))
d <- data.frame(loc = loc, x = rnorm(length(loc)))
d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] + rnorm(nrow(d), 0, 0.5)

bform <- brms::bf(y ~ x + car(W, gr = loc, type = "esicar"))
fit <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + gaussian(),
           data = d, data2 = list(W = w))
prior <- env$brms_flat_prior(bform, data = d, family = gaussian(),
                             data2 = list(W = w))
code <- brms::make_stancode(bform, data = d, family = gaussian(),
                            prior = prior, data2 = list(W = w))
sdat <- env$brms_standata(bform, data = d, family = gaussian(),
                          prior = prior, data2 = list(W = w))
cat("stan par names:", paste(env$brms_stan_par_names(code), collapse = ", "),
    "\n")
pars <- env$stan_pars_from_fit(fit, sdat, code)
cat("mapped names:", paste(names(pars), collapse = ", "), "\n")
cat("logJ:", attr(pars, "logJ"), "\n")
cat("zcar length:", length(pars$zcar), " Nloc:", sdat$Nloc, "\n")
cat("predicted const:", env$brms_car_const(sdat, "esicar"), "\n")

cat("compiling / loading the esicar program ...\n")
t0 <- Sys.time()
res <- env$brms_lp_check(bform, gaussian(), d, fit, joint = TRUE,
                         const = env$brms_car_const(sdat, "esicar"),
                         data2 = list(W = w))
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
    "s\n")
cat(sprintf("brms lp        %.15f\n", res$lp))
cat(sprintf("frmtmb joint   %.15f\n", res$ours))
cat(sprintf("measured const %.15f\n", res$measured_const))
cat(sprintf("predicted      %.15f\n", env$brms_car_const(sdat, "esicar")))
cat(sprintf("IDENTITY RESIDUAL %.6e\n",
            res$measured_const - env$brms_car_const(sdat, "esicar")))
cat(sprintf("max |grad| on inner pars %.3e\n", res$max_grad))
