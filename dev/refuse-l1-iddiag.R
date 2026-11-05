# gr(g, id = "p", cor = FALSE) against brms 2.23.0's log_prob: the suite's
# translator reads an id group as one block, and frmtmb holds these as
# separate diagonal blocks because nothing correlates them, so the map is
# written here. Also: mi() x cens() in brms reads no latent value (Y at a
# missing row is Inf). Seed 4040.
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
set.seed(4040)
n <- 200
d <- data.frame(x = rnorm(n), w = rnorm(n),
                g = factor(rep(sprintf("g%02d", 1:10), each = 20)))
d$y <- 1 + d$x + rnorm(10)[d$g] + rnorm(n, 0, exp(0.2 * rnorm(10)[d$g]))
fb <- brms::bf(y ~ x + (1 + w | gr(g, id = "p", cor = FALSE)),
               sigma ~ (1 | gr(g, id = "p", cor = FALSE)))
fit <- frm(bf(y ~ x + (1 + w | gr(g, id = "p", cor = FALSE)),
              sigma ~ (1 | gr(g, id = "p", cor = FALSE))), data = d)
ctl <- frm(bf(y ~ x + (1 + w || g), sigma ~ (1 | g)), data = d)
cat("same fit as with no key: logLik", as.numeric(logLik(fit)),
    as.numeric(logLik(ctl)), "\n")
pr <- E$brms_flat_prior(fb, data = d, family = stats::gaussian())
code <- brms::make_stancode(fb, data = d, family = stats::gaussian(), prior = pr)
sdat <- E$brms_standata(fb, data = d, family = stats::gaussian(), prior = pr)
mod <- E$brms_stan_model(code)
sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
est <- fit$estimates
fe <- frmtmb:::fixef_by_dpar(fit)
vc <- VarCorr(fit)$g$sd
print(dimnames(vc))
print(sapply(fit$frame$re_blocks, function(b) b$term_label))
# brms's r_1_k = sd_1[k] * z_1[k, ]: frmtmb's coefficient-space b per column
cv <- frmtmb:::coef_b(fit)
blk <- fit$frame$re_blocks
col_vals <- list()
sds <- numeric(0)
for (bk in blk) {
  B <- matrix(cv[bk$c_idx], nrow = bk$dim)
  S <- exp(est$theta[bk$theta_idx[frmtmb:::block_sd_idx(bk)]])
  for (j in seq_len(bk$dim)) {
    col_vals[[length(col_vals) + 1L]] <- B[j, ]
    sds <- c(sds, S[j])
  }
}
print(sds)
pars <- list(
  b = array(fe$mu[-1], 1L),
  Intercept = unname(fe$mu[1] + fe$mu[2] * mean(d$x)),
  Intercept_sigma = unname(fe$sigma[1]),
  sd_1 = array(sds, 3L),
  z_1 = do.call(rbind, Map(function(v, s) v / s, col_vals, sds)))
up <- rstan::unconstrain_pars(sf, pars)
lp <- rstan::log_prob(sf, up, adjust_transform = FALSE, gradient = FALSE)
ours <- -fit$obj$env$f(fit$obj$env$last.par.best) +
  ncol(pars$z_1) * sum(log(sds))
g <- rstan::grad_log_prob(sf, up, adjust_transform = FALSE)
idx <- E$brms_inner_index(sf, pars, "^z_")
cat("L1 id cor = FALSE: const", lp - ours, "inner grad", max(abs(g[idx])),
    "ours", ours, "\n")
# mi() x cens(): brms's Y at the missing rows
d$ym <- d$y
d$ym[1:5] <- NA
d$ce <- 0L
fb2 <- brms::bf(ym | mi() + cens(ce) ~ x)
s2 <- brms::standata(fb2, data = d)
cat("mi x cens: Y at Jmi", s2$Y[s2$Jmi], " Jmi", s2$Jmi, "\n")
code2 <- as.character(brms::stancode(fb2, data = d))
cat(grep("lpdf|lccdf|lcdf", strsplit(code2, "\n")[[1]], value = TRUE)[1:6],
    sep = "\n")
mod2 <- E$brms_stan_model(code2)
sf2 <- suppressMessages(rstan::sampling(mod2, data = s2, chains = 0))
p2 <- list(Ymi = array(rep(1, 5), 5L), b = array(1, 1L), Intercept = 1,
           sigma = 1)
u2 <- rstan::unconstrain_pars(sf2, p2)
cat("mi x cens: log_prob at a finite point",
    rstan::log_prob(sf2, u2, adjust_transform = FALSE), "\n")
