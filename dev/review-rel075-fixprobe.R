# Review short check of the release fixes: emmeans refusal, intact
# emmeans output against 0.74.0, frmtmb_covstruct_info() of residual
# terms, the sar internal row, the errorsar wording, other residual
# terms unaffected. Library given as argument 1.
#   Rscript dev/review-rel075-fixprobe.R <lib> <out.rds>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(frmtmb))
v13 <- packageVersion("frmtmb") >= "0.75.0"
cat("frmtmb", format(packageVersion("frmtmb")), "\n")
out <- list()
res <- function(expr) tryCatch({ force(expr); "ANSWERED" }, error = function(e)
  paste0("ERROR [", class(e)[1], "]: ", conditionMessage(e)))
show <- function(label, x) { cat("\n==", label, "\n"); print(x) }
# ---- 1. emmeans ----
set.seed(1)
d <- data.frame(g = factor(rep(1:20, each = 5)), x = stats::rnorm(100),
                f = factor(rep(c("a", "b"), 50)))
d$y <- 1 + d$x + (d$f == "b") + stats::rnorm(20)[d$g] + stats::rnorm(100, 0, 0.5)
fit <- frm(y ~ x + f + (1 | g), data = d)
out$emm <- summary(emmeans::emmeans(fit, ~ f))
out$rg <- summary(emmeans::ref_grid(fit))
out$pairs <- summary(pairs(emmeans::emmeans(fit, ~ f)))
out$trend <- summary(emmeans::emtrends(fit, ~ f, var = "x"))
fam <- fit$spec$responses[[1L]]$family
attr(fam, "frmtmb_package") <- "frmtmbabsenttestpkga"
fit$spec$responses[[1L]]$family <- fam
fit$frame$spec <- fit$spec
fit$cache <- new.env()
show("ref_grid, missing family package", res(emmeans::ref_grid(fit)))
show("emmeans, missing family package", res(emmeans::emmeans(fit, ~ f)))
saveRDS(out, a[2])
if (!v13) quit(save = "no")
suppressMessages(library(frmtmb.spatial))
# ---- 2. covstruct_info ----
for (n in c("sar", "fcor", "car", "ar1", "us", "nope"))
  show(paste0("frmtmb_covstruct_info('", n, "')"),
       res(print(unlist(frmtmb_covstruct_info(n)))))
# ---- 3. the internal row ----
k <- 6
gg <- expand.grid(r = seq_len(k), c = seq_len(k))
A <- outer(seq_len(k^2), seq_len(k^2), function(i, j)
  as.numeric(abs(gg$r[i] - gg$r[j]) + abs(gg$c[i] - gg$c[j]) == 1))
W <- A / rowSums(A); n <- nrow(W)
set.seed(20261024)
d2 <- data.frame(x = runif(n, 0, 1.5), z = rnorm(n), t = rep(1:6, 6),
                 s = rep(1:6, each = 6))
d2$y <- as.vector(solve(diag(n) - 0.4 * W, 1 + 0.5 * exp(d2$x) + rnorm(n, 0, 0.5)))
V <- 0.5^abs(outer(1:n, 1:n, "-"))
tab <- function(p) as.data.frame(p)[, c("prior", "class", "coef", "lb", "ub")]
for (ty in c("lag", "error")) {
  f <- as.formula(paste0("y ~ exp(x) + sar(W, type = '", ty, "')"))
  show(paste("default_prior", ty), tab(default_prior(f, data = d2, data2 = list(W = W))))
  show(paste("validate_prior theta thetaac_1,", ty),
       res(validate_prior(set_prior("normal(0, 1)", class = "theta", coef = "thetaac_1"),
                          f, data = d2, data2 = list(W = W))))
  show(paste("frm theta thetaac_1,", ty),
       res(frm(f, data = d2, data2 = list(W = W),
               prior = set_prior("normal(0, 1)", class = "theta", coef = "thetaac_1"))))
}
fl <- frm(y ~ exp(x) + sar(W), data = d2, data2 = list(W = W))
suppressMessages(library(frmtmb.sample))
show("frm_sample theta thetaac_1 on sar",
     res(suppressMessages(frm_sample(fl, chains = 1, iter = 100, warmup = 50, refresh = 0,
       prior = set_prior("normal(0, 1)", class = "theta", coef = "thetaac_1")))))
show("frm_sample lagsar prior still taken",
     res(suppressWarnings(suppressMessages(frm_sample(fl, chains = 1, iter = 100, warmup = 50,
       refresh = 0, prior = set_prior("normal(0, 0.5)", class = "lagsar"))))))
show("validate_prior lagsar on error SAR",
     res(validate_prior(set_prior("normal(0, 0.5)", class = "lagsar"),
                        y ~ x + sar(W, type = "error"), data = d2, data2 = list(W = W))))
show("validate_prior lagsar on a model with no residual term",
     res(validate_prior(set_prior("normal(0, 0.5)", class = "lagsar"), y ~ x, data = d2)))
# other residual terms
show("default_prior ar()", tab(default_prior(y ~ x + ar(t, s), data = d2)))
show("frm ar() with theta thetaac_1 prior",
     res(frm(y ~ x + ar(t, s), data = d2,
             prior = set_prior("normal(0, 1)", class = "theta", coef = "thetaac_1"))))
show("default_prior cosy()", tab(default_prior(y ~ x + cosy(t, s), data = d2)))
show("default_prior fcor()", tab(default_prior(y ~ x + fcor(V), data = d2, data2 = list(V = V))))
show("fcor fit logLik", as.numeric(logLik(frm(y ~ x + fcor(V), data = d2, data2 = list(V = V)))))
