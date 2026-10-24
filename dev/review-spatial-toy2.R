# Reviewer: myar against core's ar1 on post-fit methods, on a design
# with two rows per cell so sigma is identified; and the compat rows a
# registered structure gets with no rules of its own.
.libPaths(c("C:/Users/adf44/source/r/wt-spatial-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
ar_S <- function(theta, d) {
  rho <- tanh(theta[2])
  D <- abs(outer(seq_len(d), seq_len(d), "-"))
  exp(2 * theta[1]) * rho^D
}
frmtmb_register_covstruct(
  "myar",
  lpdf = function(b, theta, blk) {
    d <- blk[["dim"]]
    rho <- tanh(theta[2])
    D <- abs(outer(seq_len(d), seq_len(d), "-"))
    S <- RTMB::matrix(exp(2 * theta[1]) * rho^D, d, d)
    dim(b) <- c(d, length(b) %/% d)
    sum(RTMB::dmvnorm(t(b), 0, S, log = TRUE))
  },
  vcov = function(theta, blk) {
    V <- ar_S(theta, blk[["dim"]])
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  npar = 2L, sd_idx = 1L)
set.seed(42)
G <- 40; d <- 3
dat <- expand.grid(rep = 1:2, t = factor(seq_len(d)), g = factor(seq_len(G)))
S <- ar_S(c(log(0.8), atanh(0.6)), d)
B <- t(chol(S)) %*% matrix(rnorm(d * G), d)
dat$x <- rnorm(nrow(dat))
dat$y <- 1 + 0.5 * dat$x + B[cbind(as.integer(dat$t), as.integer(dat$g))] +
  rnorm(nrow(dat), sd = 0.5)
f1 <- frm(y ~ x + myar(0 + t | g), data = dat)
f0 <- frm(y ~ x + ar1(0 + t | g), data = dat)
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-8))
cat("logLik rel diff", rel(as.numeric(logLik(f1)), as.numeric(logLik(f0))),
    "\n")
cat("sigma", sigma(f1), sigma(f0), "\n")
nd <- data.frame(t = factor(1:3, levels = 1:3), g = factor("new"), x = 0)
nd2 <- rbind(nd, data.frame(t = factor(1:3, levels = 1:3), g = factor("1"),
                            x = 0))
for (pe in c(FALSE, TRUE)) {
  set.seed(5); a <- predict(f1, newdata = nd2, allow_new_levels = TRUE,
               propagate_error = pe)
  set.seed(5); b <- predict(f0, newdata = nd2, allow_new_levels = TRUE,
               propagate_error = pe)
  cat("predict new+old levels propagate_error =", pe, ": myar\n")
  print(a)
  cat("  ar1\n")
  print(b)
}
s1 <- simulate(f1, nsim = 2, seed = 3, newdata = nd2, allow_new_levels = TRUE)
s0 <- simulate(f0, nsim = 2, seed = 3, newdata = nd2, allow_new_levels = TRUE)
cat("simulate newdata myar vs ar1 max rel diff:", rel(as.matrix(s1),
                                                     as.matrix(s0)), "\n")
s1 <- simulate(f1, nsim = 3, seed = 3)
s0 <- simulate(f0, nsim = 3, seed = 3)
cat("simulate myar vs ar1 max rel diff:", rel(as.matrix(s1), as.matrix(s0)),
    "\n")
print(confint_varcorr(f1)); print(confint_varcorr(f0))
cat("variables myar:", variables(f1), "\n")
cat("variables ar1 :", variables(f0), "\n")
ce1 <- conditional_effects(f1, "x")[[1]]
ce0 <- conditional_effects(f0, "x")[[1]]
cat("CE estimate rel diff", rel(ce1$estimate__, ce0$estimate__),
    " se rel diff", rel(ce1$se__, ce0$se__), "\n")
ce1 <- conditional_effects(f1, "x", re_formula = NULL)[[1]]
ce0 <- conditional_effects(f0, "x", re_formula = NULL)[[1]]
cat("CE re_formula NULL estimate rel diff", rel(ce1$estimate__, ce0$estimate__),
    " se rel diff", rel(ce1$se__, ce0$se__), "\n")
r1 <- ranef(f1)$g; r0 <- ranef(f0)$g
cat("ranef rel diff", rel(r1, r0), "\n")
cat("vcov_cluster myar:\n")
print(tryCatch(dim(vcov_cluster(f1, cluster = dat$g)), error = conditionMessage))
cmp <- frm_compat("myar")
cat("compat rows for myar:", nrow(cmp), "\n")
print(table(cmp$status))
sel <- cmp$feature_b %in% c("quadrature", "importance", "gr_by", "student_t",
                            "dist_student", "re_id", "mm", "lkj", "REML",
                            "frm_sample") |
  grepl("quad|importance|student|by|mm|\\|ID\\||id", cmp$feature_b)
print(cmp[sel, c("feature_a", "feature_b", "status")])
fs <- frm_compat_features()
print(fs[fs$name %in% c("myar", "ar1", "quadrature", "importance"), ])
