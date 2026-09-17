# Reviewer, lane wt-conditions (recheck 1): designs that frmtmb accepts,
# run on one library arm, to find a false alarm of the round's new
# refusals: check_frame_variables() in assemble_frame(),
# check_newdata_frame() before each newdata model.frame(), and
# prior_dist_params. Each case records "OK" or the error class and
# message; the base and lane records are then compared by case.
#   Rscript dev/conditions-rev-falsealarm.R base|lane   (seed 20260917)
arm <- commandArgs(TRUE)[1L]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/conditions-lib"
            else "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE, warn = 1)
set.seed(20260917)
n <- 60
d <- data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n),
                g = factor(rep(1:10, each = 6)),
                g1 = factor(sample(1:8, n, TRUE)),
                g2 = factor(sample(1:8, n, TRUE)),
                f = factor(rep(c("a", "b", "c"), 20)),
                w = runif(n, 0.5, 1.5), expo = runif(n, 1, 3),
                w1 = 0.5, w2 = 0.5, t = rep(1:6, 10),
                loc = factor(rep(1:6, 10)))
d$k <- rbinom(n, 10, 0.4); d$ntr <- 10L
d$cnt <- rpois(n, 2 * d$expo)
d$cc <- sample(c(0, 1), n, TRUE, prob = c(0.8, 0.2))
d$ylo <- d$y - 0.5; d$yhi <- d$y + 0.5
d$ym <- ifelse(seq_len(n) %% 9 == 0, NA, d$y)
d$xm <- ifelse(seq_len(n) %% 8 == 0, NA, d$x)
d$sdx <- 0.2; d$sdy <- 0.3
d$ord <- factor(sample(1:4, n, TRUE), ordered = TRUE)
d$ordx <- sample(0:3, n, TRUE)
d$yb <- rbinom(n, 1, 0.5)
d$M <- cbind(d$x, d$z)
d$fc <- as.character(d$f)
d$dt <- as.Date("2020-01-01") + seq_len(n)
d$yp <- 2 * exp(-0.4 * abs(d$x)) + rnorm(n, 0, 0.1)
d$ych <- sample(c("p", "q", "r"), n, TRUE)
A <- diag(10); dimnames(A) <- list(levels(d$g), levels(d$g))
Wadj <- matrix(0, 6, 6); for (i in 1:5) Wadj[i, i + 1] <- Wadj[i + 1, i] <- 1
dimnames(Wadj) <- list(levels(d$loc), levels(d$loc))
xe <- rnorm(n)                      # a variable in the formula environment
kpow <- 2                           # a scalar constant in the environment

# Second batch (recheck 1): data environments with parents, scalar and
# short addition-term variables, structural matrices in data2, matrix
# data. Same record format as dev/conditions-rev-falsealarm.R.
cases <- list()
C <- function(name, expr) cases[[name]] <<- substitute(expr)
F <- function(...) frm(..., dry_run = "frame")
C("data_env_parent_var", local({
  pe <- new.env(); assign("xq", rnorm(60), envir = pe)
  de <- new.env(parent = pe); assign("y", d$y, envir = de); assign("x", d$x, envir = de)
  F(y ~ x + xq, data = de) }))
C("data_env_parent_global", local({
  de <- new.env(parent = globalenv()); assign("y", d$y, envir = de)
  F(y ~ xe, data = de) }))
C("data_list_env_var", F(y ~ x + xe, data = as.list(d[, c("y", "x")])))
C("trials_env_scalar", (function() { nt <- 10L; F(k | trials(nt) ~ x, data = d, family = binomial()) })())
C("weights_env_scalar", (function() { ws <- 2; F(y | weights(ws) ~ x, data = d) })())
C("se_env_scalar", (function() { s0 <- 0.3; F(y | se(s0) ~ x, data = d) })())
C("cens_env_scalar", (function() { c0 <- 0; F(y | cens(c0) ~ x, data = d) })())
C("fcor_data2", F(y ~ x + fcor(Vm), data = d, data2 = list(Vm = diag(60))))
C("sar_data2", F(y ~ x + sar(Wm), data = d[1:6, ], data2 = list(Wm = { m <- matrix(0, 6, 6); for (i in 1:5) m[i, i + 1] <- m[i + 1, i] <- 1; m / rowSums(m) })))
C("car_env", F(y ~ car(Wadj, gr = loc), data = d))
C("ar_cov", F(y ~ x + ar(time = t, gr = g, cov = TRUE), data = d))
C("arma_cov", F(y ~ x + arma(time = t, gr = g, cov = TRUE), data = d))
C("te_smooth", F(y ~ te(x, z), data = d))
C("gp_k", F(y ~ gp(x, k = 10), data = d))
C("s_fs", F(y ~ s(x, g, bs = "fs", k = 4), data = d))
C("matrix_data", F(y ~ x, data = as.matrix(d[, c("y", "x")])))
C("unnamed_list_data", F(y ~ x, data = unname(as.list(d[, c("y", "x")]))))
C("factor_env_var", (function() { fe <- factor(rep(c("u", "v"), 30)); F(y ~ fe, data = d) })())
C("matrix_env_var", (function() { me <- cbind(rnorm(60), rnorm(60)); F(y ~ me, data = d) })())
C("date_env_var", (function() { de <- as.Date("2020-01-01") + 1:60; suppressMessages(F(y ~ de, data = d)) })())
C("env_var_na_rows", (function() { xv <- ifelse(1:60 %% 5 == 0, NA, rnorm(60)); suppressMessages(F(y ~ xv, data = d)) })())
C("interaction_env", F(y ~ x:xe + f * xe, data = d))
C("offset_arg", F(cnt ~ x, offset = log(expo), data = d, family = poisson()))
C("nl_data_col_in_body_and_env", (function() { z <- 1; F(bf(yp ~ a * exp(-b * abs(z)), a ~ 1, b ~ 1, nl = TRUE), data = d) })())
C("hurdle_dpar", F(bf(cnt ~ x, hu ~ z), data = d, family = hurdle_poisson()))
C("zi_dpar", F(bf(cnt ~ x, zi ~ z), data = d, family = zero_inflated_poisson()))
C("mixture_fam", F(bf(y ~ x), data = d, family = mixture(gaussian, gaussian)))
C("resp_transform", F(log(expo) ~ x, data = d))
C("resp_cbind_mv", F(bf(mvbind(y, z) ~ x) + set_rescor(TRUE), data = d))
C("resp_scale", F(scale(y) ~ x, data = d))
C("rhs_scale", F(y ~ scale(x), data = d))
C("rhs_in", F(y ~ I(f %in% c("a", "b")), data = d))

res <- lapply(names(cases), function(nm) {
  out <- tryCatch({
    v <- suppressWarnings(suppressMessages(eval(cases[[nm]])))
    "OK"
  }, error = function(e) paste0("ERROR [", class(e)[1L], "] ",
                               gsub("[[:space:]]+", " ", conditionMessage(e))))
  data.frame(case = nm, result = out, stringsAsFactors = FALSE)
})
out <- do.call(rbind, res)
saveRDS(list(cases = out),
        sprintf("C:/Users/adf44/source/r/frmtmb-wt-conditions/dev/conditions-rev-log/falsealarm2-%s.rds", arm))
cat("cases", nrow(out), " OK", sum(out$result == "OK"), "\n")
