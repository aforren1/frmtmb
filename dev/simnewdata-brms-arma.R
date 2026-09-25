# brms 2.23.0 on a residual correlation at newdata, the case the fit
# method now rebuilds: ar(time, gr = g, cov = TRUE), gaussian. Two
# newdata rows of one group at adjacent times, one row of another
# group, and a time the fit never saw.
#   Rscript dev/simnewdata-brms-arma.R > dev/simnewdata-log/brms-arma.txt 2>&1
source("dev/simnewdata-prelude.R")
suppressMessages(library(brms))
set.seed(8)
n_g <- 40
tt <- 1:6
d <- expand.grid(time = tt, g = factor(seq_len(n_g)))
d$y <- 1 + unlist(lapply(seq_len(n_g), function(i) {
  as.vector(stats::arima.sim(list(ar = 0.7), length(tt)))
}))
fit <- brm(y ~ 1 + ar(time = time, gr = g, cov = TRUE), data = d,
           chains = 1, iter = 1500, refresh = 0, seed = 1, silent = 2)
print(summary(fit)$cor_pars)
nd <- data.frame(time = c(2, 3, 2), g = factor(c("a", "a", "b")),
                 y = 0)
pp <- posterior_predict(fit, newdata = nd, allow_new_levels = TRUE)
cat(sprintf("rows 1-2 (one group, lag 1): cor %.3f; rows 1-3 (two groups): cor %.3f; draws %d\n",
            stats::cor(pp[, 1], pp[, 2]), stats::cor(pp[, 1], pp[, 3]),
            nrow(pp)))
nd2 <- data.frame(time = c(2, 9), g = factor(c("a", "a")), y = 0)
r <- tryCatch({
  p2 <- posterior_predict(fit, newdata = nd2, allow_new_levels = TRUE)
  sprintf("answers: cor of the two rows %.3f", stats::cor(p2[, 1], p2[, 2]))
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("time 9, never seen:", r, "\n")
# and frmtmb on the same data
suppressMessages(library(frmtmb))
ff <- frmtmb::frm(bf(y ~ 1 + ar(time, gr = g, cov = TRUE)), data = d)
s <- simulate(ff, nsim = 2000, seed = 9, newdata = nd[, c("time", "g")])
m <- as.matrix(s)
cat(sprintf("frmtmb: rows 1-2 cor %.3f; rows 1-3 cor %.3f; fitted R[2,3] %.3f\n",
            stats::cor(m[1, ], m[2, ]), stats::cor(m[1, ], m[3, ]),
            frmtmb:::autocor_matrix(ff)[2, 3]))
r <- tryCatch(simulate(ff, nsim = 2, newdata = nd2[, c("time", "g")]),
              error = function(e) paste("ERROR:", conditionMessage(e)))
cat("frmtmb, time 9:", if (is.character(r)) r else "answers", "\n")
