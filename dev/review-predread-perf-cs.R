# Item 3: the layered route's own cost, memory first: frm_linpred(se.fit = TRUE) and
# frm_lp_basis() on a cs() fit in sample and on new data, against the shared part.
# Seed 8011. after only, since 0.78.0 has no such route.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(8011)
n <- 60000
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(sample(sprintf("g%03d", 1:300), n, TRUE)))
d$ord <- as.integer(cut(d$x + ifelse(d$z > 0, 1.2, -1.2) * rlogis(n) + rnorm(300, 0, 0.4)[d$g],
                        c(-Inf, -0.5, 0.5, 1.5, Inf), labels = FALSE))
fit <- quiet(frm(bf(ord ~ x + cs(z) + (1 | g)), family = acat(), data = d))
tm <- function(lab, e) {
  g0 <- gc(reset = TRUE)
  t <- system.time(v <- quiet(e))[["elapsed"]]
  g1 <- gc()
  cat(sprintf("%-38s %7.2f s  peak %8.1f MB  obj %8.1f MB\n", lab, t, g1[2, 6], as.numeric(object.size(v)) / 1e6))
}
tm("frm_linpred link in sample", frm_linpred(fit))
tm("frm_linpred dpar mu in sample", frm_linpred(fit, dpar = "mu"))
tm("frm_linpred dpar mu se in sample", frm_linpred(fit, dpar = "mu", se.fit = TRUE))
tm("frm_linpred layers se, 5000 rows", frm_linpred(fit, newdata = d[1:5000, ], se.fit = TRUE))
tm("frm_linpred layers se, 20000 rows", frm_linpred(fit, newdata = d[1:20000, ], se.fit = TRUE))
tm("frm_lp_basis layers, 5000 rows", frm_lp_basis(fit, newdata = d[1:5000, ]))
tm("frm_lp_basis dpar mu, 5000 rows", frm_lp_basis(fit, newdata = d[1:5000, ], dpar = "mu"))
tm("fitted linear, 20000 rows", fitted(fit, newdata = d[1:20000, ], scale = "linear"))
tm("frm_curve 500 rows", frm_curve(fit, newdata = d[1:500, ], seed = 1))
tm("frm_linpred layers se in sample", frm_linpred(fit, se.fit = TRUE))
