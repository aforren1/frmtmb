# Reviewer data and fits (lane ceplot) for the crossed-relabel attack:
# shared by dev/ceplot-rev-crossed.R and dev/ceplot-rev-crossed-brms.R.
hand <- function(fit, n, seed = 1) {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  set.seed(seed)
  M <- matrix(rep(est, each = n) + rnorm(n * length(est), 0, 0.05), n,
              dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  structure(list(stanfit = NULL,
                 draws = cbind(frmtmb.sample:::draws_to_natural(M, fit),
                               lp__ = 0),
                 fit = fit), class = "frmtmb_draws")
}
crossed_data <- function(seed = 49, drop = list(c(1, 1)), ng = 6, nh = 5) {
  set.seed(seed)
  dc <- expand.grid(g = factor(1:ng), h = factor(1:nh), r = 1:4)
  for (d in drop) dc <- dc[!(dc$g == d[1] & dc$h == d[2]), ]
  dc$x <- rnorm(nrow(dc))
  dc$y <- rnorm(nrow(dc), 1 + 0.5 * dc$x + rnorm(ng)[dc$g] +
                  rnorm(nh)[dc$h] +
                  rnorm(ng * nh, 0, 0.7)[as.integer(interaction(dc$g, dc$h))],
                0.5)
  dc
}
shape_fits <- function() {
  out <- list()
  dA <- crossed_data(49)
  out$A <- list(fit = frm(bf(y ~ x + (1 | g) + (1 | h) + (1 | g:h)),
                          family = gaussian(), data = dA),
                cond = list(g = "1", h = "1"))
  dB <- crossed_data(51, drop = list(c(1, 1), c(2, 3)))
  out$B <- list(fit = frm(bf(y ~ x + (1 | g) + (1 | h) + (1 | g:h)),
                          family = gaussian(), data = dB),
                cond = data.frame(g = c("1", "2", "3"), h = c("1", "3", "2")))
  out$D <- list(fit = frm(bf(y ~ x + (1 | g / h) + (1 | h)),
                          family = gaussian(), data = dA),
                cond = list(g = "1", h = "1"))
  dE <- crossed_data(52)
  dE$f <- factor(ifelse(as.integer(dE$g) <= 3, "a", "b"))
  out$E <- list(fit = frm(bf(y ~ x + (1 | gr(g, by = f)) + (1 | h) +
                               (1 | g:h)), family = gaussian(), data = dE),
                cond = list(g = "1", h = "1", f = "a"))
  set.seed(53)
  n <- 360
  dF <- data.frame(g1 = factor(sample(1:6, n, TRUE), levels = 1:6),
                   h = factor(sample(1:5, n, TRUE), levels = 1:5))
  dF$g2 <- factor(sample(1:6, n, TRUE), levels = 1:6)
  dF <- dF[!(dF$g1 == "1" & dF$h == "1"), ]
  dF$x <- rnorm(nrow(dF))
  u <- rnorm(6)
  dF$y <- rnorm(nrow(dF), 1 + 0.5 * dF$x + 0.5 * (u[dF$g1] + u[dF$g2]) +
                  rnorm(5)[dF$h] + rnorm(30, 0, 0.7)[
                    as.integer(interaction(dF$g1, dF$h))], 0.5)
  out$F <- list(fit = frm(bf(y ~ x + (1 | mm(g1, g2)) + (1 | h) +
                               (1 | g1:h)), family = gaussian(), data = dF),
                cond = list(g1 = "1", g2 = "2", h = "1"))
  set.seed(54)
  dG <- expand.grid(g = factor(1:5), h = factor(1:4), k = factor(1:4),
                    r = 1:2)
  dG <- dG[!(dG$g == "1" & dG$h == "1") & !(dG$g == "1" & dG$k == "1"), ]
  dG$x <- rnorm(nrow(dG))
  dG$y <- rnorm(nrow(dG), 1 + 0.5 * dG$x + rnorm(5)[dG$g] + rnorm(4)[dG$h] +
                  rnorm(4)[dG$k] +
                  rnorm(20, 0, 0.7)[as.integer(interaction(dG$g, dG$h))] +
                  rnorm(20, 0, 0.7)[as.integer(interaction(dG$g, dG$k))],
                0.5)
  out$G <- list(fit = frm(bf(y ~ x + (1 | g) + (1 | h) + (1 | k) +
                               (1 | g:h) + (1 | g:k)), family = gaussian(),
                          data = dG),
                cond = list(g = "1", h = "1", k = "1"))
  dH <- crossed_data(55)
  dH$y <- dH$y + 0.4 * dH$x * rnorm(30)[as.integer(interaction(dH$g, dH$h))]
  out$H <- list(fit = frm(bf(y ~ x + (1 | g) + (1 | h) + (1 + x | g:h)),
                          family = gaussian(), data = dH),
                cond = list(g = "1", h = "1"))
  out
}
