# Reviewer: what a user is told when a dispersion runs to its limit.
# negbinomial / negbinomial2 on Poisson data, beta_binomial on binomial
# data, dirichlet_multinomial on multinomial data. Seeds 1..8 each.
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
grab <- function(expr) {
  msgs <- character()
  m <- withCallingHandlers(expr,
    warning = function(w) {
      msgs <<- c(msgs, paste("W:", substr(conditionMessage(w), 1, 70)))
      invokeRestart("muffleWarning")
    },
    message = function(w) {
      msgs <<- c(msgs, paste("M:", substr(conditionMessage(w), 1, 70)))
      invokeRestart("muffleMessage")
    })
  list(m = m, msgs = if (length(msgs)) paste(unique(msgs), collapse = " | ")
       else "(none)")
}
disp_row <- function(fam, r, dp) {
  m <- r$m
  ci <- confint(m)
  i <- grep(paste0("^", dp, "_"), rownames(ci))
  data.frame(fam = fam, eta = unname(ci[i, "est"]),
             ci_w = unname(ci[i, "upr"] - ci[i, "lwr"]),
             conv = m$opt$convergence, msgs = r$msgs)
}
out <- list()
for (seed in 1:8) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = rnorm(n))
  d$y <- rpois(n, exp(1 + 0.3 * d$x))
  d$N <- sample(5:40, n, TRUE)
  d$yb <- rbinom(n, d$N, plogis(0.2 + 0.4 * d$x))
  E <- cbind(0, 0.3 + 0.5 * d$x, -0.2 - 0.4 * d$x)
  P <- exp(E) / rowSums(exp(E))
  d$Y <- t(vapply(seq_len(n), function(i) {
    as.numeric(rmultinom(1, d$N[i], P[i, ]))
  }, numeric(3)))
  colnames(d$Y) <- c("a", "b", "c")
  rs <- list(
    disp_row("negbinomial", grab(frm(y ~ x, family = negbinomial(),
                                     data = d)), "shape"),
    disp_row("negbinomial2", grab(frm(y ~ x,
                                      family = brmsfamily("negbinomial2"),
                                      data = d)), "sigma"),
    disp_row("beta_binomial", grab(frm(yb | trials(N) ~ x,
                                       family = beta_binomial(), data = d)),
             "phi"),
    disp_row("dirichlet_multinomial",
             grab(frm(Y | trials(N) ~ x, family = dirichlet_multinomial(),
                      data = d)), "phi"))
  out[[seed]] <- cbind(seed = seed, do.call(rbind, rs))
}
o <- do.call(rbind, out)
print(o[, 1:5], digits = 5, row.names = FALSE)
cat("\nconditions by family:\n")
print(unique(o[, c("fam", "msgs")]), row.names = FALSE)
