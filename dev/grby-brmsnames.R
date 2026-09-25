# brms 2.23.0's names for by-split terms, read off its own rename_re()
# on brm(empty = TRUE) objects: nothing is compiled or sampled.
# Run: Rscript dev/grby-brmsnames.R > dev/grby-log/brmsnames.txt
.libPaths(c("/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages(library(brms))
set.seed(1)
ng <- 12
d <- data.frame(g = factor(rep(1:ng, each = 5)))
d$f <- factor(ifelse(as.integer(d$g) <= 6, "a", "b"))
d$x <- rnorm(nrow(d)); d$y <- rnorm(nrow(d))
d$g1 <- d$g; d$g2 <- factor(sample(1:ng, nrow(d), TRUE), levels = 1:ng)
d$f1 <- d$f; d$f2 <- factor(ifelse(as.integer(d$g2) <= 6, "a", "b"))
show <- function(form, pars) {
  fit <- brm(form, data = d, empty = TRUE)
  bf <- brms:::brmsframe(brms:::brmsterms(fit$formula), data = fit$data)
  out <- brms:::rename_re(bf, pars = pars)
  for (o in out) if (!is.null(o$fnames)) print(o$fnames) else print(o)
}
show(y ~ x + (1 + x | gr(g, by = f)),
     c("sd_1[1,1]", "sd_1[2,1]", "sd_1[1,2]", "sd_1[2,2]", "cor_1_1[1]",
       "cor_1_2[1]"))
show(y ~ x + (1 | mm(g1, g2, by = cbind(f1, f2))), c("sd_1[1,1]", "sd_1[1,2]"))
show(bf(y ~ x + (1 | gr(g, by = f)), sigma ~ (1 | gr(g, by = f))),
     c("sd_1[1,1]", "sd_1[1,2]", "sd_2[1,1]", "sd_2[1,2]"))
