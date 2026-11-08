# The memory refusal on an exgaussian nested fit: with and without verbose.
source("C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev/review-nestquad-common.R")
set.seed(20261108)
n <- 3000; n_out <- n %/% 20
d <- data.frame(school = factor(rep(seq_len(n_out), each = 20)),
                cls = rep(rep(1:4, each = 5), n_out), x = rnorm(n))
d$cls_uni <- factor(paste0(d$school, "_", d$cls))
eta <- 0.3 + 0.4 * d$x + rnorm(n_out, 0, 0.5)[d$school] +
  rnorm(nlevels(d$cls_uni), 0, 0.5)[d$cls_uni]
d$y <- eta + rnorm(n) + rexp(n, 2)
fo <- bf(y ~ x + (1 | school/cls_uni)) + exgaussian()
for (v in c(0, 1)) {
  t0 <- proc.time()[[3]]
  r <- tryCatch(frm(fo, data = d, quadrature = TRUE,
                    control = frmtmb_control(verbose = v)),
                error = function(e) e)
  cat(sprintf("verbose %d: %.1fs class %s\n  %s\n", v, proc.time()[[3]] - t0,
              paste(class(r), collapse = "/"),
              if (inherits(r, "error")) conditionMessage(r) else "fit"))
  flush(stdout())
}
cat("DONE\n")
