# Prototype check of nested quadrature against the brute-force reference.
#   Rscript dev/nestquad-proto.R
wt <- "C:/Users/adf44/source/r/frmtmb-wt-nestquad"
.libPaths(c("C:/Users/adf44/source/r/wt-nestquad-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb)); cs_entry <- frmtmb:::cs_entry
source(file.path(wt, "dev/nestquad-ref.R"))

set.seed(20261020)
ng <- 12; nt <- 6; n <- ng * nt
d <- data.frame(ga = factor(rep(seq_len(ng), each = nt)),
                gb = factor(rep(c("i", "ii"), length.out = n)),
                x = rnorm(n), w = runif(n, 0.5, 2))
ab <- interaction(d$ga, d$gb, drop = TRUE)
d$eta <- 0.2 + 0.5 * d$x + rnorm(ng, 0, 0.7)[d$ga] +
  rnorm(nlevels(ab), 0, 1)[ab]
d$y <- rbinom(n, 1, plogis(d$eta))

fo <- bf(y ~ 1 + x + (1 | ga/gb)) + bernoulli()
t0 <- proc.time()[[3]]
fq <- frm(fo, data = d, quadrature = TRUE, verbose = TRUE)
cat("fit time", proc.time()[[3]] - t0, "\n")
fl <- frm(fo, data = d)
sds <- vapply(fq$frame$re_blocks, function(bk) {
  th <- fq$estimates$theta[bk$theta_idx]
  sqrt(cs_entry(bk$covstruct)$vcov(th, bk)[1, 1])
}, 1)
names(sds) <- vapply(fq$frame$re_blocks, `[[`, "", "group_name")
print(sds)
beta <- unname(fq$estimates$beta)
eta0 <- beta[1] + beta[2] * d$x
ref <- if (FALSE) ref_nested_integrate("bernoulli", d$y, eta0, d$ga, ab,
                            sds[["ga"]], sds[["ga:gb"]]) else NA
ref_gh80 <- ref_nested_gh("bernoulli", d$y, eta0, d$ga, ab, sds[["ga"]],
                        sds[["ga:gb"]], k = 80)
lap_at <- -fl$obj$fn(fq$opt$par)
cat(sprintf("quad %.10f ref_int %.10f ref_gh80 %.10f lap_at %.10f\n",
            as.numeric(logLik(fq)), ref, ref_gh80, lap_at))
cat(sprintf("err_q %.3e err_l %.3e ratio %.3e\n",
            abs(as.numeric(logLik(fq)) - ref_gh80), abs(lap_at - ref_gh80),
            abs(as.numeric(logLik(fq)) - ref_gh80) / abs(lap_at - ref_gh80)))
print(summary(fq))
