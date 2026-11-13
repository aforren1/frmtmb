# Lane gpstop, fix round: the review's sweep (dev/reviews/gpstop-review-sweep2.R,
# same rows and data) on the fixed build, arms lane and norelin interleaved
# in one process, with every call of gp_relin_regular() logged: the
# descent's code, max |gradient|, its length scales against the positions'
# span, the inner Hessian's smallest eigenvalue, and the verdict. The
# eigenvalue is printed for every candidate, so the cut can be read from
# the numbers. Args as the review's: families, ns, seeds, extra.
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/reviews/gpstop-review-lib.R")
rv_libs("lane")
rv_instrument()
ns_f <- asNamespace("frmtmb")
CUR <- new.env()
CUR$tag <- ""
real_reg <- get("gp_relin_regular", ns_f)
real_ev <- get("gp_inner_ev_min", ns_f)
rv_put("gp_relin_regular", function(obj, frame, o, tol) {
  g <- tryCatch(max(abs(obj$gr(o$par))), error = function(e) NA)
  pl <- obj$env$parList(o$par)
  th <- as.numeric(pl[["theta"]])
  rat <- NULL
  for (bk in frame[["re_blocks"]]) {
    if (!identical(bk[["covstruct"]], "gp")) next
    d2 <- Reduce(`+`, bk[["aux_D2"]])
    d <- sqrt(d2[d2 > 0])
    ls <- exp(th[bk[["theta_idx"]][-1L]])
    rat <- c(rat, min(ls) / min(d), max(ls) / max(d))
  }
  ev <- tryCatch(.GlobalEnv$real_ev(obj, o$par), error = function(e) NA)
  v <- .GlobalEnv$real_reg(obj, frame, o, tol)
  cat(sprintf(paste0("REG %s conv %d gmax %.2e ls/mingap %.3g ls/range %.3g ",
                     "inner_ev %.3e verdict %s\n"),
              .GlobalEnv$CUR$tag, as.integer(o$convergence), g, rat[1],
              rat[2], ev, v))
  v
})
fams <- strsplit(a[1], ",")[[1]]
ns <- as.integer(strsplit(a[2], ",")[[1]])
seeds <- as.integer(strsplit(a[3], ",")[[1]])
extra <- a[4]
bodies <- c("exp_a", "exp_b", "logis_a", "mm_V", "mm_K", "pow_a", "pow_b")
spell <- function(body, k = NULL) {
  g <- if (is.null(k)) quote(gp(x)) else bquote(gp(x, k = .(k)))
  switch(body,
    exp_a = bf(y ~ a * exp(b * z), as.formula(bquote(a ~ .(g))), b ~ 1,
               nl = TRUE),
    exp_b = bf(y ~ a * exp(b * z), a ~ 1, as.formula(bquote(b ~ .(g))),
               nl = TRUE),
    logis_a = bf(y ~ a / (1 + exp(-b * z)), as.formula(bquote(a ~ .(g))),
                 b ~ 1, nl = TRUE),
    mm_V = bf(y ~ V * w / (K + w), as.formula(bquote(V ~ .(g))), K ~ 1,
              nl = TRUE),
    mm_K = bf(y ~ V * w / (K + w), V ~ 1, as.formula(bquote(K ~ .(g))),
              nl = TRUE),
    pow_a = bf(y ~ a * w^b, as.formula(bquote(a ~ .(g))), b ~ 1, nl = TRUE),
    pow_b = bf(y ~ a * w^b, a ~ 1, as.formula(bquote(b ~ .(g))), nl = TRUE))
}
mk <- function(body, fam, n, seed) {
  set.seed(81100 + 1000 * match(fam, c("gaussian", "poisson", "bernoulli",
                                       "student")) +
             100 * seed + n + 10 * match(body, bodies))
  x <- sort(round(runif(n, 0, 10), 3))
  z <- round(rnorm(n), 3)
  w <- round(runif(n, 0.2, 5), 3)
  f <- if (seed %% 2) sin(x) else 2 * sin(3 * x) + x
  fs <- f / max(abs(f))
  eta <- switch(body,
    exp_a = (2 + fs) * exp(0.4 * z),
    exp_b = 2 * exp((0.3 + 0.2 * fs) * z),
    logis_a = (3 + fs) / (1 + exp(-1.5 * z)),
    mm_V = (3 + fs) * w / (1 + w),
    mm_K = 3 * w / (1 + 0.6 * fs + w),
    pow_a = (2 + fs) * w^0.5,
    pow_b = w^(0.5 + 0.3 * fs))
  y <- switch(fam,
    gaussian = eta + 0.1 * rnorm(n),
    student = eta + 0.1 * rt(n, 3),
    poisson = rpois(n, exp(0.6 * eta)),
    bernoulli = rbinom(n, 1, plogis(1.5 * (eta - median(eta)) / sd(eta))))
  data.frame(x, z, w, y)
}
famf <- function(fam) switch(fam, gaussian = gaussian(), poisson = poisson(),
                             student = student(), bernoulli = bernoulli())
rows <- list()
for (fam in fams) for (n in ns) for (seed in seeds) for (body in bodies) {
  rows[[length(rows) + 1L]] <- list(body = body, fam = fam, n = n,
                                    seed = seed, k = NULL)
}
if (extra == "big") {
  for (body in c("exp_a", "exp_b", "pow_a", "mm_K")) {
    rows[[length(rows) + 1L]] <- list(body = body, fam = "gaussian",
                                      n = 400L, seed = 2L, k = NULL)
  }
}
if (extra == "approx") {
  for (body in bodies) rows[[length(rows) + 1L]] <-
    list(body = body, fam = "gaussian", n = 120L, seed = 2L, k = 12L)
}
i <- 0L
for (r in rows) {
  i <- i + 1L
  d <- mk(r$body, r$fam, r$n, r$seed)
  tag <- sprintf("%s_%s_n%d_s%d%s", r$body, r$fam, r$n, r$seed,
                 if (is.null(r$k)) "" else paste0("_k", r$k))
  CUR$tag <- tag
  arms <- if (i %% 2) c("lane", "norelin") else c("norelin", "lane")
  for (arm in arms) {
    rv_run(tag, arm, spell(r$body, r$k), data = d, family = famf(r$fam))
  }
}
cat("DONE\n")
