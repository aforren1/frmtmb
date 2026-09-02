# Probe: two candidate routes to per-CLUSTER scores of the marginal
# (Laplace) objective. Both rest on the same identity: when every random
# effect is nested in the clustering factor, nll = sum_g nll_g, and
# setting cluster g's data weight to 1 and every other cluster's to 0
# leaves an objective that IS nll_g -- the data-free clusters integrate
# their Gaussian prior to exactly 1, contributing 0 with zero gradient.
#
#   route A  re-tape: bake the 0/1 mask into aterm_values$weights and
#            call MakeADFun once per cluster.  No objective.R change.
#   route B  weight parameter: one tape carrying an extra length-G
#            parameter `clw`; obj$gr(c(theta, e_g)) reads the score off
#            the theta block.  Needs the (gated) hook in build_objective.
#
# Run: Rscript dev/sandwich/probe-scores.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))

sink(file.path("dev", "sandwich", "probe-scores.txt"), split = TRUE)

mask_frame <- function(frame, keep) {
  for (r in names(frame$aterm_values)) {
    w0 <- frame$aterm_values[[r]]$weights %||% rep(1, frame$n_obs)
    frame$aterm_values[[r]]$weights <- w0 * keep
  }
  frame
}

random_names <- function(fit) {
  ri <- fit$obj$env$random
  if (!length(ri)) return(NULL)
  unique(names(fit$obj$env$par)[ri])
}

route_A <- function(fit, cl) {
  gs <- levels(cl)
  tpl <- fit$obj$env$parameters
  random <- random_names(fit)
  S <- matrix(NA_real_, length(gs), length(fit$opt$par))
  for (i in seq_along(gs)) {
    fr <- mask_frame(fit$frame, as.numeric(cl == gs[i]))
    o <- RTMB::MakeADFun(build_objective(fr), tpl, random = random,
                         map = fit$frame$map, silent = TRUE)
    S[i, ] <- o$gr(fit$opt$par)
  }
  colnames(S) <- names(fit$opt$par)
  S
}

route_B <- function(fit, cl) {
  gs <- levels(cl)
  fr <- fit$frame
  fr$cluster_w <- lapply(stats::setNames(nm = names(fr$aterm_values)),
                         function(r) as.integer(cl))
  tpl <- fit$obj$env$parameters
  tpl$clw <- rep(1, length(gs))
  random <- random_names(fit)
  o <- RTMB::MakeADFun(build_objective(fr), tpl, random = random,
                       map = fit$frame$map, silent = TRUE)
  keep <- names(o$par) != "clw"
  stopifnot(identical(names(o$par)[keep], names(fit$opt$par)))
  S <- matrix(NA_real_, length(gs), sum(keep))
  p <- o$par
  p[keep] <- fit$opt$par
  for (i in seq_along(gs)) {
    p[!keep] <- as.numeric(seq_along(gs) == i)
    S[i, ] <- o$gr(p)[keep]
  }
  colnames(S) <- names(fit$opt$par)
  S
}

bench <- function(lab, f) {
  t0 <- proc.time()[["elapsed"]]
  v <- f()
  cat(sprintf("%-28s %7.3f s\n", lab, proc.time()[["elapsed"]] - t0))
  v
}

report <- function(tag, fit, cl) {
  cat("\n== ", tag, " (n = ", fit$frame$n_obs, ", G = ",
      nlevels(cl), ", p = ", length(fit$opt$par), ") ==\n", sep = "")
  A <- bench("route A (re-tape per cluster)", function() route_A(fit, cl))
  B <- bench("route B (weight parameter)", function() route_B(fit, cl))
  g <- fit$obj$gr(fit$opt$par)
  cat(sprintf("max |A - B|                  %.3e\n", max(abs(A - B))))
  cat(sprintf("max |colSums(A) - obj$gr|    %.3e\n",
              max(abs(colSums(A) - drop(g)))))
  cat(sprintf("max |obj$gr| at optimum      %.3e\n", max(abs(g))))
  invisible(list(A = A, B = B))
}

## 1. GLM, one cluster per row -- the score machinery must reduce to the
##    classical per-observation scores there.
set.seed(1)
n <- 120
d1 <- data.frame(x = rnorm(n), z = runif(n))
d1$y <- rpois(n, exp(0.4 + 0.6 * d1$x - 0.3 * d1$z))
f1 <- frm(bf(y ~ x + z) + poisson(), data = d1)
r1 <- report("poisson glm, cluster = row", f1, factor(seq_len(n)))
g1 <- glm(y ~ x + z, poisson(), d1)
cat(sprintf("max |scores - sandwich::estfun(glm)| %.3e\n",
            max(abs(-r1$A - sandwich::estfun(g1)))))

## 2. Gaussian LMM, cluster = the grouping factor.
set.seed(2)
G <- 30; m <- 8
d2 <- data.frame(g = factor(rep(seq_len(G), each = m)),
                 x = rnorm(G * m))
d2$y <- 1 + 0.5 * d2$x + rnorm(G, 0, 0.9)[d2$g] + rnorm(G * m)
f2 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d2, REML = FALSE)
report("gaussian lmm, cluster = g", f2, d2$g)

## 3. Bigger G, to see how each route scales.
set.seed(3)
G <- 100; m <- 6
d3 <- data.frame(g = factor(rep(seq_len(G), each = m)),
                 x = rnorm(G * m))
d3$y <- rbinom(G * m, 1, plogis(-0.2 + 0.8 * d3$x +
                                 rnorm(G, 0, 1)[d3$g]))
f3 <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = d3)
report("bernoulli glmm, G = 100", f3, d3$g)

## 4. Cluster COARSER than the random effect (nesting, not equality).
set.seed(4)
d4 <- data.frame(school = factor(rep(1:20, each = 12)))
d4$class <- factor(paste0(d4$school, ".", rep(rep(1:3, each = 4), 20)))
d4$x <- rnorm(240)
d4$y <- 1 + 0.4 * d4$x + rnorm(20, 0, .7)[d4$school] +
  rnorm(60, 0, .5)[d4$class] + rnorm(240)
f4 <- frm(bf(y ~ x + (1 | school) + (1 | class)) + gaussian(),
          data = d4, REML = FALSE)
report("nested re, cluster = school", f4, d4$school)

sink()
