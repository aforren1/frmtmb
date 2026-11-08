# quadrature = TRUE when a row reads more than one random-effect block.
# One Gauss-Kronrod pass over every effect at once integrated the coarse
# block of (1 | ga/gb) and left the fine one to the Laplace
# approximation: the "quadrature" log-likelihood equaled the Laplace one
# to 1e-8, and summary() failed with "no 'dimnames' attribute for array"
# because sdreport() cannot mix the two methods. 0.73.0 refused every
# such model. Two nested blocks now take the iterated rule of
# R/quadrature-nested.R; crossed blocks, more than two blocks and a
# tape above the memory limit are refused by name.

quad_nest_data <- function(seed = 20261020) {
  set.seed(seed)
  ng <- 12; nt <- 6; n <- ng * nt
  d <- data.frame(ga = factor(rep(seq_len(ng), each = nt)),
                  gb = factor(rep(c("i", "ii"), length.out = n)),
                  tt = factor(rep(seq_len(nt), ng)),
                  x = rnorm(n), w = runif(n, 0.5, 2))
  ab <- interaction(d$ga, d$gb, drop = TRUE)
  d$eta <- 0.2 + 0.5 * d$x + rnorm(ng, 0, 0.7)[d$ga] +
    rnorm(nlevels(ab), 0, 1)[ab]
  d$y <- rbinom(n, 1, plogis(d$eta))
  d$yo <- factor(cut(d$eta + stats::rlogis(n), c(-Inf, -0.8, 0.8, Inf),
                     labels = FALSE), ordered = TRUE)
  d$yg <- d$eta + rnorm(n, 0, 0.5)
  # a by-level that is constant within ga, so a by-split block reads the
  # rows of its own levels only
  d$gby <- factor(c("a", "b", "c")[(as.integer(d$ga) %% 3) + 1])
  d$obs <- factor(seq_len(n))
  d
}

# Gauss-Hermite rule for the standard normal weight, written here so the
# reference shares no code with the package
quad_nest_gh <- function(k) {
  i <- seq_len(k - 1)
  J <- matrix(0, k, k)
  J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, lw = log(e$vectors[1, ]^2))
}
quad_nest_lse <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

# The marginal log-likelihood of a weighted logit model over
# (1 | ga/gb), from the model alone: per ga level, a Gauss-Hermite
# integral over the ga effect of the product of the integrals over its
# ga:gb effects, k nodes per level on the effects' own scale.
quad_nest_ref <- function(d, beta, sd_a, sd_ab, k) {
  g <- quad_nest_gh(k)
  eta0 <- beta[1] + beta[2] * d$x
  tot <- 0
  for (a in levels(d$ga)) {
    ra <- d$ga == a
    lu <- vapply(sd_a * g$x, function(u) {
      s <- 0
      for (b in unique(d$gb[ra])) {
        r <- ra & d$gb == b
        lv <- vapply(sd_ab * g$x, function(v) {
          sum(d$w[r] * stats::dbinom(d$y[r], 1,
                                     stats::plogis(eta0[r] + u + v),
                                     log = TRUE))
        }, 1)
        s <- s + quad_nest_lse(g$lw + lv)
      }
      s
    }, 1)
    tot <- tot + quad_nest_lse(g$lw + lu)
  }
  tot
}

# The standard deviation of each block at outer parameters `par`
quad_nest_sds <- function(fit, par) {
  th <- imp_par_list(fit$frame$par_template, par, fit$frame$map)$theta
  s <- vapply(fit$frame$re_blocks, function(bk) {
    sqrt(cs_entry(bk$covstruct)$vcov(th[bk$theta_idx], bk)[1, 1])
  }, 1)
  stats::setNames(s, vapply(fit$frame$re_blocks, `[[`, "", "group_name"))
}

# The nested rule's log-likelihood at outer parameters `par`, anchored
# there, and the Laplace one, both from the fit's own frame
quad_nest_at <- function(fit, par, nodes) {
  fr <- fit$frame
  lap <- RTMB::MakeADFun(build_objective(fr), fr$par_template,
                         random = "b", map = fr$map, silent = TRUE)
  plan <- quad_check_blocks(fr)
  pl <- imp_par_list(fr$par_template, par, fr$map)
  qo <- quad_nest_objective(fr, plan, quad_nest_anchor(lap, plan, par),
                            nodes, pl)
  qo$prime(pl)
  c(nested = -qo$fn(pl), laplace = -lap$fn(par))
}

quad_nest_fo <- bf(y | weights(w) ~ 1 + x + (1 | ga/gb)) + bernoulli()

test_that("quadrature over (1 | ga/gb) is the iterated integral", {
  d <- quad_nest_data()
  fq <- allow_boundary(frm(quad_nest_fo, data = d, quadrature = TRUE))
  fl <- allow_boundary(frm(quad_nest_fo, data = d))
  # every effect went to the rule, none to Laplace
  expect_length(fq$obj$env$random, 0L)
  s <- quad_nest_sds(fq, fq$opt$par)
  b <- unname(fq$estimates$beta)
  ref <- quad_nest_ref(d, b, s[["ga"]], s[["ga:gb"]], 90)
  # the reference has converged: 60 and 90 nodes agree far inside the
  # Laplace error
  lap_at <- -fl$obj$fn(fq$opt$par)
  err_l <- abs(lap_at - ref)
  expect_lt(abs(quad_nest_ref(d, b, s[["ga"]], s[["ga:gb"]], 60) - ref) /
              err_l, 1e-4)
  err_q <- abs(as.numeric(logLik(fq)) - ref)
  # the conditional centering of the inner rule is what gets this close:
  # 8e-6 against a Laplace error of 0.85 at the default nodes
  expect_lt(err_q / err_l, 5e-5)
  sm <- allow_boundary(summary(fq))
  expect_s3_class(sm, "summary.frmtmb_fit")
  expect_true(all(is.finite(sqrt(diag(vcov(fq))))))
  expect_false(anyNA(fq$estimates$b))
  expect_length(fq$estimates$b, 36L)
})

test_that("the nested rule converges to the reference as the nodes rise", {
  d <- quad_nest_data()
  fq <- allow_boundary(frm(quad_nest_fo, data = d, quadrature = TRUE))
  par <- fq$opt$par
  s <- quad_nest_sds(fq, par)
  ref <- quad_nest_ref(d, unname(fq$estimates$beta), s[["ga"]],
                       s[["ga:gb"]], 90)
  v <- lapply(c(5, 11, 21), function(k) quad_nest_at(fq, par, c(k, k)))
  err <- vapply(v, function(x) abs(x[["nested"]] - ref), 1)
  err_l <- abs(v[[1L]][["laplace"]] - ref)
  expect_true(all(diff(err) < 0))
  expect_lt(err[[3L]] / err_l, 1e-5)
})

test_that("with the inner sd near zero the rule is the one-block rule", {
  d <- quad_nest_data()
  fs <- allow_boundary(frm(bf(y ~ 1 + x + (1 | ga)) + bernoulli(), data = d,
                           quadrature = TRUE))
  fsl <- allow_boundary(frm(bf(y ~ 1 + x + (1 | ga)) + bernoulli(),
                            data = d))
  fn <- allow_boundary(frm(bf(y ~ 1 + x + (1 | ga/gb)) + bernoulli(),
                           data = d, quadrature = TRUE))
  gk <- as.numeric(logLik(fs))
  err_l <- abs(-fsl$obj$fn(fs$opt$par) - gk)
  at <- function(sd_in) {
    p <- fn$opt$par
    p[names(p) == "beta"] <- fs$opt$par[names(fs$opt$par) == "beta"]
    p[names(p) == "theta"] <- c(fs$opt$par[names(fs$opt$par) == "theta"],
                                log(sd_in))
    # theta is the log sd of a dimension-1 block
    expect_equal(quad_nest_sds(fn, p)[["ga:gb"]], sd_in)
    abs(quad_nest_at(fn, p, c(31, 15))[["nested"]] - gk)
  }
  d2 <- at(1e-2)
  d4 <- at(1e-4)
  # the gap closes like the inner variance, to the one-block value
  expect_lt(d4 / d2, 1e-3)
  expect_lt(d4 / err_l, 1e-5)
})

test_that("blocks on two predictors are refused, and the reason named", {
  d <- quad_nest_data()
  # one grouping factor, so the blocks nest in each other; the refusal is
  # about the sigma effect changing the width of the mu integrand
  fo <- bf(yg ~ 1 + x + (1 | ga), sigma ~ (1 | ga)) + gaussian()
  err <- tryCatch(frm(fo, data = d, quadrature = TRUE),
                  error = conditionMessage)
  expect_match(err, "enter the same predictors", fixed = TRUE)
  expect_match(err, "sigma: (1 | ga)", fixed = TRUE)
  expect_match(err, "quadrature = FALSE", fixed = TRUE)
})

test_that("quadrature refuses crossed blocks, and says why", {
  d <- quad_nest_data()
  err <- tryCatch(frm(bf(y ~ 1 + x + (1 | ga) + (1 | tt)) + bernoulli(),
                      data = d, quadrature = TRUE),
                  error = conditionMessage)
  expect_match(err, "are crossed", fixed = TRUE)
  expect_match(err, "does not factor", fixed = TRUE)
  expect_match(err, "1 | ga", fixed = TRUE)
  expect_match(err, "1 | tt", fixed = TRUE)
  expect_match(err, "quadrature = FALSE", fixed = TRUE)
  # the same model fits under Laplace, which the refusal names
  fl <- allow_boundary(frm(bf(y ~ 1 + x + (1 | ga) + (1 | tt)) +
                             bernoulli(), data = d))
  expect_true(is.finite(as.numeric(logLik(fl))))
})

test_that("quadrature refuses three nested blocks, and says why", {
  d <- quad_nest_data()
  d$gc <- factor(rep(c("p", "p", "q", "q", "r", "r"), 12))
  err <- tryCatch(frm(bf(y ~ 1 + x + (1 | ga/gb/gc)) + bernoulli(),
                      data = d, quadrature = TRUE),
                  error = conditionMessage)
  expect_match(err, "has 3 blocks", fixed = TRUE)
  expect_match(err, "a row reads 3 of them", fixed = TRUE)
})

test_that("each refusal of a block layout gives its own reason", {
  d <- quad_nest_data()
  # two responses, each with nested blocks: four blocks, two per row
  err_mv <- tryCatch(
    frm(mvbf(bf(y ~ x + (1 | ga/gb)) + bernoulli(),
             bf(yg ~ x + (1 | ga/gb)) + gaussian()),
        data = d, quadrature = TRUE),
    error = conditionMessage)
  expect_match(err_mv, "for one response", fixed = TRUE)
  # a by-split outer term is one block per by-level
  err_by <- tryCatch(
    frm(bf(y ~ x + (1 | gr(ga, by = gby)) + (1 | ga:gb)) + bernoulli(),
        data = d, quadrature = TRUE),
    error = conditionMessage)
  expect_match(err_by, "is a by-split term", fixed = TRUE)
  # (1 + x || g) is one block of two coefficients, and the message says
  # how to write it as the two blocks quadrature integrates
  err_dd <- tryCatch(frm(bf(y ~ x + (1 + x || ga)) + bernoulli(), data = d,
                         quadrature = TRUE),
                     error = conditionMessage)
  expect_match(err_dd, "(1 | g) + (0 + x | g)", fixed = TRUE)
  fs <- allow_boundary(frm(bf(y ~ x + (1 | ga) + (0 + x | ga)) +
                             bernoulli(), data = d, quadrature = TRUE))
  expect_length(fs$obj$env$random, 0L)
  # the crossed message gives each direction its own clause
  err_x <- tryCatch(frm(bf(y ~ 1 + x + (1 | ga) + (1 | tt)) + bernoulli(),
                        data = d, quadrature = TRUE),
                    error = conditionMessage)
  expect_match(err_x, "so tt does not nest in ga", fixed = TRUE)
  expect_match(err_x, "so ga does not nest in tt", fixed = TRUE)
})

test_that("the memory guard refuses a tape above quad_max_gb", {
  d <- quad_nest_data()
  fq <- allow_boundary(frm(quad_nest_fo, data = d, quadrature = TRUE))
  # the guard's own predictions for this model at two node counts
  fr <- fq$frame
  plan <- quad_check_blocks(fr)
  pred <- function(nodes) {
    quad_nest_guard(fr, plan, nodes, fr$par_template, "b", Inf, NULL)$bytes
  }
  small <- pred(c(3L, 3L))
  big <- pred(c(11L, 15L))
  expect_gt(big, small)
  # a limit between the two: the small rule fits, the default refuses
  lim <- sqrt(small * big) / 2^30
  # three nodes are a crude rule, which may not settle; the guard is the
  # point here
  ok <- allow_warnings(
    allow_boundary(frm(quad_nest_fo, data = d, quadrature = TRUE,
                       control = frmtmb_control(quad_nodes = 3,
                                                quad_max_gb = lim))),
    c("did not settle the nested rule", "Large maximum absolute gradient",
      "did not report convergence"))
  expect_s3_class(ok, "frmtmb_fit")
  msgs <- character(0)
  err <- tryCatch(
    withCallingHandlers(
      allow_boundary(frm(quad_nest_fo, data = d, quadrature = TRUE,
                         control = frmtmb_control(quad_max_gb = lim,
                                                  verbose = TRUE))),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }),
    error = conditionMessage)
  expect_match(err, "above the limit of", fixed = TRUE)
  expect_match(err, "quad_max_gb", fixed = TRUE)
  expect_match(err, "quad_nodes", fixed = TRUE)
  expect_match(err, "quadrature = FALSE", fixed = TRUE)
  # the refusal comes before the Laplace warm start
  expect_false(any(grepl("quadrature warm start", msgs, fixed = TRUE)))
  # without TMB's tape sizes the guard reads a per-row cost, and does not
  # refuse this design at the default limit
  local_mocked_bindings(quad_tape_size = function(obj) NULL)
  nb <- quad_nest_guard(fr, plan, c(11L, 15L), fr$par_template, "b", Inf,
                        NULL)
  expect_null(nb$tape)
  expect_lt(nb$bytes, 4 * 2^30)
})

test_that("a one-row inner block composes with the weights check", {
  d <- quad_nest_data()
  # (1 | obs) on a gaussian response is confounded with sigma, which the
  # olre check says; the point here is the weights check
  ctl <- frmtmb_control(check_olre = "ignore")
  # every row owns an effect: with weights summing above the rows a
  # gaussian likelihood has no maximum, exact integral or not
  d$wbig <- d$w * 2
  expect_error(frm(bf(yg | weights(wbig) ~ 1 + x + (1 | ga) + (1 | obs)) +
                     gaussian(), data = d, quadrature = TRUE,
                   control = ctl),
               "weights")
  # below the rows it fits, and a gaussian integral is exact either way
  d$wsmall <- d$w / 2
  fo <- bf(yg | weights(wsmall) ~ 1 + x + (1 | ga) + (1 | obs)) + gaussian()
  fq <- allow_boundary(frm(fo, data = d, quadrature = TRUE, control = ctl))
  fl <- allow_boundary(frm(fo, data = d, control = ctl))
  expect_length(fq$obj$env$random, 0L)
  ll_q <- as.numeric(logLik(fq))
  ll_l <- as.numeric(logLik(fl))
  expect_lt(abs(ll_q - ll_l) / abs(ll_l), 1e-7)
})

test_that("one block per row is still integrated", {
  d <- quad_nest_data()
  # a by-split gives three blocks over ga, but each row reads one
  fb <- allow_boundary(frm(bf(y ~ x + (1 | gr(ga, by = gby))) +
                             bernoulli(), data = d, quadrature = TRUE))
  expect_length(fb$frame$re_blocks, 3L)
  expect_length(fb$obj$env$random, 0L)
  # a single block: every effect integrated, and summary() works
  fq <- allow_boundary(frm(bf(y | weights(w) ~ x + (1 | ga)) +
                             bernoulli(), data = d, quadrature = TRUE))
  expect_length(fq$obj$env$random, 0L)
  sm <- allow_boundary(summary(fq))
  expect_s3_class(sm, "summary.frmtmb_fit")
})

# Unbalanced schools and classes, as the review of this route built them
# (dev/review-nestquad-acc.R): 16 schools of 1 to 10 classes, classes of
# 1 to 6 rows, given rows per class, uniform draws as there.
quad_nest_mk <- function(seed, n_s, ncl, nr, sd_a, sd_i, gen) {
  set.seed(seed)
  d <- do.call(rbind, lapply(seq_len(n_s), function(s) {
    do.call(rbind, lapply(seq_len(ncl[(s - 1) %% length(ncl) + 1]),
                          function(c) {
      k <- nr[sample.int(length(nr), 1)]
      data.frame(school = s, cls = c, x = rnorm(k))
    }))
  }))
  d$school <- factor(d$school)
  d$cls_uni <- factor(paste0(d$school, "_", d$cls))
  ua <- rnorm(nlevels(d$school), 0, sd_a)
  ui <- rnorm(nlevels(d$cls_uni), 0, sd_i)
  d$ua <- ua[d$school]
  d$ui <- ui[d$cls_uni]
  gen(d)
}

# log int exp(f(v)) dv by Gauss-Hermite centered on the mode of f and
# scaled by its curvature, both found numerically
quad_nest_adapt <- function(f, k, lo, hi) {
  m <- stats::optimize(f, c(lo, hi), maximum = TRUE, tol = 1e-10)$maximum
  h <- 1e-4
  c2 <- -(f(m + h) - 2 * f(m) + f(m - h)) / h^2
  if (!is.finite(c2) || c2 <= 0) c2 <- 1
  sc <- 1 / sqrt(c2)
  g <- quad_nest_gh(k)
  v <- m + sc * g$x
  log(sc) + quad_nest_lse(g$lw + g$x^2 / 2 + 0.5 * log(2 * pi) +
                            vapply(v, f, 1))
}

test_that("the inner rule's conditional centering is what makes it exact", {
  d <- quad_nest_mk(20261106, 16, c(1, 10, 2, 1, 5, 3), c(1, 1, 2, 4, 6),
                    0.8, 1.0, function(d) {
    d$off <- stats::runif(nrow(d), -0.5, 0.5)
    d$w <- stats::runif(nrow(d), 0.5, 2)
    d$y <- stats::rbinom(nrow(d), 1, stats::plogis(0.2 + 0.6 * d$x + d$off +
                                                     d$ua + d$ui))
    d
  })
  fo <- bf(y | weights(w) ~ x + offset(off) + (1 | school/cls_uni)) +
    bernoulli()
  fq <- allow_boundary(frm(fo, data = d, quadrature = TRUE))
  fl <- allow_boundary(frm(fo, data = d))
  s <- quad_nest_sds(fq, fq$opt$par)
  b <- unname(fq$estimates$beta)
  eta0 <- b[1] + b[2] * d$x + d$off
  lp <- function(r, eta) {
    sum(d$w[r] * stats::dbinom(d$y[r], 1, stats::plogis(eta), log = TRUE))
  }
  ref_k <- function(k) {
    tot <- 0
    for (a in levels(d$school)) {
      ra <- which(d$school == a)
      fa <- function(u) {
        out <- stats::dnorm(u, 0, s[[1L]], log = TRUE)
        for (cc in unique(d$cls_uni[ra])) {
          rc <- ra[d$cls_uni[ra] == cc]
          out <- out + quad_nest_adapt(function(v) {
            stats::dnorm(v, 0, s[[2L]], log = TRUE) + lp(rc, eta0[rc] + u + v)
          }, k, -20, 20)
        }
        out
      }
      tot <- tot + quad_nest_adapt(fa, k, -20, 20)
    }
    tot
  }
  ref <- ref_k(40)
  err_l <- abs(-fl$obj$fn(fq$opt$par) - ref)
  expect_lt(abs(ref_k(30) - ref) / err_l, 1e-9)
  # 2.4e-8 against a Laplace error of 0.85; an inner rule centered on the
  # inner block's marginal instead gave 9.8e-7 (dev/nestquad-findings.md)
  expect_lt(abs(as.numeric(logLik(fq)) - ref) / err_l, 2e-7)
})

# A near-separable Bernoulli, the review's design: few rows per class and
# a steep slope, where the frozen rule's optimum lies far from its anchor
quad_nest_sep_data <- function() {
  quad_nest_mk(20261106, 10, c(2, 3), c(2, 3), 0.5, 0.5, function(d) {
    d$w <- 1
    d$y <- as.integer(4 * d$x + d$ua + d$ui +
                        stats::rlogis(nrow(d), scale = 0.3) > 0)
    d
  })
}

test_that("the re-anchoring loop does not cycle, and says when it stops", {
  d <- quad_nest_sep_data()
  fo <- bf(y ~ x + (1 | school/cls_uni)) + bernoulli()
  seen <- character(0)
  msgs <- character(0)
  fq <- withCallingHandlers(
    frm(fo, data = d, quadrature = TRUE),
    warning = function(w) {
      seen <<- c(seen, class(w)[1L])
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  # either the rule settled, or a warning of its own class says it did
  # not; a stop at the round cap with nothing said is the defect
  expect_true(isTRUE(fq$quad_nested$settled) ||
                "frmtmb_quad_unsettled" %in% seen)
  # every other warning here is the generic convergence one
  expect_true(all(seen %in% c("frmtmb_quad_unsettled", "frmtmb_warning")))
  if (!isTRUE(fq$quad_nested$settled)) {
    # the warning gives the rule's own error at the estimate, against a
    # larger rule, and says more nodes are not a sure remedy
    um <- msgs[seen == "frmtmb_quad_unsettled"]
    expect_match(um, "the rule's own error there is of the order of",
                 fixed = TRUE)
    expect_match(um, "against 21 x 21", fixed = TRUE)
    expect_match(um, "not a sure remedy", fixed = TRUE)
    re <- fq$quad_nested$rule_error
    expect_true(is.finite(re) && re > 0)
    expect_match(um, format(re, digits = 2), fixed = TRUE)
  }
  ll <- as.numeric(logLik(fq))
  # a probability is at most 1; the frozen rule of the first build
  # reported +29.6 in a round
  expect_lt(ll, 0)
  # the value is the rule anchored at the estimate, and no worse than
  # the rule anchored at the Laplace optimum it started from
  fl <- allow_boundary(frm(fo, data = d))
  at <- quad_nest_at(fq, fq$opt$par, c(11, 15))[["nested"]]
  expect_lt(abs(ll - at) / abs(ll), 1e-8)
  start <- quad_nest_at(fq, fl$opt$par, c(11, 15))[["nested"]]
  expect_gte(ll, start)
})

test_that("a log-likelihood above the family's bound is not reported", {
  d <- quad_nest_data()
  fq <- allow_boundary(frm(quad_nest_fo, data = d, quadrature = TRUE))
  fr <- fq$frame
  p <- fq$opt$par
  tpl <- fr$par_template
  # a Bernoulli log-likelihood of +1 is impossible
  expect_error(quad_nest_check_bound(fr, -1, NULL, tpl, p, fr$map),
               class = "frmtmb_quad_bound")
  expect_silent(quad_nest_check_bound(fr, 1, NULL, tpl, p, fr$map))
  # a density can exceed 1, so a gaussian response has no such bound
  fg <- allow_boundary(frm(bf(yg ~ 1 + x + (1 | ga/gb)) + gaussian(),
                           data = d, quadrature = TRUE))
  expect_silent(quad_nest_check_bound(fg$frame, -1, NULL,
                                      fg$frame$par_template, fg$opt$par,
                                      fg$frame$map))
})

test_that("profile intervals and condVar are refused on a nested fit", {
  d <- quad_nest_data()
  fq <- allow_boundary(frm(quad_nest_fo, data = d, quadrature = TRUE))
  cls <- "frmtmb_quad_nested_unsupported"
  expect_error(confint(fq, method = "profile", parm = "x"), class = cls)
  expect_error(confint(fq, method = "uniroot", parm = "x"), class = cls)
  expect_error(profile(fq, "x"), class = cls)
  # with no parameter named, the refusal still comes first
  expect_error(profile(fq), class = cls)
  expect_error(ranef(fq, condVar = TRUE), class = cls)
  # Wald intervals and the modes work
  w <- confint(fq)
  expect_true(all(is.finite(w)))
  expect_false(anyNA(unlist(ranef(fq))))
  # a one-block quadrature fit keeps its profile
  f1 <- allow_boundary(frm(bf(y ~ 1 + x + (1 | ga)) + bernoulli(), data = d,
                           quadrature = TRUE))
  pr <- confint(f1, method = "profile", parm = "x")
  expect_true(all(is.finite(pr)))
})
