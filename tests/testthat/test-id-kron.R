# Per-component Kronecker blocks for |ID|-merged groups (v0.32).
#
# Terms sharing an |ID| key merge into ONE covariance block. Until v0.32
# that block was always unstructured, so a gr(cov = A) / gr(prec = Q)
# term whose key was shared lost its matrix; v0.29 refused the construct
# rather than fit it wrong. It is now built as a gr_cov / gr_prec block
# of the total merged dimension, which is the SAME joint density as the
# long-format spelling with one gr() term of that dimension.
#
# The equivalence below is the load-bearing check: two spellings of one
# model, fitted independently, agreeing on the log likelihood, the
# variance components and the fixed effects.

# The pedigree generator of vignette("case-studies") section 1b, under
# names of its own so the two test files can share a session.
kron_pedigree <- function(nsire, ndam_per, noff) {
  nfound <- nsire + nsire * ndam_per
  n <- nfound + nsire * ndam_per * noff
  ped <- data.frame(id = seq_len(n), sire = NA_integer_, dam = NA_integer_)
  k <- nfound
  for (s in seq_len(nsire)) {
    for (j in seq_len(ndam_per)) {
      dam <- nsire + ndam_per * (s - 1L) + j
      for (o in seq_len(noff)) {
        k <- k + 1L
        ped$sire[k] <- s
        ped$dam[k] <- dam
      }
    }
  }
  ped
}

kron_relmat <- function(ped) {
  n <- nrow(ped)
  A <- diag(n)
  for (i in seq_len(n)) {
    s <- ped$sire[i]
    d <- ped$dam[i]
    if (!is.na(s)) {
      A[i, i] <- 1 + 0.5 * A[s, d]
      for (j in seq_len(i - 1)) {
        A[i, j] <- A[j, i] <- 0.5 * (A[j, s] + A[j, d])
      }
    }
  }
  dimnames(A) <- list(as.character(ped$id), as.character(ped$id))
  A
}

#' Two correlated traits over a pedigree, in both spellings: `wide` has
#' one row per animal with y1 and y2 side by side (mvbf + |ID|), `long`
#' stacks them with a trait factor (one gr() term).
kron_traits <- function(seed = 5, nsire = 8, ndam_per = 2, noff = 4,
                        G = matrix(c(1.0, 0.5, 0.5, 0.8), 2, 2),
                        s1 = 0.7, s2 = 0.9) {
  set.seed(seed)
  ped <- kron_pedigree(nsire, ndam_per, noff)
  A <- kron_relmat(ped)
  n <- nrow(ped)
  U <- t(chol(A)) %*% matrix(stats::rnorm(n * 2), n, 2) %*% chol(G)
  ids <- factor(as.character(ped$id), levels = as.character(ped$id))
  y1 <- 3 + U[, 1] + stats::rnorm(n, 0, s1)
  y2 <- 1 + U[, 2] + stats::rnorm(n, 0, s2)
  list(
    A = A,
    wide = data.frame(id = ids, y1 = y1, y2 = y2),
    long = data.frame(
      id = factor(rep(as.character(ped$id), times = 2),
                  levels = as.character(ped$id)),
      trait = factor(rep(c("y1", "y2"), each = n)),
      value = c(y1, y2))
  )
}

kron_fit_wide <- function(d, ...) {
  frm(mvbf(bf(y1 ~ 1 + (1 | q | gr(id, cov = A))) + gaussian(),
           bf(y2 ~ 1 + (1 | q | gr(id, cov = A))) + gaussian()),
      data = d$wide, data2 = list(A = d$A), ...)
}

kron_fit_long <- function(d, ...) {
  frm(bf(value ~ 0 + trait + (0 + trait | gr(id, cov = A)),
         sigma ~ 0 + trait) + gaussian(),
      data = d$long, data2 = list(A = d$A), ...)
}


## the equivalence ------------------------------------------------------

test_that("|ID|-merged gr(cov =) is the long-format model", {
  skip_on_cran()
  d <- kron_traits()
  fw <- kron_fit_wide(d)
  fl <- kron_fit_long(d)

  bk <- fw$frame$re_blocks[[1]]
  expect_length(fw$frame$re_blocks, 1L)
  expect_equal(bk$covstruct, "gr_cov")
  expect_equal(bk$dim, 2L)
  expect_equal(bk$n_levels, nrow(d$A))
  # the index maps are rebuilt at the MERGED dimension, not at the
  # dimension of the first component (which is 1 here)
  expect_equal(length(bk$aux_kron$ia), (2L * nrow(d$A))^2)

  # the same joint density, so the same maximum and the same df
  expect_equal(as.numeric(logLik(fw)), as.numeric(logLik(fl)),
               tolerance = 1e-8)
  expect_equal(attr(logLik(fw), "df"), attr(logLik(fl), "df"))

  # the same genetic covariance, up to optimizer tolerance
  Vw <- unname(VarCorr(fw)[[1]])
  Vl <- unname(VarCorr(fl)[[1]])
  expect_equal(Vw, Vl, tolerance = 1e-4)

  # the same fixed effects: two trait means and two residual sds. The
  # long format keeps them in one linear predictor each, the wide one
  # splits them per response, so they are compared by value.
  expect_equal(unname(c(fixef(fw)$y1_mu, fixef(fw)$y2_mu)),
               unname(fixef(fl)$mu), tolerance = 1e-4)
  expect_equal(unname(c(fixef(fw)$y1_sigma, fixef(fw)$y2_sigma)),
               unname(fixef(fl)$sigma), tolerance = 1e-4)
})

test_that("the merged block actually reads the relationship matrix", {
  skip_on_cran()
  d <- kron_traits()
  fw <- kron_fit_wide(d)
  # the v0.29 wrong answer was a us block, which is exactly cov = I.
  # Refitting with the identity has to move the log likelihood.
  I <- diag(nrow(d$A))
  dimnames(I) <- dimnames(d$A)
  f_I <- frm(mvbf(bf(y1 ~ 1 + (1 | q | gr(id, cov = I))) + gaussian(),
                  bf(y2 ~ 1 + (1 | q | gr(id, cov = I))) + gaussian()),
             data = d$wide, data2 = list(I = I))
  expect_gt(as.numeric(logLik(fw)) - as.numeric(logLik(f_I)), 1)

  # and a plain us |ID| merge is that identity fit
  f_us <- frm(mvbf(bf(y1 ~ 1 + (1 | q | id)) + gaussian(),
                   bf(y2 ~ 1 + (1 | q | id)) + gaussian()),
              data = d$wide)
  expect_equal(f_us$frame$re_blocks[[1]]$covstruct, "us")
  expect_equal(as.numeric(logLik(f_us)), as.numeric(logLik(f_I)),
               tolerance = 1e-5)
})

test_that("the merged Kronecker equivalence holds under REML", {
  skip_on_cran()
  d <- kron_traits()
  fw <- kron_fit_wide(d, REML = TRUE)
  fl <- kron_fit_long(d, REML = TRUE)
  expect_equal(as.numeric(logLik(fw)), as.numeric(logLik(fl)),
               tolerance = 1e-6)
  expect_equal(unname(VarCorr(fw)[[1]]), unname(VarCorr(fl)[[1]]),
               tolerance = 1e-4)
})

test_that("|ID|-merged gr(prec =) matches the cov = spelling", {
  skip_on_cran()
  d <- kron_traits()
  Q <- solve(d$A)
  dimnames(Q) <- dimnames(d$A)
  fq <- frm(mvbf(bf(y1 ~ 1 + (1 | q | gr(id, prec = Q))) + gaussian(),
                 bf(y2 ~ 1 + (1 | q | gr(id, prec = Q))) + gaussian()),
            data = d$wide, data2 = list(Q = Q))
  bk <- fq$frame$re_blocks[[1]]
  expect_equal(bk$covstruct, "gr_prec")
  expect_equal(bk$dim, 2L)
  # one sparse piece per (a, b) entry of the 2 x 2 within-level
  # precision, built at the merged dimension
  expect_length(bk$aux_Qk, 3L)

  fw <- kron_fit_wide(d)
  expect_equal(as.numeric(logLik(fq)), as.numeric(logLik(fw)),
               tolerance = 1e-6)
  expect_equal(unname(VarCorr(fq)[[1]]), unname(VarCorr(fw)[[1]]),
               tolerance = 1e-4)
})

test_that("merging across dpars of one response takes the same path", {
  skip_on_cran()
  d <- kron_traits()
  # the sigma component is weakly identified on one row per animal, so
  # the merged correlation runs to the boundary and the optimizer says
  # so; the block structure is what this test is about
  f <- suppressWarnings(
    frm(bf(y1 ~ 1 + (1 | q | gr(id, cov = A)),
           sigma ~ 1 + (1 | q | gr(id, cov = A))) + gaussian(),
        data = d$wide, data2 = list(A = d$A),
        control = frmtmb_control(check_olre = "ignore")))
  bk <- f$frame$re_blocks[[1]]
  expect_length(f$frame$re_blocks, 1L)
  expect_equal(bk$covstruct, "gr_cov")
  expect_equal(bk$dim, 2L)
  expect_equal(bk$cnms, c("y1.mu:(Intercept)", "y1.sigma:(Intercept)"))
  # draw_b() reaches the merged block through simulate()
  set.seed(3)
  s <- simulate(f, nsim = 2)
  expect_equal(dim(as.matrix(s)), c(nrow(d$wide), 2L))
  expect_true(all(is.finite(as.matrix(s))))
})


## data2, persistence ---------------------------------------------------

test_that("the merged matrix rides data2 through saveRDS and refit", {
  skip_on_cran()
  d <- kron_traits()
  fw <- kron_fit_wide(d)
  f <- tempfile(fileext = ".rds")
  saveRDS(fw, f)
  # a fresh environment: nothing named A is reachable except through
  # the data2 stored on the fit
  reloaded <- local({
    rm(list = ls())
    readRDS(f)
  })
  expect_equal(as.numeric(logLik(reloaded)), as.numeric(logLik(fw)))
  expect_equal(unname(VarCorr(reloaded)[[1]]), unname(VarCorr(fw)[[1]]))
  again <- update(reloaded, data = d$wide)
  expect_equal(as.numeric(logLik(again)), as.numeric(logLik(fw)),
               tolerance = 1e-6)
  expect_equal(again$frame$re_blocks[[1]]$covstruct, "gr_cov")
})

test_that("linked gr() terms resolving to different matrices are refused", {
  skip_on_cran()
  d <- kron_traits(nsire = 4, ndam_per = 2, noff = 3)
  # the merge key is the |ID| label plus the deparsed grouping call, so
  # both formulas write gr(id, cov = A). Give each its own formula
  # environment, binding A to a different matrix: the check that catches
  # this compares the RESOLVED matrices in frame assembly, which is why
  # deparsed sameness is not enough on its own.
  I <- diag(nrow(d$A))
  dimnames(I) <- dimnames(d$A)
  e1 <- new.env(parent = globalenv()); e1$A <- d$A
  e2 <- new.env(parent = globalenv()); e2$A <- I
  f1 <- eval(quote(bf(y1 ~ 1 + (1 | q | gr(id, cov = A))) + gaussian()), e1)
  f2 <- eval(quote(bf(y2 ~ 1 + (1 | q | gr(id, cov = A))) + gaussian()), e2)
  expect_error(frm(mvbf(f1, f2), data = d$wide),
               "must resolve to the same matrix")
})


## what stays refused ---------------------------------------------------

test_that("mixed structures under one |ID| label stay refused", {
  skip_on_cran()
  d <- kron_traits(nsire = 4, ndam_per = 2, noff = 3)
  A <- d$A
  Q <- solve(A)
  dimnames(Q) <- dimnames(A)
  B <- diag(nrow(A))
  dimnames(B) <- dimnames(A)

  # us against gr_cov
  expect_error(
    frm(mvbf(bf(y1 ~ 1 + (1 | q | id)) + gaussian(),
             bf(y2 ~ 1 + (1 | q | gr(id, cov = A))) + gaussian()),
        data = d$wide, data2 = list(A = A)),
    "more than one grouping specification")
  # cov against prec
  expect_error(
    frm(mvbf(bf(y1 ~ 1 + (1 | q | gr(id, cov = A))) + gaussian(),
             bf(y2 ~ 1 + (1 | q | gr(id, prec = Q))) + gaussian()),
        data = d$wide, data2 = list(A = A, Q = Q)),
    "more than one grouping specification")
  # two different matrices
  expect_error(
    frm(mvbf(bf(y1 ~ 1 + (1 | q | gr(id, cov = A))) + gaussian(),
             bf(y2 ~ 1 + (1 | q | gr(id, cov = B))) + gaussian()),
        data = d$wide, data2 = list(A = A, B = B)),
    "more than one grouping specification")
  # an explicit covariance structure on the bar is still refused at
  # parse time, before the gr() rewrite can reach it
  expect_error(
    frm(mvbf(bf(y1 ~ 1 + ar1(0 + o | q | id)) + gaussian(),
             bf(y2 ~ 1 + (1 | q | id)) + gaussian()),
        data = d$wide),
    "only supported for default")
})

test_that("an unshared |ID| key on a gr() term is still a no-op", {
  skip_on_cran()
  d <- kron_traits(nsire = 4, ndam_per = 2, noff = 3)
  f_id <- frm(bf(y1 ~ 1 + (1 | q | gr(id, cov = A))) + gaussian(),
              data = d$wide, data2 = list(A = d$A))
  f_plain <- frm(bf(y1 ~ 1 + (1 | gr(id, cov = A))) + gaussian(),
                 data = d$wide, data2 = list(A = d$A))
  expect_equal(f_id$frame$re_blocks[[1]]$covstruct, "gr_cov")
  expect_equal(f_id$frame$re_blocks[[1]]$dim, 1L)
  expect_null(f_id$frame$re_blocks[[1]]$aux_kron)
  expect_equal(as.numeric(logLik(f_id)), as.numeric(logLik(f_plain)),
               tolerance = 1e-10)
})


## methods on the merged block ------------------------------------------

test_that("the merged block reports sane named summaries", {
  skip_on_cran()
  d <- kron_traits()
  fw <- kron_fit_wide(d)
  nm <- c("y1.mu:(Intercept)", "y2.mu:(Intercept)")

  V <- VarCorr(fw)[[1]]
  expect_equal(dimnames(V), list(nm, nm))
  expect_true(all(diag(V) > 0))
  expect_lt(abs(stats::cov2cor(V)[1, 2]), 1)
  expect_output(print(VarCorr(fw)), "y1.mu")

  re <- ranef(fw, condVar = TRUE)[[1]]
  expect_equal(dim(re), c(nrow(d$A), 2L))
  expect_equal(colnames(re), nm)
  expect_equal(rownames(re), rownames(d$A))
  expect_true(all(is.finite(attr(re, "condSD"))))

  cv <- confint_varcorr(fw)
  # two sds and one correlation, all from the one merged block
  expect_equal(nrow(cv), 3L)
  expect_equal(cv$type, c("sd", "sd", "cor"))
  expect_equal(length(unique(cv$block)), 1L)
  expect_true(all(cv$lwr < cv$estimate & cv$estimate < cv$upr))

  # hypothesis() sees the merged coefficients under the
  # correlated-slopes naming; the correlation is the genetic one
  vn <- variables(fw)
  expect_true(all(c("sd_id__y1.muIntercept", "sd_id__y2.muIntercept",
                    "cor_id__y1.muIntercept__y2.muIntercept") %in% vn))
  h <- hypothesis(fw, "cor_id__y1.muIntercept__y2.muIntercept = 0")
  expect_equal(h$estimate, stats::cov2cor(V)[1, 2], tolerance = 1e-6)
  expect_true(is.finite(h$se))
})

test_that("predict() on the merged block behaves as for a single gr()", {
  skip_on_cran()
  d <- kron_traits(nsire = 4, ndam_per = 2, noff = 3)
  fw <- kron_fit_wide(d)
  nd <- d$wide[1:4, ]
  p <- predict(fw, newdata = nd, resp = "y1")
  expect_length(p, 4L)
  expect_true(all(is.finite(p)))
  # the levels ARE the structure, so an unseen one gets the population
  # value and no block contribution of its own
  nd2 <- nd
  levels(nd2$id) <- c(levels(nd2$id), "zzz")
  nd2$id[] <- "zzz"
  pn <- predict(fw, newdata = nd2, resp = "y1", allow_new_levels = TRUE)
  expect_true(all(is.finite(pn)))
  expect_equal(length(unique(round(pn, 10))), 1L)
})
