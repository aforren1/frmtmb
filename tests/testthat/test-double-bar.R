# Every `||` term as brms 2.23.0 reads it: (x | gr(g, cor = FALSE)), one
# term coded with the formula's own contrasts, an sd per column and no
# correlation (brms:::split_re_terms()). Lane dblbar's probe compared 35
# spellings with brms's standata() column for column
# (dev/dblbar-compare.R, dev/dblbar-findings.md). Seen to fail on the
# 0.70.0 build: a factor had lme4's column per level, (1 + x || g) was
# two blocks, a diag block named zero correlations, gr(cor = ) and
# mm(cor = ) were refused, and (x || p | g) read p as a variable.

dbar_data <- function(n = 240, seed = 7) {
  set.seed(seed)
  g <- factor(sample(sprintf("g%02d", 1:16), n, TRUE))
  f <- factor(sample(c("a", "b", "c"), n, TRUE))
  d <- data.frame(
    g = g, f = f,
    hg = factor(ifelse(as.integer(g) <= 8L, "u", "v")),
    g1 = factor(sample(sprintf("p%d", 1:5), n, TRUE)),
    g2 = factor(sample(sprintf("q%d", 1:4), n, TRUE)),
    m1 = factor(sample(sprintf("m%02d", 1:10), n, TRUE)),
    m2 = factor(sample(sprintf("m%02d", 1:10), n, TRUE)),
    x = stats::rnorm(n), z = stats::runif(n, -1, 1))
  d$y <- 1 + 0.5 * d$x - 0.3 * d$z + c(0, 0.4, -0.2)[as.integer(f)] +
    stats::rnorm(16, 0, 0.8)[as.integer(g)] +
    matrix(stats::rnorm(48, 0, 0.6), 16)[cbind(as.integer(g),
                                               as.integer(f))] +
    stats::rnorm(n, 0, 0.7)
  d$y2 <- d$y + stats::rnorm(n)
  d
}

# The cases, with the sd names brms 2.23.0 gives them, from
# dev/dblbar-log/brms-probe.txt
dbar_cases <- list(
  list(f = bf(y ~ f + (1 + f || g)),
       sd = c("g__Intercept", "g__fb", "g__fc")),
  list(f = bf(y ~ f + (f || g)), sd = c("g__Intercept", "g__fb", "g__fc")),
  list(f = bf(y ~ f + (0 + f || g)), sd = c("g__fa", "g__fb", "g__fc")),
  list(f = bf(y ~ z + f + (1 + z + f || g)),
       sd = c("g__Intercept", "g__z", "g__fb", "g__fc")),
  list(f = bf(y ~ x * f + (1 + x * f || g)),
       sd = c("g__Intercept", "g__x", "g__fb", "g__fc", "g__x:fb",
              "g__x:fc")),
  list(f = bf(y ~ x + f:x + (1 + f:x || g)),
       sd = c("g__Intercept", "g__fa:x", "g__fb:x", "g__fc:x")),
  list(f = bf(y ~ x + (1 + x || g)), sd = c("g__Intercept", "g__x")),
  list(f = bf(y ~ f + (1 + f || gr(g, by = hg))),
       sd = paste0("g__", c("Intercept", "fb", "fc"), ":hg",
                   rep(c("u", "v"), each = 3))),
  list(f = bf(y ~ f + (1 + f || g1:g2)),
       sd = c("g1:g2__Intercept", "g1:g2__fb", "g1:g2__fc")),
  list(f = bf(y ~ f + (1 + f || g1/g2)),
       sd = c("g1__Intercept", "g1__fb", "g1__fc", "g1:g2__Intercept",
              "g1:g2__fb", "g1:g2__fc")),
  list(f = bf(y ~ f + (1 + f || mm(m1, m2))),
       sd = c("mmm1m2__Intercept", "mmm1m2__fb", "mmm1m2__fc")),
  list(f = bf(y ~ x + (1 | g1), sigma ~ (1 + f || g)),
       sd = c("g1__Intercept", "g__sigma_Intercept", "g__sigma_fb",
              "g__sigma_fc")),
  list(f = bf(y ~ a + b * x, a ~ 1 + (1 + f || g), b ~ 1, nl = TRUE),
       sd = c("g__a_Intercept", "g__a_fb", "g__a_fc")),
  list(f = bf(y ~ 0 + f + (0 + f || g), cmc = FALSE),
       sd = c("g__fb", "g__fc")),
  list(f = mvbf(bf(y ~ f + (1 + f || g)), bf(y2 ~ x + (1 + x || g)),
                rescor = FALSE),
       sd = c("g__y_Intercept", "g__y_fb", "g__y_fc", "g__y2_Intercept",
              "g__y2_x"))
)

test_that("every || term has brms's columns and names, and no cor_", {
  dd <- dbar_data()
  for (cs in dbar_cases) {
    fit <- allow_warnings(frm(cs$f, family = gaussian(), data = dd),
                          c("Boundary", "singular"))
    v <- variables(fit)
    lab <- paste(cs$sd, collapse = " ")
    expect_setequal(grep("^sd_", v, value = TRUE), paste0("sd_", cs$sd))
    expect_identical(grep("^cor_", v, value = TRUE), character(0),
                     info = lab)
    expect_true(all(vapply(fit$frame$re_blocks, `[[`, "", "covstruct") %in%
                      c("diag", "us")), info = lab)
  }
})

test_that("a || term is one block, addressed by its factor", {
  dd <- dbar_data()
  fit <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd)
  expect_length(fit$frame$re_blocks, 1L)
  expect_identical(fit$frame$re_blocks[[1L]]$cnms,
                   c("(Intercept)", "fb", "fc"))
  expect_identical(colnames(ranef(fit)$g), c("Intercept", "fb", "fc"))
  vc <- VarCorr(fit)$g
  expect_identical(names(vc), "sd")
  expect_identical(rownames(vc$sd), c("Intercept", "fb", "fc"))
  expect_identical(rownames(summary(fit)$random$g),
                   c("sd(Intercept)", "sd(fb)", "sd(fc)"))
  cv <- confint_varcorr(fit)
  expect_false(any(cv$type == "cor"))
  num <- frm(bf(y ~ x + (1 + x || g)), family = gaussian(), data = dd)
  expect_identical(colnames(ranef(num)$g), c("Intercept", "x"))
})

test_that("gr(cor = FALSE), mm(cor = FALSE) and diag() are the || term", {
  dd <- dbar_data()
  fit_with <- function(form) frm(bf(form), family = gaussian(), data = dd)
  dbar <- fit_with(y ~ f + (1 + f || g))
  for (form in list(y ~ f + (1 + f | gr(g, cor = FALSE)),
                    y ~ f + diag(1 + f | g),
                    # || overrides gr()'s cor, as in brms
                    y ~ f + (1 + f || gr(g, cor = TRUE)))) {
    other <- fit_with(form)
    expect_identical(other$estimates, dbar$estimates)
    expect_identical(as.numeric(logLik(other)), as.numeric(logLik(dbar)))
  }
  # cor = TRUE is the default: the correlated term
  us <- fit_with(y ~ f + (1 + f | g))
  expect_identical(fit_with(y ~ f + (1 + f | gr(g, cor = TRUE)))$estimates,
                   us$estimates)
  mmd <- fit_with(y ~ f + (1 + f || mm(m1, m2)))
  mmc <- fit_with(y ~ f + (1 + f | mm(m1, m2, cor = FALSE)))
  expect_identical(mmc$estimates, mmd$estimates)
  expect_error(fit_with(y ~ f + (1 + f | gr(g, cor = "no"))),
               "must be TRUE or FALSE")
})

test_that("a factor in a || term is not lme4's column per level", {
  dd <- dbar_data()
  fit <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd)
  lme4 <- frm(bf(y ~ f + (1 | g) + diag(0 + f | g)), family = gaussian(),
              data = dd)
  ll <- as.numeric(logLik(fit))
  # different models: one sd more, and a likelihood the data tell apart
  expect_identical(length(grep("^sd_", variables(lme4))),
                   length(grep("^sd_", variables(fit))) + 1L)
  expect_gt(abs(ll - as.numeric(logLik(lme4))), 1e-6 * abs(ll))
})

test_that("an |ID| key beside || is refused, as brms refuses it", {
  dd <- dbar_data()
  f <- function(form) frm(bf(form), family = gaussian(), data = dd,
                          dry_run = "frame")
  expect_error(f(y ~ x + (1 + x || p | g)), "Invalid syntax")
  expect_error(f(y ~ x + (1 + x | p || g)), "Invalid syntax")
  expect_error(f(y ~ x + (1 + x | p | gr(g, cor = FALSE))),
               "drops the key")
  # the absent case: a key with correlation, and || beside a key
  expect_silent(f(y ~ x + (1 + x | p | gr(g, cor = TRUE))))
  expect_silent(f(y ~ x + (1 + x | p | g) + (1 + f || g1)))
})

test_that("|| over gr(cov = ) keeps numeric columns and refuses a factor", {
  dd <- dbar_data()
  A <- diag(16) * 0.5 + 0.5
  dimnames(A) <- list(levels(dd$g), levels(dd$g))
  f <- function(form) frm(bf(form), family = gaussian(), data = dd,
                          data2 = list(A = A), dry_run = "frame")
  # a known covariance correlates the columns of a term within a level,
  # so a numeric || term is one block per column, brms's model
  fr <- f(y ~ x + (1 + x || gr(g, cov = A)))
  expect_identical(vapply(fr$re_blocks, `[[`, 1L, "dim"), c(1L, 1L))
  expect_identical(vapply(fr$re_blocks, `[[`, "", "covstruct"),
                   c("gr_cov", "gr_cov"))
  expect_silent(f(y ~ x + (1 + x | gr(g, cov = A))))
  expect_error(f(y ~ f + (1 + f || gr(g, cov = A))), "own sd")
  expect_error(f(y ~ x + diag(1 + x | gr(g, cov = A))), "own sd")
  expect_error(f(y ~ x + ar1(f + 0 | gr(g, cov = A))),
               "do not combine with another structure")
  expect_error(f(y ~ x + homdiag(0 + x + z | gr(g, cov = A))),
               "do not combine with another structure")
  # one column needs no structure within the level, so any wrapper is
  # the plain known-covariance term (review of lane dblbar, item 3)
  one <- f(y ~ x + homdiag(0 + x | gr(g, cov = A)))
  ref <- f(y ~ x + (0 + x | gr(g, cov = A)))
  expect_identical(vapply(one$re_blocks, `[[`, "", "covstruct"), "gr_cov")
  expect_identical(one$re_blocks[[1L]]$dim, 1L)
  expect_identical(one$re_blocks[[1L]]$cnms, ref$re_blocks[[1L]]$cnms)
})

test_that("a Student-t || term has one latent scale per level", {
  # brms 2.23.0 multiplies every column of a cor = FALSE Student-t term
  # by one dfm_1 per level (dev/dblbar-log/brms-id.txt); the base split
  # (1 + x || gr(g, dist = "student")) into two diag_t blocks of one
  dd <- dbar_data()
  for (form in list(y ~ x + (1 + x || gr(g, dist = "student")),
                    y ~ f + (1 + f || gr(g, dist = "student")))) {
    fr <- frm(bf(form), family = gaussian(), data = dd, dry_run = "frame")
    expect_length(fr$re_blocks, 1L)
    expect_identical(fr$re_blocks[[1L]]$covstruct, "diag_t")
  }
  expect_identical(fr$re_blocks[[1L]]$cnms, c("(Intercept)", "fb", "fc"))
  fit <- allow_warnings(
    frm(bf(y ~ x + (1 + x || gr(g, dist = "student"))), family = gaussian(),
        data = dd), c("Boundary", "singular"))
  expect_setequal(grep("^sd_|^cor_", variables(fit), value = TRUE),
                  c("sd_g__Intercept", "sd_g__x"))
  expect_identical(names(VarCorr(fit)$g), "sd")
})

test_that("messages name a || term as it is written", {
  dd <- dbar_data()
  f <- function(form) frm(bf(form), family = gaussian(), data = dd,
                          dry_run = "frame")
  expect_error(f(y ~ f + (1 + f || g) + (1 | g)), "(1 + f || g)",
               fixed = TRUE)
  expect_error(f(y ~ f + (1 + f | gr(g, cor = FALSE)) + (1 | g)),
               "(1 + f | gr(g, cor = FALSE))", fixed = TRUE)
  # two spellings that rewrite to one bar keep their own names, in the
  # order written (review round 2: the last one named both sides)
  expect_error(f(y ~ f + (1 + f || g) + (1 + f | gr(g, cor = FALSE))),
               paste("appears in both (1 + f || g) and",
                     "(1 + f | gr(g, cor = FALSE))"), fixed = TRUE)
  expect_error(f(y ~ f + (1 + f || g) + diag(1 + f | g)),
               "appears in both (1 + f || g) and diag(1 + f | g)",
               fixed = TRUE)
  # a refused term over a known matrix is named as written, not by its
  # lme4 piece, and a structure keeps its arguments
  A <- diag(16) * 0.5 + 0.5
  dimnames(A) <- list(levels(dd$g), levels(dd$g))
  fa <- function(form) frm(bf(form), family = gaussian(), data = dd,
                           data2 = list(A = A), dry_run = "frame")
  expect_error(fa(y ~ f + (1 + f || gr(g, cov = A))),
               "(1 + f || gr(g, cov = A)) has a term with 3 columns",
               fixed = TRUE)
  expect_error(fa(y ~ x + rr(0 + f | gr(g, cov = A), d = 1)),
               "rr(0 + f | gr(g, cov = A), d = 1) has 3 columns",
               fixed = TRUE)
})

test_that("re_formula reads a || term's columns as the | term's", {
  dd <- dbar_data()
  fit <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd)
  nd <- dd[1:6, ]
  expect_identical(fitted(fit, newdata = nd, re_formula = ~ (1 + f || g)),
                   fitted(fit, newdata = nd))
  expect_identical(fitted(fit, newdata = nd, re_formula = ~ (1 + f | g)),
                   fitted(fit, newdata = nd))
  # the intercept alone: the population prediction plus the intercept's
  # conditional mode
  int <- fitted(fit, newdata = nd, re_formula = ~ (1 | g))[, 1]
  pop <- fitted(fit, newdata = nd, re_formula = NA)[, 1]
  u <- unname(ranef(fit)$g[as.character(nd$g), "Intercept"])
  expect_equal(int, pop + u, tolerance = 1e-12)
  # lme4's column fa is not a column of this term
  expect_error(fitted(fit, newdata = nd, re_formula = ~ (0 + f || g)),
               "matches no group-level term")
})

# brms's Z of block `id`, one n x levels matrix per column stacked in
# brms's column order; a multi-membership row is the weighted sum over
# its members, as brms's Stan code adds them
brms_block_Z <- function(sd, id) {
  n <- sd$N
  L <- sd[[paste0("N_", id)]]
  out <- list()
  k <- 1L
  repeat {
    z <- sd[[paste0("Z_", id, "_", k)]]
    M <- matrix(0, n, L)
    if (!is.null(z)) {
      M[cbind(seq_len(n), sd[[paste0("J_", id)]])] <- z
    } else if (!is.null(sd[[paste0("Z_", id, "_", k, "_1")]])) {
      m <- 1L
      while (!is.null(z <- sd[[paste0("Z_", id, "_", k, "_", m)]])) {
        A <- matrix(0, n, L)
        A[cbind(seq_len(n), sd[[paste0("J_", id, "_", m)]])] <-
          sd[[paste0("W_", id, "_", m)]] * z
        M <- M + A
        m <- m + 1L
      }
    } else {
      break
    }
    out[[k]] <- M
    k <- k + 1L
  }
  do.call(rbind, out)
}

# frmtmb's Z of one block in the same layout
frm_block_Z <- function(fr, bk) {
  Z <- as.matrix(fr$linpreds[[bk$components[[1L]]$lp_key]]$Z)
  do.call(rbind, lapply(seq_len(bk$dim), function(j) {
    Z[, bk$b_idx[(seq_len(bk$n_levels) - 1L) * bk$dim + j], drop = FALSE]
  }))
}

# equal up to the order of the levels, which the two packages sort
# apart for an interaction group
same_up_to_levels <- function(A, B) {
  if (!identical(dim(A), dim(B))) return(FALSE)
  key <- function(M) sort(apply(M, 2L, function(v) {
    paste(format(v, digits = 15L), collapse = ",")
  }))
  identical(key(A), key(B))
}

test_that("every || case has brms 2.23.0's Z columns and sd prior rows", {
  skip_unless_brms()
  dd <- dbar_data()
  cases <- c("y ~ f + (1 + f || g)", "y ~ f + (0 + f || g)",
             "y ~ f + (f || g)", "y ~ z + f + (1 + z + f || g)",
             "y ~ x * f + (1 + x * f || g)", "y ~ x + (1 + x || g)",
             "y ~ f + (1 + f || g1:g2)", "y ~ f + (1 + f || g1/g2)",
             "y ~ f + (1 + f || mm(m1, m2))",
             "y ~ f + (1 + f | gr(g, cor = FALSE))")
  for (txt in cases) {
    form <- stats::as.formula(txt)
    sd_b <- brms_standata(brms::bf(form), data = dd, family = gaussian())
    fr <- frm(bf(form), family = gaussian(), data = dd, dry_run = "frame")
    ids <- seq_along(grep("^M_", names(sd_b)))
    expect_length(fr$re_blocks, length(ids))
    hit <- vapply(fr$re_blocks, function(bk) {
      Zf <- frm_block_Z(fr, bk)
      m <- which(vapply(ids, function(id) {
        same_up_to_levels(Zf, brms_block_Z(sd_b, id))
      }, NA))
      if (length(m) == 1L) m else NA_integer_
    }, 1L)
    expect_false(anyNA(hit), info = txt)
    expect_false(anyDuplicated(hit) > 0L, info = txt)
    pb <- as.data.frame(brms::default_prior(brms::bf(form), data = dd,
                                            family = gaussian()))
    pf <- as.data.frame(default_prior(bf(form), data = dd,
                                      family = gaussian()))
    key <- function(p) {
      p <- p[p$class %in% c("sd", "cor") & !nzchar(p$coef), ]
      sort(paste(p$class, p$group, p$dpar, sep = "|"))
    }
    # a nested or multi-membership group is spelled apart in the prior
    # table on every term, | as well as || (dev/dblbar-findings.md)
    if (!grepl("/|mm[(]", txt)) {
      expect_identical(key(pf), key(pb), info = txt)
    }
    # an uncorrelated term has no correlation in either
    expect_false(any(pb$class == "cor"), info = txt)
    expect_false(any(pf$class == "cor"), info = txt)
  }
})

