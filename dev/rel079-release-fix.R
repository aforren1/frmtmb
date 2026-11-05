# The release fixes of the seven reviews' "Short check of fix round 2",
# applied by dev/rel079-resolve.R after the hand resolutions. Each edit
# matches once or the run stops. NEWS text is in dev/rel079-news-split.R,
# because the NEWS files are written after this. Per fix the review's
# exact text is used; every place it did not apply as written is listed in
# dev/round-20261104.md, "Where a review's text did not apply".

# ---- refuse --------------------------------------------------------------
# RF1: the Tierney-Kadane gate, then a second sweep per parameter. A fit under
# the gate pays 0.1 percent of its time and a nearly constant error stays quiet.
edit_file("R/frame.R",
          paste0("mi_laplace_row_errors <- function(fit, at) {\n",
                 "  env <- fit$obj$env\n",
                 "  par <- env$last.par.best\n"),
          paste0("mi_laplace_row_errors <- function(fit, at) {\n",
                 "  env <- fit$obj$env\n",
                 "  # last.par, so that this reads a moved parameter after fit$obj$fn(q)\n",
                 "  par <- env$last.par\n"))
edit_file("R/frame.R",
          paste0("#' The most read values `mi_laplace_note()` integrates: each is about 100\n",
                 "#' sweeps of the whole objective, so beyond this the rest are scaled from an\n",
                 "#' even sample. The share of a fit this costs is measured in\n",
                 "#' `dev/refuse-f2-cost.R`.\n",
                 "#'\n",
                 "#' @noRd\n",
                 "mi_laplace_max_rows <- 30L\n"),
          paste0("#' The most read values `mi_laplace_note()` integrates per parameter: each is\n",
                 "#' about 100 sweeps of the whole objective, so beyond this the rest are scaled\n",
                 "#' from an even sample. With 10 rows the second stage is about 37 percent of a\n",
                 "#' fit (dev/reviews/2026-11-02-refuse.md, RF1).\n",
                 "#'\n",
                 "#' @noRd\n",
                 "mi_laplace_max_rows <- 10L\n\n",
                 "#' Tierney-Kadane's first term over every read value: `f4 / (8 f2^2) - 5 f3^2 /\n",
                 "#' (24 f2^3)`, from three taped inner Hessians. It tracks `integrate()` on smooth\n",
                 "#' families and fails near a kink, so it gates and does not decide.\n",
                 "#'\n",
                 "#' @noRd\n",
                 "mi_laplace_tk <- function(fit) {\n",
                 "  env <- fit$obj$env\n",
                 "  par <- env$last.par.best\n",
                 "  mp <- which(names(par) == \"miss\")\n",
                 "  ix <- match(mp, env$random)\n",
                 "  d2 <- function(p) Matrix::diag(env$spHess(p, random = TRUE))[ix]\n",
                 "  f2 <- d2(par)\n",
                 "  # in units of each value's own width, so that a well-pinned value is read\n",
                 "  # on its own scale\n",
                 "  st <- 0.2 / sqrt(pmax(f2, 1e-8))\n",
                 "  up <- par\n",
                 "  up[mp] <- up[mp] + st\n",
                 "  dn <- par\n",
                 "  dn[mp] <- dn[mp] - st\n",
                 "  fp <- d2(up)\n",
                 "  fm <- d2(dn)\n",
                 "  f3 <- (fp - fm) / (2 * st)\n",
                 "  f4 <- (fp - 2 * f2 + fm) / st^2\n",
                 "  e <- f4 / (8 * f2^2) - 5 * f3^2 / (24 * f2^3)\n",
                 "  e[!is.finite(e)] <- 0\n",
                 "  sum(e)\n",
                 "}\n\n",
                 "#' The share of a `circle` chart's read rows near the end of the window, or\n",
                 "#' `NA`. There the value has two modes, which the gate's Hessians barely see: a\n",
                 "#' move of y's parameters hardly moves a two-mode total, and the second stage\n",
                 "#' alone misses a shift of 7.34 standard errors. It reads 0.97 to 1.00 on 36 fits\n",
                 "#' near the end, 0 on 12 at mean 0 (dev/reviews/2026-11-02-refuse.md, RF2).\n",
                 "#'\n",
                 "#' @noRd\n",
                 "mi_laplace_window_share <- function(fit, r, rows) {\n",
                 "  tryCatch({\n",
                 "    dp <- eval_dpars(fit)[[r]]\n",
                 "    n <- max(rows)\n",
                 "    mu <- rep_len(as.numeric(dp[[\"mu\"]]), max(n, length(dp[[\"mu\"]])))[rows]\n",
                 "    k <- rep_len(as.numeric(dp[[\"kappa\"]]), max(n, length(dp[[\"kappa\"]])))[rows]\n",
                 "    mean(pi - abs(mu) < 2 / sqrt(k))\n",
                 "  }, error = function(e) NA_real_)\n",
                 "}\n\n",
                 "#' The outer parameters a read response's own density carries: each `betad`\n",
                 "#' of it, or its `mu` intercept if it has none, as exponential.\n",
                 "#'\n",
                 "#' @noRd\n",
                 "mi_laplace_own_pars <- function(fit, who) {\n",
                 "  mp <- outer_par_map(fit)\n",
                 "  out <- integer()\n",
                 "  for (r in who) {\n",
                 "    pre <- paste0(r, \"_\")\n",
                 "    own <- which(mp$comp == \"betad\" & startsWith(mp$names, pre))\n",
                 "    if (!length(own)) {\n",
                 "      own <- which(mp$comp == \"beta\" &\n",
                 "                     mp$names == paste0(pre, \"(Intercept)\"))\n",
                 "    }\n",
                 "    out <- c(out, own)\n",
                 "  }\n",
                 "  stats::setNames(out, mp$names[out])\n",
                 "}\n"))
x <- rd("R/frame.R")
a <- which(startsWith(x, "mi_laplace_note <- function(fit) {"))
b <- a + which(x[(a + 1L):length(x)] == "}")[1L]
stopifnot(length(a) == 1L)
note <- c(
  "mi_laplace_note <- function(fit) {",
  "  act <- fit$control$check_mi %||% \"warning\"",
  "  if (identical(act, \"ignore\")) return(fit)",
  "  mp <- fit$frame[[\"mi_map\"]] %||% list()",
  "  n_read <- sum(vapply(mp, function(m) length(m[[\"rows\"]]), 0L))",
  "  if (!n_read || mi_laplace_exact(fit$frame)) return(fit)",
  "  who <- names(mp)[vapply(mp, function(m) length(m[[\"rows\"]]) > 0L, NA)]",
  "  head <- paste0(\"The Laplace approximation of the \", n_read,",
  "                 \" missing value(s) of \", paste0(\"`\", who, \"`\", collapse = \", \"),",
  "                 \" that another density reads \")",
  "  tail <- paste0(\"frmtmb.sample::frm_sample() has no such error. If the \",",
  "                 \"predictor that reads the value carries little about it, \",",
  "                 \"deleting those rows gives nearly the same model without the \",",
  "                 \"error. See ?bf, Missing values under mi()\")",
  "  say <- function(msg) {",
  "    if (identical(act, \"stop\")) {",
  "      frm_stop(msg, call. = FALSE, class = \"frmtmb_mi_laplace_error\")",
  "    }",
  "    frm_warning(msg, call. = FALSE, class = \"frmtmb_mi_laplace_error\")",
  "  }",
  "  # the window's end has two modes, a reason a measurement of one mode misses",
  "  for (r in who) {",
  "    if (!identical(mp[[r]][[\"chart\"]], \"circle\")) next",
  "    sh <- mi_laplace_window_share(fit, r, mp[[r]][[\"rows\"]])",
  "    if (isTRUE(sh >= 0.1)) {",
  "      fit$mi_laplace <- list(n_read = n_read, window = sh)",
  "      say(paste0(head, \"expands about one mode, and \",",
  "                 format(100 * sh, digits = 2), \" percent of the rows of `\", r,",
  "                 \"` it reads have a mean within 2 / sqrt(kappa) of the end of \",",
  "                 \"the window, where the value has two. The error there is -1 \",",
  "                 \"to -4.2 per value, and estimates can be off by up to 17 \",",
  "                 \"standard errors. \", tail))",
  "      return(fit)",
  "    }",
  "  }",
  "  tk <- tryCatch(mi_laplace_tk(fit), error = function(e) Inf)",
  "  fit$mi_laplace <- list(tk = tk, n_read = n_read)",
  "  if (abs(tk) <= 1) return(fit)",
  "  at <- if (n_read > mi_laplace_max_rows) {",
  "    unique(round(seq(1, n_read, length.out = mi_laplace_max_rows)))",
  "  } else {",
  "    seq_len(n_read)",
  "  }",
  "  total <- function(q) {",
  "    suppressWarnings(fit$obj$fn(q))",
  "    err <- tryCatch(mi_laplace_row_errors(fit, at),",
  "                    error = function(e) rep(NA_real_, length(at)))",
  "    if (any(is.finite(err))) mean(err, na.rm = TRUE) * n_read else NA_real_",
  "  }",
  "  p <- fit$opt$par",
  "  # the inner mode and last.par come back whatever a moved parameter leaves",
  "  on.exit(suppressWarnings(fit$obj$fn(fit$opt$par)), add = TRUE)",
  "  own <- mi_laplace_own_pars(fit, who)",
  "  se <- tryCatch(",
  "    sqrt(diag(solve(stats::optimHess(p, fit$obj$fn, fit$obj$gr)))),",
  "    error = function(e) rep(NA_real_, length(p)))",
  "  se[!is.finite(se)] <- 0.1",
  "  t0 <- total(p)",
  "  moves <- vapply(own, function(j) {",
  "    q <- p",
  "    q[j] <- q[j] + se[j]",
  "    abs(total(q) - t0)",
  "  }, 0)",
  "  fit$mi_laplace <- list(tk = tk, n_read = n_read, rows = at, total = t0,",
  "                         moves = moves)",
  "  if (any(is.finite(moves)) && max(moves, na.rm = TRUE) >= mi_laplace_limit) {",
  "    k <- which.max(moves)",
  "    say(paste0(head, \"moves `\", names(moves)[k], \"` by about \",",
  "               format(moves[k], digits = 2), \" standard errors, by integrate() \",",
  "               \"over the values at the estimate and at one standard error \",",
  "               \"above it. \", tail))",
  "  }",
  "  fit",
  "}")
x <- c(x[seq_len(a - 1L)], note, x[-seq_len(b)])
wr(x, "R/frame.R")
cat("rewrote mi_laplace_note()\n")
edit_file("R/frame.R",
          paste0("#' Calibrated in `dev/refuse-findings.md` against the shifts in exact\n",
                 "#' standard errors of a read value's approximation.\n",
                 "#'\n",
                 "#' @noRd\n",
                 "mi_laplace_limit <- 1\n"),
          paste0("#' A move of 1 standard error changes the objective's slope by 1, so a total\n",
                 "#' that moves by d predicts a shift of d standard errors: 75 of 76 shifts of 1\n",
                 "#' or more, and 3 of 39 below 0.5 (dev/reviews/2026-11-02-refuse.md, RF1).\n",
                 "#'\n",
                 "#' @noRd\n",
                 "mi_laplace_limit <- 0.5\n"))
edit_file("R/fit.R",
          "#' @param check_se What to do when a converged fit has parameters with\n",
          paste0("#' @param check_mi What to do when the Laplace approximation of a read `mi()`\n",
                 "#'   value moves an estimate: `\"warning\"` (default), `\"ignore\"`, or `\"stop\"`.\n",
                 "#'   See [bf()], section \"Missing values under `mi()`\". `\"ignore\"` skips the\n",
                 "#'   check.\n",
                 "#' @param check_se What to do when a converged fit has parameters with\n"))
edit_file("R/fit.R",
          "                           check_se = c(\"warning\", \"ignore\", \"stop\"),\n",
          paste0("                           check_se = c(\"warning\", \"ignore\", \"stop\"),\n",
                 "                           check_mi = c(\"warning\", \"ignore\", \"stop\"),\n"))
edit_file("R/fit.R",
          "       check_se = frm_match_arg(check_se),\n",
          "       check_se = frm_match_arg(check_se),\n       check_mi = frm_match_arg(check_mi),\n")
x <- rd("R/bf.R")
a <- grep("`von_mises` loses about$", x)
b <- grep("^#' the same model without the error\\.$", x)
stopifnot(length(a) == 1L, length(b) == 1L, b > a)
st <- grep("^#' \\+0\\.008 `frechet`, -0\\.125 `beta`, -0\\.079$", x)
cat("?bf old paragraph:\n", paste0(x[a:b], "\n"), sep = "")
a1 <- sub("`von_mises` loses about$", "`von_mises` loses -1 to", x[a])
x <- c(x[seq_len(a - 1L)], a1,
  "#' -4.2 per value at a mean near the end of its window, where the value has two",
  "#' modes, and estimates there can be off by up to 17 standard errors.",
  "#'",
  "#' After the fit, [frm()] warns, with class `frmtmb_mi_laplace_error`, when the",
  "#' approximation moves an estimate. Three inner Hessians over every read value",
  "#' give Tierney-Kadane's term, a gate that costs 0.1 percent of a fit. Past it, an",
  "#' `integrate()` over 10 of the values at the estimate and at one standard error",
  "#' above it, per parameter of the read response, tells how far a parameter moves:",
  "#' a move of 0.5 or more warns, and `fit$mi_laplace` holds the moves. That stage",
  "#' costs about 37 percent of the fit. For a `von_mises` value, 10 percent or more",
  "#' of read rows with `pi - abs(mu) < 2 / sqrt(kappa)` warn with no integral.",
  "#' `check_mi` in [frmtmb_control()] sets what happens. The shift from an exact",
  "#' maximum, in standard errors, over 144 fits: a shift of 1 or more warned on 75 of",
  "#' 76, and a shift below 0.5 on 3 of 39, and the largest silent shift was 1.01.",
  "#'",
  "#' `frmtmb.sample::frm_sample()` holds every value and has no such error. Its",
  "#' draws name a value as brms does, `Ymi_<resp>[k]`, on the response's scale. If",
  "#' the reader carries little about the value, deleting those rows gives nearly",
  x[b:length(x)][1L],
  x[-seq_len(b)])
wr(x, "R/bf.R")
cat("rewrote the ?bf paragraph\n")
# RF3
edit_file("R/predict.R",
          paste0("#' against brms's 0.3106. A mixture and a `mu` with a `zi` that read one term do\n",
                 "#' the same.\n"),
          paste0("#' against brms's 0.3106. A mixture and a `mu` with a `zi` that read one term do\n",
                 "#' the same.\n#'\n",
                 "#' At a row whose `mi()` value left the objective, a reader's fit reads the\n",
                 "#' value's mean. Under a nonlinear link this is below brms's average: a median\n",
                 "#' of 25% for a log link on a Gamma value of shape 1.2.\n"))
# RF4: brms reads the first flag of a key written id =; a |p| term with cor =
# FALSE has dropped its key in brms as well, and builds two blocks.
edit_file("R/parse.R",
          paste0("      if (!is.null(attr(read, \"id_no_cor\"))) {\n",
                 "        id_no_cor[[length(id_no_cor) + 1L]] <- list(\n",
                 "          label = attr(read, \"id_no_cor\"), group = bar_key_group(read[[3L]], env),\n",
                 "          by_id = bar_has_id_arg(inner),\n",
                 "          text = orig)\n"),
          paste0("      if (!is.null(attr(read, \"id_no_cor\")) && bar_has_id_arg(inner)) {\n",
                 "        id_no_cor[[length(id_no_cor) + 1L]] <- list(\n",
                 "          label = attr(read, \"id_no_cor\"), group = bar_key_group(read[[3L]], env),\n",
                 "          by_id = TRUE,\n",
                 "          text = orig)\n"))
edit_file("R/parse.R",
          "          text = orig)\n        attr(read, \"id_no_cor\") <- NULL\n      }\n",
          "          text = orig)\n      }\n      attr(read, \"id_no_cor\") <- NULL\n")
edit_file("tests/testthat/test-mm-id-cov.R",
          paste0("  expect_error(f(bf(y ~ x + (1 | p | mm(g1, g2)),\n",
                 "                    sigma ~ 1 + (1 | p | mm(g1, g2, cor = FALSE)))), msg)\n",
                 "  expect_error(f(bf(y ~ x + (1 | p | mm(g1, g2, cor = FALSE)),\n",
                 "                    sigma ~ 1 + (1 | p | mm(g1, g2)))), msg)\n"),
          paste0("  # the |p| spelling drops its key under cor = FALSE, as brms's program does,\n",
                 "  # with two uncorrelated blocks in either order\n",
                 "  fp <- f(bf(y ~ x + (1 | p | mm(g1, g2)),\n",
                 "             sigma ~ 1 + (1 | p | mm(g1, g2, cor = FALSE))))\n",
                 "  expect_length(fp$re_blocks, 2L)\n",
                 "  fp <- f(bf(y ~ x + (1 | p | mm(g1, g2, cor = FALSE)),\n",
                 "             sigma ~ 1 + (1 | p | mm(g1, g2))))\n",
                 "  expect_length(fp$re_blocks, 2L)\n",
                 "  fg <- f(bf(y ~ x + (1 | p | gr(g1, cor = FALSE)),\n",
                 "             sigma ~ 1 + (1 | p | g1)))\n",
                 "  expect_length(fg$re_blocks, 2L)\n"))
# the tests of the old rule's fields, on the new rule
tf <- "tests/testthat/test-mi-rows-read.R"
edit_file(tf,
          paste0("  expect_gt(abs(r$fit$mi_laplace$total), frmtmb:::mi_laplace_limit)\n",
                 "  expect_identical(r$fit$mi_laplace$n_read, 60L)\n",
                 "  expect_length(r$fit$mi_laplace$error, frmtmb:::mi_laplace_max_rows)\n"),
          paste0("  expect_gt(abs(r$fit$mi_laplace$tk), 1)\n",
                 "  expect_gte(max(r$fit$mi_laplace$moves), frmtmb:::mi_laplace_limit)\n",
                 "  expect_identical(r$fit$mi_laplace$n_read, 60L)\n",
                 "  expect_length(r$fit$mi_laplace$rows, frmtmb:::mi_laplace_max_rows)\n"))
edit_file(tf,
          paste0("  expect_match(msg, \"deleting those rows\", fixed = TRUE)\n",
                 "})\n"),
          paste0("  expect_match(msg, \"deleting those rows\", fixed = TRUE)\n",
                 "  expect_match(msg, \"standard errors\", fixed = TRUE)\n",
                 "  expect_match(msg, \"`y_\", fixed = TRUE)\n",
                 "})\n\n",
                 "test_that(\"check_mi sets what the check does\", {\n",
                 "  d <- laplace_data(0.05, 60)\n",
                 "  f <- bf(y | mi() ~ x) + exgaussian() + bf(z ~ mi(y) + x) + gaussian()\n",
                 "  off <- suppressWarnings(\n",
                 "    frm(f, data = d, control = frmtmb_control(check_mi = \"ignore\")))\n",
                 "  expect_null(off$mi_laplace)\n",
                 "  expect_error(\n",
                 "    suppressWarnings(\n",
                 "      frm(f, data = d, control = frmtmb_control(check_mi = \"stop\"))),\n",
                 "    class = \"frmtmb_mi_laplace_error\")\n",
                 "  expect_error(frmtmb_control(check_mi = \"no\"))\n",
                 "})\n\n",
                 "test_that(\"the parameters come back after the second stage\", {\n",
                 "  r <- laplace_fit(laplace_data(0.05, 60))\n",
                 "  env <- r$fit$obj$env\n",
                 "  expect_equal(unname(env$last.par[env$lfixed()]),\n",
                 "               unname(r$fit$opt$par))\n",
                 "})\n\n",
                 "test_that(\"a von_mises value near the window's end warns by name\", {\n",
                 "  set.seed(11600)\n",
                 "  n <- 200\n",
                 "  x <- stats::rnorm(n)\n",
                 "  y <- 3.14 + 0.05 * x + stats::rnorm(n, 0, 0.2)\n",
                 "  y <- ((y + pi) %% (2 * pi)) - pi\n",
                 "  d <- data.frame(x = x, y = y)\n",
                 "  d$z <- 0.5 + 0.7 * y + 0.3 * x + stats::rnorm(n, 0, 0.7)\n",
                 "  d$y[sample(n, 15)] <- NA\n",
                 "  said <- character()\n",
                 "  fit <- withCallingHandlers(\n",
                 "    frm(bf(y | mi() ~ x) + von_mises() + bf(z ~ mi(y) + x) + gaussian(),\n",
                 "        data = d),\n",
                 "    frmtmb_mi_laplace_error = function(w) {\n",
                 "      said <<- c(said, conditionMessage(w))\n",
                 "      invokeRestart(\"muffleWarning\")\n",
                 "    },\n",
                 "    warning = function(w) invokeRestart(\"muffleWarning\"))\n",
                 "  expect_length(said, 1L)\n",
                 "  expect_match(said, \"frm_sample()\", fixed = TRUE)\n",
                 "  expect_match(said, \"two\", fixed = TRUE)\n",
                 "  expect_gte(fit$mi_laplace$window, 0.1)\n",
                 "})\n"))
edit_file(tf,
          paste0("  expect_false(r$warned)\n",
                 "  expect_lt(abs(r$fit$mi_laplace$total), frmtmb:::mi_laplace_limit)\n"),
          paste0("  expect_false(r$warned)\n",
                 "  ml <- r$fit$mi_laplace\n",
                 "  expect_true(abs(ml$tk) <= 1 ||\n",
                 "                max(ml$moves) < frmtmb:::mi_laplace_limit)\n"))

# ---- impsamp -------------------------------------------------------------
# RF1 is NEWS (dev/rel079-news-split.R). The findings are the lane's record.
edit_file("R/fit.R",
          "\"in the same place\"))",
          paste0("\"in the same place. If neither settles it, check whether the \",\n",
                 "        \"Laplace fit converged and whether the groups have enough rows to \",\n",
                 "        \"identify every variance component\"))"))
edit_file("vignettes/diagnostics.Rmd",
          "sizes, and raise the count only if they ask for it.\n",
          paste0("sizes, and raise the count only if they ask for it.\n\n",
                 "Center a column whose mean is more than about 1e4 times its standard\n",
                 "deviation. At 2e4 times, six of eight ordinal fits warned of false\n",
                 "convergence, with their standard errors not available, and one moved a\n",
                 "slope by 0.13 standard errors with no warning.\n"))
x <- rd("R/families.R")
a <- grep("^#' @param extra_pars Optional function", x)
stopifnot(length(a) == 1L)
b <- a + which(grepl("^#' @param ", x[(a + 1L):length(x)]) |
                 x[(a + 1L):length(x)] == "#'")[1L] - 1L
x <- c(x[seq_len(b)],
       "#'   Keep each on a scale near 1: `frm(importance =)` differentiates by steps of",
       "#'   1e-4 on it.",
       x[-seq_len(b)])
wr(x, "R/families.R")
cat("edited R/families.R\n")

# ---- latres --------------------------------------------------------------
# RF1: the Tierney-Kadane gate below the row limit; the lane's condition stays
# above it, where the dense v of the correction is not paid.
x <- rd("R/autocor-latent.R")
a <- which(startsWith(x, "latres_laplace_reason <- function(fit, bk) {"))
stopifnot(length(a) == 1L)
tkfn <- c(
  "#' The first-order Laplace correction of a fit at outer parameters `par`:",
  "#' `-1/8 sum(h4 v^2) + 1/8 sum(h3 h3 v v v) + 1/12 sum(h3 h3 v^3)`, with `h3`, `h4`",
  "#' the row's derivatives of `-log density` in `eta` and `v = A H^-1 A'` the",
  "#' covariance of the latents at the inner mode. Read at the fit it is the exact",
  "#' value's distance from Laplace: median 0.340 to 0.012 on designs no side had",
  "#' seen (dev/reviews/2026-11-02-latres.md, RF1).",
  "#'",
  "#' @noRd",
  "latres_tk <- function(fit, par, bk, map) {",
  "  env <- fit$obj$env",
  "  if (!is.finite(suppressWarnings(fit$obj$fn(par)))) return(NA_real_)",
  "  full <- env$last.par",
  "  at <- fit_set_outer(fit, par, map)",
  "  at$estimates[[\"b\"]] <- as.numeric(full[env$random])",
  "  a <- bk[[\"aux\"]]",
  "  r <- a[[\"ac\"]][[\"resp\"]]",
  "  fam <- at$spec$responses[[r]]$family",
  "  dp <- eval_dpars(at)[[r]]",
  "  lk <- at$frame[[\"linpreds\"]][[linpred_key(r, \"mu\")]][[\"link\"]]",
  "  eta <- lk$linkfun(as.numeric(dp[[\"mu\"]]))",
  "  y <- at$frame[[\"y\"]][[r]]",
  "  av <- at$frame[[\"aterm_values\"]][[r]]",
  "  ex <- fit_extras(at)",
  "  ll <- function(e) {",
  "    dp[[\"mu\"]] <- lk$linkinv(e)",
  "    as.numeric(row_lpdf(fam, y, y, dp, av, ex))",
  "  }",
  "  h <- 0.02 * latres_marginal_sd(at$estimates[[\"theta\"]][bk[[\"theta_idx\"]]], a)",
  "  l0 <- ll(eta)",
  "  lp <- ll(eta + h)",
  "  ln <- ll(eta - h)",
  "  lp2 <- ll(eta + 2 * h)",
  "  ln2 <- ll(eta - 2 * h)",
  "  h3 <- -(lp2 - 2 * lp + 2 * ln - ln2) / (2 * h^3)",
  "  h4 <- -(lp2 - 4 * lp + 6 * l0 - 4 * ln + ln2) / h^4",
  "  pj <- latres_jacobian(bk, at$estimates[[\"theta\"]][bk[[\"theta_idx\"]]])",
  "  A <- Matrix::sparseMatrix(i = pj$i, j = pj$j, x = pj$x,",
  "                            dims = c(a[[\"n\"]], length(env$random)))",
  "  H <- env$spHess(full, random = TRUE)",
  "  # TMB fills one template's values, and Matrix keeps a factor in its slot, so",
  "  # Cholesky() on it hands back the factor of an earlier point: a residual of",
  "  # 0.037 here, where C read 0.97 for 1.10",
  "  H@factors <- list()",
  "  Ch <- Matrix::Cholesky(H, LDL = FALSE)",
  "  v <- as.matrix(A %*% Matrix::solve(Ch, Matrix::t(A)))",
  "  dv <- diag(v)",
  "  -sum(h4 * dv^2) / 8 + sum(outer(h3 * dv, h3 * dv) * v) / 8 +",
  "    sum(outer(h3, h3) * v^3) / 12",
  "}",
  "",
  "#' How far that correction moves log `sderr` and `ar[1]`: `V grad(C)`, the",
  "#' gradient by differences of 0.05 in each outer parameter but `beta`. `Inf` where",
  "#' the fit has no covariance to read it by. The tape's last point goes back on",
  "#' exit, or a reader of `env$last.par` reads the last step: 1.5 nats on \"the",
  "#' readers agree with the tape\".",
  "#'",
  "#' @noRd",
  "latres_shift_tk <- function(fit, bk) {",
  "  par <- fit$opt$par",
  "  map <- outer_par_map(fit)",
  "  env <- fit$obj$env",
  "  lp0 <- env$last.par",
  "  lb0 <- env$last.par.best",
  "  on.exit({",
  "    env$last.par <- lp0",
  "    env$last.par.best <- lb0",
  "  })",
  "  V <- tryCatch(as.matrix(sdr_of(fit)$cov.fixed), error = function(e) NULL)",
  "  if (is.null(V) || nrow(V) != length(par)) return(Inf)",
  "  g <- numeric(length(par))",
  "  for (k in which(names(par) != \"beta\")) {",
  "    pp <- par",
  "    pp[k] <- pp[k] + 0.05",
  "    pm <- par",
  "    pm[k] <- pm[k] - 0.05",
  "    g[k] <- (latres_tk(fit, pp, bk, map) - latres_tk(fit, pm, bk, map)) / 0.1",
  "  }",
  "  g[!is.finite(g)] <- 0",
  "  th <- which(names(par) == \"theta\")[bk[[\"theta_idx\"]][1:2]]",
  "  if (anyNA(th) || !all(is.finite(V[th, ]))) return(Inf)",
  "  d <- as.numeric(V[th, ] %*% g)",
  "  t2 <- par[th[2]]",
  "  max(abs(d[1]),",
  "      abs((t2 + d[2]) / sqrt(1 + (t2 + d[2])^2) - t2 / sqrt(1 + t2^2)))",
  "}",
  "")
# the doc block above latres_laplace_reason stays with it
da <- a
while (startsWith(x[da - 1L], "#'")) da <- da - 1L
x <- c(x[seq_len(da - 1L)], tkfn, x[da:(a - 1L)],
       "latres_laplace_reason <- function(fit, bk) {",
       "  r <- bk[[\"aux\"]][[\"ac\"]][[\"resp\"]]",
       "  fam <- fit$spec$responses[[r]]$family",
       "  # above the limit the correction's dense v is the cost the limit avoids, so",
       "  # the row-level condition stays there",
       "  if (fit$frame[[\"n_obs\"]] > getOption(\"frmtmb.latres_se_rows\", 3000L)) {",
       "    return(latres_share_reason(fit, bk, fam, r))",
       "  }",
       "  # its error is not local: the shift is 0.02 to 0.06 with the exact MLE far off",
       "  # along a flat ridge, so the correction cannot read it",
       "  if (identical(fam[[\"type\"]], \"ordinal\")) return(\"each row is one category\")",
       "  s <- tryCatch(latres_shift_tk(fit, bk), error = function(e) Inf)",
       "  if (s <= 0.15) return(NULL)",
       "  if (!is.finite(s)) {",
       "    return(\"the fit has no curvature to read its first-order correction by\")",
       "  }",
       "  paste0(\"the first-order correction moves log(sderr) or ar[1] by \",",
       "         format(s, digits = 2), \", above 0.15\")",
       "}",
       "",
       "#' The row-level condition above the row limit.",
       "#'",
       "#' @noRd",
       "latres_share_reason <- function(fit, bk, fam, r) {",
       x[(a + 3L):length(x)])
wr(x, "R/autocor-latent.R")
cat("edited R/autocor-latent.R: the gate\n")
# RF4
edit_file("R/autocor-latent.R",
          "    other <- !state$na & !is.finite(state$se)\n",
          paste0("    other <- !state$na & !is.finite(state$se)\n",
                 "    # beside a chain at its limit a row can read NA with a finite, large standard\n",
                 "    # error, and no sentence named it\n",
                 "    unsaid <- state$na & !lim\n"))
edit_file("R/autocor-latent.R",
          paste0("        paste0(paste(nm[other], collapse = \", \"),\n",
                 "               \" has no standard error either. \")\n",
                 "      },\n"),
          paste0("        paste0(paste(nm[other], collapse = \", \"),\n",
                 "               \" has no standard error either. \")\n",
                 "      },\n",
                 "      if (chain && any(unsaid)) {\n",
                 "        paste0(paste(nm[unsaid], collapse = \", \"),\n",
                 "               \" is not known there either. \")\n",
                 "      },\n"))
# RF2: reading the joint covariance of a latent residual is n^2 dense
edit_file("R/predict.R",
          paste0("  if (!is.null(cache$Vjoint)) return(cache$Vjoint)\n",
                 "  Q <- joint_precision(fit)\n"),
          paste0("  if (!is.null(cache$Vjoint)) return(cache$Vjoint)\n",
                 "  if (length(frame_latres(fit$frame))) {\n",
                 "    m <- fit$frame[[\"n_obs\"]]\n",
                 "    lim <- getOption(\"frmtmb.latres_dense_rows\", 20000L)\n",
                 "    if (m > lim) {\n",
                 "      frm_stop(\"This standard error reads the whole joint covariance of a \",\n",
                 "               \"latent residual as a dense matrix, about \",\n",
                 "               format(16 * m^2 / 1e9, digits = 2), \" GB at \", m, \" rows. \",\n",
                 "               \"options(frmtmb.latres_dense_rows = ) moves the limit of \",\n",
                 "               lim, \" rows; frm_linpred(), fitted() on new rows, emmeans() \",\n",
                 "               \"and predict() do not read it\", call. = FALSE)\n",
                 "    }\n",
                 "  }\n",
                 "  Q <- joint_precision(fit)\n"))
# RF3: the factor by construction, so that a chain at the numeric edge still
# has a point in brms's coordinates
edit_file("tests/testthat/helper-brms.R",
          paste0("brms_latres_chol <- function(a, th, k, J = NULL) {\n",
                 "  cf <- latres_coefs(th, a)\n"),
          paste0("brms_latres_chol <- function(a, th, k, J = NULL) {\n",
                 "  if (identical(a$kind, \"unstr\") && identical(as.integer(J), seq_len(a$d))) {\n",
                 "    return(as.matrix(us_chol_L(th, a$d)))\n",
                 "  }\n",
                 "  cf <- latres_coefs(th, a)\n"))
edit_file("tests/testthat/test-latres-brms.R",
          paste0("  at_edge <- 0L\n",
                 "  for (tm in lrb_terms) {\n"),
          "  for (tm in lrb_terms) {\n")
edit_file("tests/testthat/test-latres-brms.R",
          paste0("      # a chain at the numeric edge of its range has no Cholesky factor in\n",
                 "      # brms's coordinates, so its program has no point to compare at\n",
                 "      r <- tryCatch(\n",
                 "        brms_lp_check(brms::bf(f), fm[[1]], d, fit, joint = TRUE),\n",
                 "        error = function(e) {\n",
                 "          if (!grepl(\"leading minor\", conditionMessage(e), fixed = TRUE)) stop(e)\n",
                 "          NULL\n",
                 "        })\n",
                 "      if (is.null(r)) {\n",
                 "        at_edge <- at_edge + 1L\n",
                 "        next\n",
                 "      }\n"),
          paste0("      r <- brms_lp_check(brms::bf(f), fm[[1]], d, fit, joint = TRUE)\n"))
edit_file("tests/testthat/test-latres-brms.R",
          paste0("      expect_lt(abs(r$measured_const), 1e-9 * max(1, abs(r$ours)))\n"),
          paste0("      # unstr() on Beta() walks to a nearly singular correlation (6e6 on the chart),\n",
                 "      # where rounding reads -3.3e-8 on the reference BLAS and 8.6e-8 on OpenBLAS\n",
                 "      edge <- grepl(\"unstr\", tm, fixed = TRUE) && inherits(fm[[2]], \"frmtmb_family\") &&\n",
                 "        identical(fm[[2]][[\"family\"]], \"beta\")\n",
                 "      tol <- if (edge) 1e-8 else 1e-9\n",
                 "      expect_lt(abs(r$measured_const), tol * max(1, abs(r$ours)))\n"))
edit_file("tests/testthat/test-latres-brms.R",
          paste0("  expect_lte(at_edge, 2L)\n"), "")
# RF1's tests, by class
tf <- "tests/testthat/test-latres.R"
edit_file(tf,
          paste0("  # three trials has a row's share of a poisson of mean 1 that reads true, so\n",
                 "  # the condition reads the trials\n",
                 "  b3 <- lr_series(3L, 15L, 10L, 0.2, 0.8, tr = 3L, y = \"bin\")\n",
                 "  said <- lr_laplace_warned(frm(bf(ft), family = binomial(), data = b3))\n",
                 "  expect_length(said, 1L)\n",
                 "  expect_match(said, \"3 trials\", fixed = TRUE)\n"),
          paste0("  # above the row limit: three trials has a row's share of a poisson of mean 1\n",
                 "  # that reads true, so the condition reads the trials\n",
                 "  b3 <- lr_series(3L, 15L, 10L, 0.2, 0.8, tr = 3L, y = \"bin\")\n",
                 "  said <- withr::with_options(list(frmtmb.latres_se_rows = 10L),\n",
                 "    lr_laplace_warned(frm(bf(ft), family = binomial(), data = b3)))\n",
                 "  expect_length(said, 1L)\n",
                 "  expect_match(said, \"3 trials\", fixed = TRUE)\n"))
# the stand-aside's mocks read the row limit's branch
edit_file(tf,
          paste0("  # the stand-aside by construction: a thin row with no third derivative reads no\n",
                 "  # reason, and the same row with one does, so no other condition answers for it\n"),
          paste0("  # the stand-aside by construction, above the row limit: a thin row with no\n",
                 "  # third derivative reads no reason, and the same row with one does, so no\n",
                 "  # other condition answers for it\n",
                 "  withr::local_options(frmtmb.latres_se_rows = 10L)\n"))
# this design's skew reads 3.45 with a standard error of 2.8, so the shift is 2.16,
# on the lane's build as well: what the test holds is that units do not move it
edit_file(tf,
          paste0("  d$y <- e + z\n",
                 "  expect_equal(lr_warned_on(f, skew_normal(), d), 0L)\n",
                 "  d$y <- 10 * (e + z)\n",
                 "  expect_equal(lr_warned_on(f, skew_normal(), d), 0L)\n"),
          paste0("  d$y <- e + z\n",
                 "  w1 <- lr_warned_on(f, skew_normal(), d)\n",
                 "  d$y <- 10 * (e + z)\n",
                 "  expect_equal(lr_warned_on(f, skew_normal(), d), w1)\n"))
# negbinomial is the row rule's case above the limit, and below it the gate reads
# 16 of 38 where the fits are wrong on 16 (the review's `nbg.log`)
edit_file(tf,
          "  expect_equal(lr_warned_on(f, negbinomial(), d), 1L)\n",
          paste0("  expect_equal(withr::with_options(list(frmtmb.latres_se_rows = 10L),\n",
                 "                                  lr_warned_on(f, negbinomial(), d)), 1L)\n"))
# a chain of a cosy() beside ar() on a shuffled second response is read too
edit_file(tf,
          paste0("    frm(bf(cnt ~ x + ar(t, g, cov = TRUE)) + bf(cnt2 ~ x + cosy(t, g)),\n",
                 "        family = poisson(), data = d),\n",
                 "    c(\"Optimizer did not report convergence\",\n"),
          paste0("    frm(bf(cnt ~ x + ar(t, g, cov = TRUE)) + bf(cnt2 ~ x + cosy(t, g)),\n",
                 "        family = poisson(), data = d),\n",
                 "    c(\"The Laplace approximation\",\n",
                 "      \"Optimizer did not report convergence\",\n"))
# S4c: M3's rule. Est.Error holds the draw beside the coefficients'
edit_file(tf,
          paste0("  expect_equal(unname(v2[1L, \"Est.Error\"]),\n",
                 "               sqrt(drop(crossprod(g, V %*% g))), tolerance = 0.05)\n"),
          paste0("  draw <- m * sqrt(exp(vv) - 1)\n",
                 "  got <- unname(v2[1L, \"Est.Error\"])\n",
                 "  expect_gte(got, 0.999 * draw)\n",
                 "  expect_lt(got, 1.05 * sqrt(draw^2 + drop(crossprod(g, V %*% g))))\n"))
# RF4's other pin
edit_file("tests/testthat/test-refuse-pairs.R",
          paste0("  expect_error(f(bf(y ~ x + (1 | p | gr(g, cor = FALSE)),\n",
                 "                    sigma ~ (1 | p | g))),\n",
                 "               \"will not pick one silently\", fixed = TRUE)\n"),
          paste0("  # the |p| spelling drops its key under cor = FALSE, as brms's program does\n",
                 "  fr <- f(bf(y ~ x + (1 | p | gr(g, cor = FALSE)), sigma ~ (1 | p | g)))\n",
                 "  expect_length(fr$re_blocks, 2L)\n"))
# RF5: the help of frm_lp_basis() states the cost
edit_file("R/predict.R",
          "#'     contributes), a dense matrix otherwise.}\n",
          paste0("#'     contributes), a dense matrix otherwise. Beside a latent residual, `A` and\n",
                 "#'     `V` are dense, about `16 n^2` bytes for `n` rows: 1 GB at 8000 rows.}\n"))

# ---- gpgap ---------------------------------------------------------------
# RF1 holds on ncgp's gp_factor_ordered(); its test row is in test-gp-noncentered.R
# (below). RF2, RF3 were done on ncgp's tree. RF4: the factored route's clock
# is ncgp's sentence, and RF5 is in the page written above.

# ---- nlfix ---------------------------------------------------------------
edit_file("extensions/frmtmb.learn/R/family.R",
          "counterfactual = NULL, session_expr = NULL) {",
          "counterfactual = NULL, session_expr = NULL,\n                      type = \"discrete\") {")
edit_file("extensions/frmtmb.learn/R/rlddm.R",
          "            sim = FALSE,\n",
          "            sim = FALSE, type = \"continuous\",\n")
x <- rd("extensions/frmtmb.learn/R/family.R")
k <- grep("type = \"discrete\",?$", x)
k <- k[k > grep("^ln_family <- function", x)]
k <- k[!grepl("^ +type = \"discrete\") \\{", x[k])]
print(x[k])
stopifnot(length(k) == 1L)
x[k] <- sub("type = \"discrete\"", "type = type", x[k], fixed = TRUE)
wr(x, "extensions/frmtmb.learn/R/family.R")
cat("edited extensions/frmtmb.learn/R/family.R: type\n")
# RF2: a first agreement stops only an exact alias, where a restart cannot gain
edit_file("R/par-template.R",
          paste0("    if (from_flat && same(ll, floor)) break\n",
                 "    agrees <- any(vapply(seen, function(s) {\n",
                 "      same(ll, s$ll) && same_par(r$opt$par, s$par)\n",
                 "    }, NA))\n",
                 "    if (agrees && !(saddle_floor && same(ll, floor))) break\n",
                 "    seen[[length(seen) + 1L]] <- list(ll = ll, par = r$opt$par)\n"),
          "    if (from_flat && same(ll, floor)) break\n")
edit_file("R/par-template.R",
          paste0("  same_par <- function(a, b) {\n",
                 "    length(a) == length(b) && all(abs(a - b) <= 1e-3 * (1 + abs(a)))\n",
                 "  }\n"),
          "")
edit_file("R/par-template.R", "  best_held <- list()\n  seen <- list()\n", "  best_held <- list()\n")
edit_file("R/par-template.R",
          "                              from_flat = !on_saddle)\n",
          paste0("                              from_flat = !on_saddle && any(vapply(held, function(w) {\n",
                 "                                grepl(\"are not identified\", conditionMessage(w))\n",
                 "                              }, NA)))\n"))
tf <- "tests/testthat/test-nl-family.R"
x <- rd(tf)
a <- grep("^test_that\\(\"restarts that agree end the round\"", x)
b <- a + which(x[(a + 1L):length(x)] == "})")[1L]
stopifnot(length(a) == 1L)
new <- c(
  "test_that(\"a restart past an agreement finds the flip of the sign\", {",
  "  # acat, seed 4, slope -0.6: a stop at the first agreement read -376.2397 at b1",
  "  # -15200, the ridge, at code 0, and the patterns past it read -370.2536",
  "  # (dev/reviews/2026-11-02-nlfix.md, release fix 2; dev/review-nlfix-d7.R)",
  "  set.seed(401)",
  "  n <- 300",
  "  d <- data.frame(x = round(stats::rnorm(n), 3), z = round(stats::rnorm(n), 3))",
  "  lat <- 1.5 * exp(-0.6 * d$x) + 0.3 * d$z",
  "  d$y <- ordered(cut(lat + stats::rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))",
  "  body <- bf(y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ 1, nl = TRUE)",
  "  fit <- suppressMessages(suppressWarnings(frm(body, family = acat(), data = d)))",
  "  best <- -Inf",
  "  for (b2 in c(-0.5, 0.5)) for (b1 in c(-1, 1)) {",
  "    r <- suppressMessages(suppressWarnings(frm(",
  "      body, family = acat(), data = d,",
  "      start = list(beta = c(b1_Intercept = b1, b2_Intercept = b2)))))",
  "    best <- max(best, as.numeric(logLik(r)))",
  "  }",
  "  # 6 nats separate the ridge from the flip, on a scale the fit reads to 1e-6",
  "  expect_lt(best - as.numeric(logLik(fit)), 1e-5 * abs(best))",
  "})")
x <- c(x[seq_len(a - 1L)], new, x[-seq_len(b)])
wr(x, tf)
cat("edited", tf, "\n")
# NIT: both sentences
edit_file("R/fit.R",
          paste0("#'   marks a second one. A `start` with other signs can show it, and\n",
                 "#'   [frm_allfit()] cannot, as it starts every optimizer from one point.\n"),
          paste0("#'   marks a second one. A `start` with other signs, or a located prior, can\n",
                 "#'   show it, and [frm_allfit()] cannot, as it starts every optimizer from one\n",
                 "#'   point.\n"))
edit_file("R/par-template.R",
          "#' another sign of `b1` can show it. [frm_allfit()] cannot, as its optimizers\n",
          "#' another sign of `b1`, or a located prior, can show it. [frm_allfit()] cannot, as its optimizers\n")

# ---- predread ------------------------------------------------------------
# F1 to F4 of dev/review-predread-r2-release.R, each pattern matching once. F1's
# `at` goes as the review's text asks; the script's `at <- NULL` is the same path.
pr_ed <- function(f, old, new) edit_file(f, gsub("\r", "", old, fixed = TRUE),
                                         gsub("\r", "", new, fixed = TRUE))
x <- rd("R/predict.R")
a <- which(x == "  at <- if (!is.null(mm)) matrix(mm$plug, dim(out)[1L], dim(out)[3L])")
stopifnot(length(a) == 1L)
cm <- a
while (grepl("^  #", x[cm - 1L])) cm <- cm - 1L
cat("F1 dropped:\n", paste0(x[cm:a], "\n"), sep = "")
x <- c(x[seq_len(cm - 1L)],
       "  if (!is.null(mm)) return(prob_band_draw(out, probs, mm, up))",
       x[-seq_len(a)])
wr(x, "R/predict.R")
pr_ed("R/predict.R",
      "pk <- if (is.null(at)) mk else as.numeric(at[, k]) / up",
      "pk <- mk")
pr_ed("R/predict.R",
      "prob_band_quantiles <- function(out, probs, mm = NULL, upper = 1) {",
      paste0(
        "prob_band_draw <- function(out, probs, mm, up) {\n",
        "  q <- setdiff(seq_len(dim(out)[2L]), 1:2)\n",
        "  # a middle category has a hump in the draw, and no interval about one point\n",
        "  # reads that, so each row reads the draw's own quantile, on 99 points of\n",
        "  # equal weight; the logit carries the coefficients' share of the variance\n",
        "  vals <- lapply(stats::qnorm(stats::ppoints(99L)), mm$val)\n",
        "  vr <- matrix(mm$var, dim(out)[1L], dim(out)[3L])\n",
        "  for (k in seq_len(dim(out)[3L])) {\n",
        "    V <- vapply(vals, function(v) as.matrix(v)[, k], numeric(dim(out)[1L]))\n",
        "    V <- matrix(V, dim(out)[1L], length(vals))\n",
        "    # each row is sorted once and read at its two ends, where quantile() sorts per end\n",
        "    V <- matrix(t(apply(V, 1L, sort, na.last = TRUE)), nrow(V), ncol(V))\n",
        "    m <- as.numeric(out[, 1L, k]) / up\n",
        "    se <- as.numeric(out[, 2L, k]) / up\n",
        "    v2 <- pmax(vr[, k] / up^2, .Machine$double.eps)\n",
        "    infl <- sqrt(1 + pmax(se^2 - v2, 0) / v2)\n",
        "    inside <- m > 0 & m < 1\n",
        "    for (i in seq_along(q)) {\n",
        "      h <- (ncol(V) - 1L) * probs[i] + 1\n",
        "      v <- (V[, floor(h)] + (h - floor(h)) * (V[, ceiling(h)] - V[, floor(h)])) / up\n",
        "      r <- stats::plogis(\n",
        "        stats::qlogis(m) + (stats::qlogis(v) - stats::qlogis(m)) * infl)\n",
        "      # a probability at 0 or 1 has no logit, and no spread either\n",
        "      r[which(!inside)] <- m[which(!inside)]\n",
        "      out[, q[i], k] <- r * up\n",
        "    }\n",
        "  }\n",
        "  out\n",
        "}\n\n",
        "prob_band_quantiles <- function(out, probs, mm = NULL, upper = 1) {"))
pr_ed("R/curve-cov.R",
      "    if (!length(E@x) || !any(E@x != 0)) return(S)\n",
      "    if (!length(E@x) || !any(E@x != 0, na.rm = TRUE)) return(S)\n")
pr_ed("R/curve-cov.R",
      "    if (!length(E@x) || !any(E@x != 0)) return(NULL)\n",
      "    if (!length(E@x) || !any(E@x != 0, na.rm = TRUE)) return(NULL)\n")
pr_ed("R/curve-cov.R",
      "  if (!any(E != 0)) NULL else unname(E)",
      "  if (!any(E != 0, na.rm = TRUE)) NULL else unname(E)")
pr_ed("R/curve-cov.R",
      "  ev <- eigen((S + t(S)) / 2, symmetric = TRUE)\n  L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow = m)",
      paste0(
        "  # a row of NA in the grid has no place in the band; the others do\n",
        "  gone <- !is.finite(diag(S))\n",
        "  if (any(gone)) {\n",
        "    S <- S[!gone, !gone, drop = FALSE]\n",
        "    keep <- keep[!gone]\n",
        "    div <- div[!gone]\n",
        "    m <- nrow(S)\n",
        "  }\n",
        "  ev <- eigen((S + t(S)) / 2, symmetric = TRUE)\n  L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow = m)"))
pr_ed("R/interop.R",
      paste0(
        "  frm_warning(\"marginaleffects forms its standard errors from the population-\",\n",
        "              \"level coefficients alone, so the standard error of a level \",\n",
        "              \"from predictions() holds no variance of \",\n",
        "              paste(omits, collapse = \", or of \"),\n",
        "              \". A slope, contrast or test that reads those coefficients \",\n",
        "              \"alone is not affected. frm_linpred(newdata = , se.fit = TRUE) \",\n",
        "              \"and fitted(newdata = ) carry it. vcov = FALSE gives the \",\n",
        "              \"estimates alone, with no standard error\", call. = FALSE,\n"),
      paste0(
        "  frm_warning(\"marginaleffects forms its standard errors from the population-\",\n",
        "              \"level coefficients alone. It holds fixed \",\n",
        "              paste(omits, collapse = \", and \"),\n",
        "              \". The standard error of a level from predictions(), or of an \",\n",
        "              \"average of levels, can be too wide or too narrow, often by a \",\n",
        "              \"factor of 2 or more. On a log or logit link, or with a \",\n",
        "              \"group-level slope, a slope or contrast can be off too. A test \",\n",
        "              \"of a coefficient is not affected. frm_linpred(newdata = , \",\n",
        "              \"se.fit = TRUE) and fitted(newdata = ) give a level's standard \",\n",
        "              \"error. vcov = FALSE gives the estimates alone, with no \",\n",
        "              \"standard error\", call. = FALSE,\n"))
pr_ed("R/interop.R",
      "\"the group-level effects (re_formula = NA drops them, and then it agrees)\"",
      "\"the group-level effects (re_formula = NA drops them, and then a level agrees)\"")
pr_ed("R/compat.R",
      "every row reads as new, and get_vcov() gets no rows. A slope, contrast or test of a coefficient is right. frm_linpred",
      "every row reads as new, and get_vcov() gets no rows. A test of a coefficient is right, and so is a slope or contrast on an identity link. frm_linpred")
pr_ed("R/me.R",
      paste0(
        "  memo <- new.env(parent = emptyenv())\n",
        "  function(z, d = 0) {\n",
        "    key <- paste(c(z, d), collapse = \" \")\n",
        "    old <- memo[[key]]\n",
        "    if (!is.null(old)) return(old)\n"),
      paste0(
        "  # a list: z has one value per row, and an environment's name holds 10000\n",
        "  # bytes, which about 500 rows pass. A fit reads a few hundred nodes at most\n",
        "  memo <- list()\n",
        "  function(z, d = 0) {\n",
        "    for (m in memo) if (m$d == d && identical(m$z, z)) return(m$value)\n"))
pr_ed("R/me.R",
      paste0(
        "    memo[[key]] <- what(sh)\n",
        "    memo[[key]]\n"),
      paste0(
        "    value <- what(sh)\n",
        "    memo[[length(memo) + 1L]] <<- list(z = z, d = d, value = value)\n",
        "    value\n"))
x <- rd("R/predict.R")
a <- grep("The `Q` columns of a category probability", x)
stopifnot(length(a) == 1L)
print(x[a:(a + 7L)])

# ---- ncgp ----------------------------------------------------------------
# R1: the warm start's first read takes six length scales; gp_marginal_opt()'s
# own default of 4 is what the working start reads, and at 6 it loses there.
edit_file("R/gp-start.R",
          "gp_marginal_max <- function(frame, keep = 3L, distinct = 4L) {",
          "gp_marginal_max <- function(frame, keep = 3L, distinct = 6L) {")
edit_file("tests/testthat/test-gp-start.R",
          "local_mocked_bindings(gp_marginal_max = function(frame, keep = 3L, distinct = 4L) {",
          "local_mocked_bindings(gp_marginal_max = function(frame, keep = 3L, distinct = 6L) {")
edit_file("R/fit.R",
          paste0("#'   the maximum of the model's marginal, read from a coarse grid. Another\n"),
          paste0("#'   the maximum of the model's marginal, read from the best cell of six\n",
                 "#'   length scales on a coarse grid. Another\n"))
edit_file("R/fit.R",
          paste0("#'   384 rows of 16 curves no test had seen, the fit stopped more than 0.001\n",
                 "#'   short of the maximum on 5, the worst 9.4, with no warning. A fit from other\n"),
          paste0("#'   744 rows of curves no test had seen, the fit stopped more than 0.001\n",
                 "#'   short of the maximum on 1, the worst 0.56, with no warning, where 0.78.0\n",
                 "#'   read 40 and 62 short on the two tables of 192. A fit from other\n"))
edit_file("R/gp-start.R",
          "#' `distinct` different length scales and from the best `keep` cells overall, the\n",
          paste0("#' `distinct` different length scales (six for the warm start, 1 of 744 rows\n",
                 "#' short, and four for the working start, where six lose on low-count Poisson and\n",
                 "#' binomial) and from the best `keep` cells overall, the\n"))
# R3
x <- rd("extensions/frmtmb.sample/R/sample.R")
a <- grep("^#' small `n`: at n 15, Rhat on `sdgp` read 1\\.04 and 1\\.69 at `1e-3`, 1\\.36 at$", x)
stopifnot(length(a) == 1L)
cat("R3 old:\n", paste0(x[a:(a + 2L)], "\n"), sep = "")
ol <- x[a + 1L]
if (!grepl("^#' `1e-4` and", ol)) stop("R3 second line: ", ol)
x[a] <- "#' small `n`: at n 15, Rhat on `sdgp` read 1.04 to 2.37 at `1e-3`, 1.36 to 1.95 at"
x[a + 1L] <- "#' `1e-4` and 1.56 to 2.47 at `1e-6`, and `sigma`'s median up to 394 times the"
cat("R3 new line 2 replaced:", ol, "\n")
wr(x, "extensions/frmtmb.sample/R/sample.R")

# ---- cross-lane tests --------------------------------------------------
source("dev/rel079-tests.R", local = FALSE)
cat("RELEASE FIXES DONE\n")
