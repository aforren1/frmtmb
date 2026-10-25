# Added at the 0.75.0 consolidation (dev/rel075-resolve.R, item 6): the
# combined refusal above and the residual correlation seam
# (frmtmb_register_autocor(), sar() in frmtmb.spatial) meet only in the
# merged tree. One model cannot hold cox() and sar() (sar() takes a
# gaussian or student response, and cox() is refused in a multivariate
# fit), so the family, the structure and the residual term come from
# three absent packages in process, and a real car() and sar() fit is
# read in child processes.

mx_absent_c <- "frmtmbabsenttestpkgc"

# a fixed-matrix residual term, as fcor(), owned by `pkg`
mx_register_fix <- function(pkg) {
  frmtmb_register_autocor(
    "mxfix", args = "V",
    build = function(spec, rows, n_all, lookup) {
      V <- lookup(spec$V, "mxfix()")[rows, rows, drop = FALSE]
      L <- t(chol(V))
      list(L = L, Linv = forwardsolve(L, diag(nrow(L))),
           logdet = sum(log(diag(L))))
    },
    whiten = function(y, mu, sigma, theta, aux) {
      n <- length(y)
      lsig <- if (length(sigma) == 1L) n * log(sigma) else sum(log(sigma))
      list(u = as.vector(aux$Linv %*% ((y - mu) / sigma)),
           logjac = -aux$logdet - lsig)
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      mu + sigma * (aux$L %*% u)
    },
    package = pkg, replace = TRUE)
}

# a fit with mxiid() (package `struct_pkg`; left out when `struct` is
# FALSE), mxfix() (package `term_pkg`) and a family said to be
# `fam_pkg`'s; the registrations are undone on return, as in a session
# that reads a saved fit
mx_fit3 <- function(fam_pkg, struct_pkg, term_pkg, struct = TRUE) {
  local_frmtmb_registries()
  do.call(frmtmb_register_covstruct,
          c(list(name = "mxiid"), mx_iid_slots(),
            list(package = struct_pkg, replace = TRUE)))
  mx_register_fix(term_pkg)
  set.seed(2)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = stats::rnorm(100))
  V <- 0.5^abs(outer(1:100, 1:100, "-"))
  d$y <- 1 + d$x + stats::rnorm(20)[d$g] +
    as.vector(t(chol(V)) %*% stats::rnorm(100))
  fo <- if (struct) y ~ x + mxiid(1 | g) + mxfix(V) else y ~ x + mxfix(V)
  fit <- frm(bf(fo), family = gaussian(), data = d, data2 = list(V = V))
  fit$spec$responses[[1L]]$family <-
    mx_owned_by(fit$spec$responses[[1L]]$family, fam_pkg)
  fit$frame$spec <- fit$spec
  fit$cache <- new.env()
  fit
}

test_that("a family, a structure and a residual term from three absent
           packages are named in one refusal", {
  three <- mx_fit3(mx_absent_a, mx_absent_b, mx_absent_c)
  expect_false(any(c("mxiid", "mxfix") %in%
                     c(cs_names(), frmtmb:::ac_entry_names())))
  msg3 <- paste0(
    "This model needs 3 extension packages that are not installed or ",
    "cannot be loaded in this session, so their functions are not ",
    "available: ", mx_absent_a, " (the gaussian family of this model), ",
    mx_absent_b, " (the mxiid() structure of this model) and ",
    mx_absent_c, " (the mxfix() term of this model). Install them and ",
    "run library(", mx_absent_a, "); library(", mx_absent_b,
    "); library(", mx_absent_c, ")")
  for (nm in names(mx_calls)) {
    e <- mx_error(three, nm)
    expect_s3_class(e, "frmtmb_extension_missing")
    expect_identical(conditionMessage(e %||% simpleError("")), msg3,
                     info = nm)
  }
  # the term alone: its package in the one-package words, whichever way
  # the method reaches it (the fit-level check, or the block's entry)
  term_only <- mx_fit3("frmtmb", NULL, mx_absent_c, struct = FALSE)
  e <- mx_error(term_only, "summary")
  expect_s3_class(e, "frmtmb_extension_missing")
  expect_identical(conditionMessage(e %||% simpleError("")),
                   mx_one("The mxfix() term of this model", mx_absent_c))
  ac <- term_only$frame$autocor[[1L]]
  e2 <- tryCatch(frmtmb:::ac_block_entry(ac), error = identity)
  expect_s3_class(e2, "frmtmb_extension_missing")
  expect_identical(conditionMessage(e2),
                   mx_one("The mxfix() term of this model", mx_absent_c))
  # a term registered in a session names its own seam
  sess <- mx_fit3("frmtmb", NULL, NULL, struct = FALSE)
  want_sess <- paste0(
    "The mxfix() term of this model was registered in an R session, not ",
    "by a package, and is not registered in this one. Register it again ",
    "with frmtmb_register_autocor() before using the fit")
  ac <- sess$frame$autocor[[1L]]
  e3 <- tryCatch(frmtmb:::ac_block_entry(ac), error = identity)
  expect_s3_class(e3, "frmtmb_extension_missing")
  expect_identical(conditionMessage(e3), want_sess)
  e4 <- mx_error(sess, "summary")
  expect_s3_class(e4, "frmtmb_extension_missing")
  expect_identical(conditionMessage(e4 %||% simpleError("")), want_sess)
  # the guard absent: with every part registered, the fit answers
  local_frmtmb_registries()
  ok <- mx_fit3("frmtmb", mx_absent_b, mx_absent_c)
  do.call(frmtmb_register_covstruct,
          c(list(name = "mxiid"), mx_iid_slots(),
            list(package = mx_absent_b, replace = TRUE)))
  mx_register_fix(mx_absent_c)
  for (nm in names(mx_calls)) expect_null(mx_error(ok, nm), info = nm)
})

test_that("a saved car() and sar() fit names frmtmb.spatial once", {
  skip_on_cran()
  if (!nzchar(system.file(package = "frmtmb.spatial"))) {
    skip(paste("the saved-fit arm needs these packages installed:",
               "frmtmb.spatial"))
  }
  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(rds), add = TRUE)
  made <- child_with_libs(.libPaths(), c(
    "suppressMessages({library(frmtmb); library(frmtmb.spatial)})",
    "g <- expand.grid(r = 1:3, c = 1:3)",
    "W <- 1 * (as.matrix(dist(g, method = 'manhattan')) == 1)",
    "dimnames(W) <- list(paste0('L', 1:9), paste0('L', 1:9))",
    "set.seed(20261024)",
    "dd <- data.frame(loc = factor(rep(rownames(W), each = 12),",
    "                              levels = rownames(W)), x = rnorm(108))",
    "Wr <- W[dd$loc, dd$loc]",
    "Wr <- Wr / rowSums(Wr)",
    "dd$y <- as.vector(solve(diag(108) - 0.4 * Wr,",
    "                        1 + 0.5 * dd$x + rnorm(108)))",
    "fit <- suppressWarnings(suppressMessages(frm(",
    "  bf(y ~ x + car(W, gr = loc, type = 'escar') + sar(Wr)),",
    "  family = gaussian(), data = dd, data2 = list(W = W, Wr = Wr))))",
    sprintf("saveRDS(fit, %s)", deparse(rds)),
    "cat('SAVED\\n')"))
  expect_true(any(made == "SAVED"), info = paste(made, collapse = "\n"))
  read_in <- function(libs) {
    child_with_libs(libs, c(
      "suppressMessages(library(frmtmb))",
      "cat('HAVE', requireNamespace('frmtmb.spatial', quietly = TRUE),",
      "    '\\n')",
      sprintf("f <- suppressWarnings(readRDS(%s))", deparse(rds)),
      "calls <- list(print = function() capture.output(print(f)),",
      "  summary = function() summary(f), logLik = function() logLik(f),",
      "  vcov = function() vcov(f), fitted = function() fitted(f),",
      "  residuals = function() residuals(f),",
      "  autocor_matrix = function() autocor_matrix(f),",
      "  simulate = function() simulate(f, nsim = 1, seed = 1),",
      "  update = function() update(f), fixef = function() fixef(f),",
      "  hypothesis = function() hypothesis(f, 'x > 0'))",
      "for (m in names(calls)) {",
      "  e <- tryCatch({ calls[[m]](); NULL }, error = function(e) e)",
      "  cat('CALL', m, '::', if (is.null(e)) 'OK' else",
      "      paste(inherits(e, 'frmtmb_extension_missing'), '::',",
      "            gsub('\\n', ' ', conditionMessage(e))), '\\n')",
      "}"))
  }
  want <- paste0(
    "The car() structure of this model and the sar() term of this model ",
    "are from the extension package frmtmb.spatial, which is not ",
    "installed or cannot be loaded in this session, so its functions are ",
    "not available. Install it and run library(frmtmb.spatial)")
  arms <- list(unloadable = function() {
    unloadable_package_libs("frmtmb.spatial")
  }, absent = function() absent_package_libs("frmtmb.spatial"))
  for (a in names(arms)) {
    out <- read_in(arms[[a]]())
    expect_true(any(out == "HAVE FALSE "), info = paste(out, collapse = "\n"))
    cl <- grep("^CALL ", out, value = TRUE)
    expect_length(cl, 11L)
    # every method but logLik(), which reads the stored optimum
    rf <- grep(":: TRUE ::", cl, value = TRUE, fixed = TRUE)
    expect_length(rf, 10L)
    for (r in rf) expect_match(r, want, fixed = TRUE, label = paste(a, r))
    expect_true(any(startsWith(cl, "CALL logLik :: OK")))
  }
})
