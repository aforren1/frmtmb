# A fit whose family is from one extension and whose structure is from
# another (backlog "A saved fit that needs two extensions", 0.74.0).
# Each check fires only for its own part, so a method that needs
# nothing of an absent package still answers, but a refusal names every
# part this session cannot use, so one install covers them all. Until
# then summary() of a saved cox() and car() fit named frmtmb.survival,
# print() named frmtmb.spatial, and the user met the second package only
# after installing the first.
#
# The first two blocks own the family and the structure by packages no
# session has. The third reads a real cox() and car() fit in child
# processes (helper-absent-package.R), and runs where both extensions
# are installed.

mx_absent_a <- "frmtmbabsenttestpkga"
mx_absent_b <- "frmtmbabsenttestpkgb"

mx_owned_by <- function(fam, pkg) {
  attr(fam, "frmtmb_package") <- pkg
  fam
}

mx_iid_slots <- function() {
  list(lpdf = function(b, theta, blk) {
         sum(RTMB::dnorm(b, 0, exp(theta[1]), log = TRUE))
       },
       vcov = function(theta, blk) {
         matrix(exp(2 * theta[1]), 1, 1,
                dimnames = list(blk[["cnms"]], blk[["cnms"]]))
       },
       npar = 1L, sd_idx = 1L)
}

# a fit with the structure mxiid(), registered by `struct_pkg`, whose
# family is then said to be `fam_pkg`'s, with a fresh cache as
# readRDS() gives. The registration is undone on return when `keep` is
# FALSE, as in a session that reads a saved fit.
mx_fit <- function(fam_pkg, struct_pkg, keep = FALSE,
                   env = parent.frame()) {
  if (keep) local_frmtmb_registries(env) else local_frmtmb_registries()
  do.call(frmtmb_register_covstruct,
          c(list(name = "mxiid"), mx_iid_slots(),
            list(package = struct_pkg, replace = TRUE)))
  set.seed(1)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = stats::rnorm(100))
  d$y <- 1 + d$x + stats::rnorm(20)[d$g] + stats::rnorm(100, 0, 0.5)
  fit <- frm(bf(y ~ x + mxiid(1 | g)), family = gaussian(), data = d)
  fit$spec$responses[[1L]]$family <-
    mx_owned_by(fit$spec$responses[[1L]]$family, fam_pkg)
  fit$frame$spec <- fit$spec
  fit$cache <- new.env()
  fit
}

mx_calls <- list(
  print = function(f) utils::capture.output(print(f)),
  summary = function(f) summary(f), vcov = function(f) vcov(f),
  VarCorr = function(f) VarCorr(f), ranef = function(f) ranef(f),
  predict = function(f) predict(f), update = function(f) update(f),
  simulate = function(f) simulate(f, nsim = 1),
  confint_varcorr = function(f) confint_varcorr(f),
  # these reached one part on their own, or answered: fixef() read the
  # refusal from vcov() as "no covariance" and gave Est.Error NA
  fixef = function(f) fixef(f),
  hypothesis = function(f) hypothesis(f, "x > 0"),
  residuals = function(f) residuals(f), drop1 = function(f) drop1(f),
  ngrps = function(f) ngrps(f), variables = function(f) variables(f),
  # emmeans runs recover_data() inside try() and replaced a refusal
  # raised there with "Perhaps a 'data' or 'params' argument is needed"
  # (the 0.75.0 release review, finding 1)
  ref_grid = function(f) emmeans::ref_grid(f),
  emmeans = function(f) emmeans::emmeans(f, ~ x))

mx_error <- function(f, nm) {
  tryCatch({ mx_calls[[nm]](f); NULL }, error = identity)
}

mx_one <- function(what, pkg) {
  paste0(what, " is from the extension package ", pkg, ", which is not ",
         "installed or cannot be loaded in this session, so its functions ",
         "are not available. Install it and run library(", pkg, ")")
}

test_that("a fit that needs two absent packages names both at once", {
  both <- mx_fit(mx_absent_a, mx_absent_b)
  struct_only <- mx_fit("frmtmb", mx_absent_b)
  expect_false("mxiid" %in% cs_names())
  msg2 <- paste0(
    "This model needs 2 extension packages that are not installed or ",
    "cannot be loaded in this session, so their functions are not ",
    "available: ", mx_absent_a, " (the gaussian family of this model) ",
    "and ", mx_absent_b, " (the mxiid() structure of this model). ",
    "Install them and run library(", mx_absent_a, "); library(",
    mx_absent_b, ")")
  for (nm in names(mx_calls)) {
    e <- mx_error(both, nm)
    expect_s3_class(e, "frmtmb_extension_missing")
    expect_identical(conditionMessage(e %||% simpleError("")), msg2,
                     info = nm)
    # one absent: that package alone, in the one-package words
    e1 <- mx_error(struct_only, nm)
    expect_s3_class(e1, "frmtmb_extension_missing")
    expect_identical(conditionMessage(e1 %||% simpleError("")),
                     mx_one("The mxiid() structure of this model",
                            mx_absent_b), info = nm)
  }
})

test_that("a refusal for the family alone names the family's package", {
  # the structure is registered in this session, so only the family's
  # package is missing; a method that needs nothing of it answers
  fam_only <- mx_fit(mx_absent_a, mx_absent_b, keep = TRUE)
  expect_true("mxiid" %in% cs_names())
  n_refused <- 0L
  for (nm in names(mx_calls)) {
    e <- mx_error(fam_only, nm)
    if (is.null(e)) next
    n_refused <- n_refused + 1L
    expect_s3_class(e, "frmtmb_extension_missing")
    expect_identical(conditionMessage(e),
                     mx_one("The gaussian family of this model", mx_absent_a),
                     info = nm)
  }
  expect_gt(n_refused, 0L)
  # the guard absent: neither package is missing, and every method
  # answers
  ok <- mx_fit("frmtmb", mx_absent_b, keep = TRUE)
  for (nm in names(mx_calls)) {
    expect_null(mx_error(ok, nm), info = nm)
  }
})

test_that("a family an older frmtmb built names the other absent packages", {
  # a core family under a moved family's name stands for a cox() fit
  # saved by frmtmb 0.73.0 (test-extension-home.R), here in a fit whose
  # block needs a package this session lacks. The instruction installs
  # every package that is not installed and refits with both loaded
  old <- mx_owned_by(as_frmtmb_family(gaussian()), "frmtmb")
  old[["family"]] <- "cox"
  x <- list(spec = list(responses = list(list(family = old))),
            frame = list(re_blocks = list(
              list(covstruct = "mxabsent", covstruct_package = mx_absent_b))))
  surv <- nzchar(system.file(package = "frmtmb.survival"))
  want <- paste0(
    "The cox family of this model was built by an older frmtmb, when the ",
    "family was in frmtmb. It is in the extension package frmtmb.survival ",
    "now, and the functions this fit calls are no longer in frmtmb. ",
    "The mxabsent() structure of this model is from the extension ",
    "package ", mx_absent_b, ", which is not installed or cannot be ",
    "loaded in this session, so its functions are not available. ",
    "Install ", if (surv) mx_absent_b else
      paste0("frmtmb.survival and ", mx_absent_b),
    " and refit the model with library(frmtmb.survival); library(",
    mx_absent_b, ") loaded")
  for (f in list(require_model_packages, require_structure_packages)) {
    e <- tryCatch(f(x), error = identity)
    expect_s3_class(e, "frmtmb_extension_missing")
    expect_identical(conditionMessage(e), want)
  }
  # the guard absent: a block of core's own structure passes, and the
  # family's refusal then names no structure
  x$frame$re_blocks[[1L]] <- list(covstruct = "us")
  expect_true(require_structure_packages(x))
  e <- tryCatch(require_model_packages(x), error = identity)
  expect_s3_class(e, "frmtmb_extension_missing")
  expect_no_match(conditionMessage(e), "structure", fixed = TRUE)
})

test_that("a saved cox() and car() fit names every package it lacks", {
  skip_on_cran()
  need <- c("frmtmb.survival", "frmtmb.spatial")
  have <- vapply(need, function(p) nzchar(system.file(package = p)), NA)
  if (!all(have)) {
    # R CMD check hides undeclared packages, and core declares neither
    skip(paste("the saved-fit arm needs these packages installed:",
               paste(need[!have], collapse = ", ")))
  }
  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(rds), add = TRUE)
  # saved in a child, so that this session loads neither extension
  made <- child_with_libs(.libPaths(), c(
    "suppressMessages({library(frmtmb); library(frmtmb.survival)",
    "                  library(frmtmb.spatial)})",
    "g <- expand.grid(r = 1:3, c = 1:3)",
    "n <- nrow(g)",
    "W <- 1 * (as.matrix(dist(g, method = 'manhattan')) == 1)",
    "dimnames(W) <- list(paste0('L', 1:n), paste0('L', 1:n))",
    "set.seed(20261024)",
    "dd <- data.frame(loc = factor(rep(rownames(W), each = 12),",
    "                              levels = rownames(W)), x = rnorm(12 * n))",
    "u <- rnorm(n, 0, 0.5)",
    "dd$time <- rexp(nrow(dd), exp(-0.5 + 0.6 * dd$x + u[dd$loc]))",
    "dd$out <- rbinom(nrow(dd), 1, 0.25)",
    "fit <- suppressWarnings(suppressMessages(frm(",
    "  bf(time | cens(out) ~ x + car(W, gr = loc, type = 'escar')),",
    "  family = cox(df = 4), data = dd, data2 = list(W = W))))",
    sprintf("saveRDS(list(fit = fit, nd = dd[1:3, ]), %s)", deparse(rds)),
    "cat('SAVED\\n')"))
  expect_true(any(made == "SAVED"), info = paste(made, collapse = "\n"))
  read_in <- function(libs) {
    out <- child_with_libs(libs, c(
      "suppressMessages(library(frmtmb))",
      "cat('HAVE', vapply(c('frmtmb.survival', 'frmtmb.spatial'),",
      "    requireNamespace, NA, quietly = TRUE), '\\n')",
      sprintf("s <- suppressWarnings(readRDS(%s))", deparse(rds)),
      "f <- s$fit",
      "calls <- list(",
      "  print = function() capture.output(print(f)),",
      "  summary = function() summary(f), logLik = function() logLik(f),",
      "  vcov = function() vcov(f), VarCorr = function() VarCorr(f),",
      "  ranef = function() ranef(f), fitted = function() fitted(f),",
      "  predict_new = function() frm_linpred(f, newdata = s$nd,",
      "                                       se.fit = TRUE),",
      "  simulate = function() simulate(f, nsim = 1, seed = 1),",
      "  confint_varcorr = function() confint_varcorr(f),",
      "  update = function() update(f, . ~ . - x),",
      "  frm_curve = function() frm_curve(f,",
      "    newdata = data.frame(x = c(-1, 0, 1))),",
      "  fixef = function() fixef(f),",
      "  hypothesis = function() hypothesis(f, 'x > 0'),",
      "  residuals = function() residuals(f),",
      "  drop1 = function() drop1(f))",
      "for (m in names(calls)) {",
      "  e <- tryCatch({ calls[[m]](); NULL }, error = function(e) e)",
      "  cat('CALL', m, '::', if (is.null(e)) 'OK' else",
      "      paste(inherits(e, 'frmtmb_extension_missing'), '::',",
      "            gsub('\\n', ' ', conditionMessage(e))), '\\n')",
      "}"))
    list(have = grep("^HAVE ", out, value = TRUE),
         calls = grep("^CALL ", out, value = TRUE),
         all = paste(out, collapse = "\n"))
  }
  refused <- function(r) {
    grep(":: TRUE ::", r$calls, value = TRUE, fixed = TRUE)
  }
  pair <- paste0("frmtmb.survival (the cox family of this model) and ",
                 "frmtmb.spatial (the car() structure of this model)")
  # Both absent: every refusal names both, with one instruction.
  # Installed but unloadable first, then not installed (which skips
  # where frmtmb cannot load without the libraries that hold them)
  arms <- list(unloadable = function() unloadable_package_libs(need),
               absent = function() absent_package_libs(need))
  for (a in names(arms)) {
    r <- read_in(arms[[a]]())
    expect_identical(r$have, "HAVE FALSE FALSE ", info = r$all)
    expect_length(r$calls, 16L)
    rf <- refused(r)
    # every method but logLik(), which reads the stored optimum
    expect_length(rf, 15L)
    for (cl in rf) {
      expect_match(cl, pair, fixed = TRUE, label = paste(a, cl))
      expect_match(cl, paste0("Install them and run ",
                              "library(frmtmb.survival); ",
                              "library(frmtmb.spatial)"),
                   fixed = TRUE, label = paste(a, cl))
    }
    expect_false(any(grepl("could not find function", r$calls,
                           fixed = TRUE)), info = r$all)
  }
  # One absent: a refusal names that package alone. The other stays
  # installed and loadable, so this arm makes the one unloadable: a
  # stack without it alone cannot be built where both share a library
  for (p in need) {
    r <- read_in(unloadable_package_libs(p))
    expect_identical(r$have, if (identical(p, need[1L])) {
      "HAVE FALSE TRUE "
    } else "HAVE TRUE FALSE ", info = r$all)
    rf <- refused(r)
    expect_gt(length(rf), 0L)
    for (cl in rf) {
      expect_match(cl, paste0("is from the extension package ", p,
                              ", which is not installed"),
                   fixed = TRUE, label = cl)
      expect_no_match(cl, setdiff(need, p), fixed = TRUE)
    }
    expect_false(any(grepl("could not find function", r$calls,
                           fixed = TRUE)), info = r$all)
  }
})

test_that("fixef() refuses a fit with a missing part, not NA errors", {
  # fixef() read the refusal from vcov() as "no covariance" and answered
  # with Est.Error NA, in 0.74.0 too. A cox() fit saved by frmtmb 0.73.0
  # stands here as a core family under the moved name (as in
  # test-extension-home.R): it is refused whether or not
  # frmtmb.survival is installed
  set.seed(4)
  d <- data.frame(x = stats::runif(80))
  d$y <- 1 + 2 * d$x + stats::rnorm(80, 0, 0.3)
  fit <- frm(y ~ x, family = gaussian(), data = d)
  old <- fit
  fam <- mx_owned_by(old$spec$responses[[1L]]$family, "frmtmb")
  fam[["family"]] <- "cox"
  old$spec$responses[[1L]]$family <- fam
  old$frame$spec <- old$spec
  old$cache <- new.env()
  for (nm in c("fixef", "vcov", "summary", "hypothesis")) {
    e <- tryCatch({ mx_calls[[nm]](old); NULL }, error = identity)
    expect_s3_class(e, "frmtmb_extension_missing")
    expect_match(conditionMessage(e %||% simpleError("")),
                 "The cox family of this model was built by an older frmtmb",
                 fixed = TRUE, info = nm)
  }
  # the guard absent: the fit as it was has its standard errors
  expect_true(all(is.finite(fixef(fit)[, "Est.Error"])))
})

test_that("the instruction names the package after another sentence", {
  # review finding 3: with a session structure's sentence between the
  # missing package's sentence and the instruction, "Install it" read as
  # the structure
  f <- mx_fit(mx_absent_a, NULL)
  e <- mx_error(f, "summary")
  expect_s3_class(e, "frmtmb_extension_missing")
  expect_identical(conditionMessage(e), paste0(
    mx_one("The gaussian family of this model", mx_absent_a) |>
      sub(pattern = " Install it and run library[(].*$", replacement = ""),
    " The mxiid() structure of this model was registered in an R session, ",
    "not by a package, and is not registered in this one. Register it ",
    "again with frmtmb_register_covstruct() before using the fit. Install ",
    mx_absent_a, " and run library(", mx_absent_a, ")"))
})

test_that("an old family's home that cannot be loaded is to be reinstalled", {
  skip_on_cran()
  # review finding 4: the home of a family an older frmtmb built is not
  # loaded by the check, so an installed copy that cannot load was told
  # to "Refit the model with library(frmtmb.survival) loaded", which
  # fails. An installed package without Meta/package.rds is not whole
  code <- c(
    "suppressMessages(library(frmtmb))",
    "fam <- as_frmtmb_family(gaussian())",
    "attr(fam, 'frmtmb_package') <- 'frmtmb'",
    "fam[['family']] <- 'cox'",
    "e <- tryCatch(frmtmb:::require_model_packages(fam), error = identity)",
    "cat('MSG', gsub('\\n', ' ', conditionMessage(e)), '\\n')",
    "cat('LOADED', isNamespaceLoaded('frmtmb.survival'), '\\n')")
  msg <- function(out) sub("^MSG ", "", grep("^MSG ", out, value = TRUE))
  stub <- child_with_libs(unloadable_package_libs("frmtmb.survival"), code)
  expect_match(msg(stub), paste0(
    "frmtmb.survival is installed but cannot be loaded in this session. ",
    "Reinstall frmtmb.survival and refit the model with ",
    "library(frmtmb.survival) loaded"), fixed = TRUE,
    info = paste(stub, collapse = "\n"))
  expect_true(any(stub == "LOADED FALSE "))
  # the guard absent: a whole installed copy is loaded by the user
  if (nzchar(system.file(package = "frmtmb.survival"))) {
    whole <- child_with_libs(.libPaths(), code)
    expect_match(msg(whole), paste0(
      "no longer in frmtmb. Refit the model with ",
      "library(frmtmb.survival) loaded"), fixed = TRUE)
    expect_true(any(whole == "LOADED FALSE "))
  }
  # not installed: skipped where frmtmb cannot load without the library
  # that holds frmtmb.survival (R CMD check)
  gone <- child_with_libs(absent_package_libs("frmtmb.survival"), code)
  expect_match(msg(gone), paste0(
    "Install frmtmb.survival and refit the model with ",
    "library(frmtmb.survival) loaded"), fixed = TRUE)
})

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
