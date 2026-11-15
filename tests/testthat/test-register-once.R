# A frame check or a prior-default provider registered again replaces
# its own entry. `pkgload::load_all()` of an extension runs its
# .onLoad() a second time with new closures, and the two list registries
# used to append them, so each of its frame checks ran twice per fit for
# the rest of the session (backlog, testiso).
#
# "Again" has to be provable, because a wrong merge drops a check with no
# word while a wrong duplicate only runs one twice: the same code, and
# the same values in every environment between the function and the
# first named one (a namespace, the global environment), with no
# environment or external pointer among those values. Anything else is
# kept. The factories below live in a namespace, as a package's own
# would; one made inside this test would reach the test's environment,
# which holds environments, and so would never merge.

in_ns <- function(f) {
  environment(f) <- asNamespace("stats")
  f
}

test_that("a frame check registered again runs once, in its place", {
  local_frmtmb_registries()
  ns <- asNamespace("frmtmb")
  n0 <- length(ns$frmtmb_frame_checks$fns)
  # the same code and the same (no) captured values in two closures is
  # what a reload hands over
  mk <- in_ns(function() function(spec, frame) invisible(NULL))
  first <- mk()
  frmtmb_register_frame_check(first)
  other <- in_ns(function(spec, frame, ...) invisible(NULL))
  frmtmb_register_frame_check(other)
  expect_length(ns$frmtmb_frame_checks$fns, n0 + 2L)
  second <- mk()
  frmtmb_register_frame_check(second)
  expect_length(ns$frmtmb_frame_checks$fns, n0 + 2L)
  # the new closure took the old one's place, ahead of `other`: its
  # environment is the one held there, and the first closure's is gone
  # (on a build that appends, the first one is still in that slot)
  envs <- lapply(ns$frmtmb_frame_checks$fns[n0 + 1:2], environment)
  expect_identical(envs[[1L]], environment(second))
  expect_false(any(vapply(envs, identical, NA, environment(first))))
  # a namespace function is known by its namespace and its code
  frmtmb_register_frame_check(stats::model.frame)
  frmtmb_register_frame_check(stats::model.frame)
  expect_length(ns$frmtmb_frame_checks$fns, n0 + 3L)
})

test_that("two checks one factory made for different values stay two", {
  # the absent case the review constructed: a key of code alone dropped
  # the refusing check and the fit ran
  local_frmtmb_registries()
  ns <- asNamespace("frmtmb")
  n0 <- length(ns$frmtmb_frame_checks$fns)
  make_check <- in_ns(function(tag) {
    force(tag)
    function(spec, frame) {
      if (tag == "refuse") stop("refused by check ", tag, call. = FALSE)
      invisible(NULL)
    }
  })
  frmtmb_register_frame_check(make_check("refuse"))
  frmtmb_register_frame_check(make_check("log_only"))
  expect_length(ns$frmtmb_frame_checks$fns, n0 + 2L)
  set.seed(1)
  d <- data.frame(x = stats::rnorm(40), y = stats::rnorm(40))
  expect_error(frm(y ~ x, data = d), "refused by check refuse")
  # the same factory value again is a reload of the first, not a third
  frmtmb_register_frame_check(make_check("refuse"))
  expect_length(ns$frmtmb_frame_checks$fns, n0 + 2L)
})

test_that("a value one frame further up keeps two checks apart", {
  # the short check's cases: the value that differs is not in the
  # closure's own environment, so a key of that environment alone merged
  # them and the refusing check was lost
  local_frmtmb_registries()
  ns <- asNamespace("frmtmb")
  n0 <- length(ns$frmtmb_frame_checks$fns)
  outer <- in_ns(function(mode) {
    force(mode)
    function(label) {
      force(label)
      function(spec, frame) {
        if (mode == "refuse") stop("refused, nested", call. = FALSE)
      }
    }
  })
  frmtmb_register_frame_check(outer("refuse")("l"))
  frmtmb_register_frame_check(outer("log")("l"))
  expect_length(ns$frmtmb_frame_checks$fns, n0 + 2L)
  # the same values all the way up are one check
  frmtmb_register_frame_check(outer("log")("l"))
  expect_length(ns$frmtmb_frame_checks$fns, n0 + 2L)
  # local() inside one factory: the closure's own frame is empty
  via_local <- in_ns(function(mode) {
    force(mode)
    local(function(spec, frame) {
      if (mode == "refuse") stop("refused, local", call. = FALSE)
    })
  })
  frmtmb_register_frame_check(via_local("refuse"))
  frmtmb_register_frame_check(via_local("log"))
  expect_length(ns$frmtmb_frame_checks$fns, n0 + 4L)
  set.seed(1)
  d <- data.frame(x = stats::rnorm(40), y = stats::rnorm(40))
  expect_error(frm(y ~ x, data = d), "refused, nested")
})

test_that("an environment or a pointer among the values is never merged", {
  # serialize() writes an environment by its contents and a pointer
  # without its address, so equal bytes do not prove the same object
  local_frmtmb_registries()
  ns <- asNamespace("frmtmb")
  n0 <- length(ns$frmtmb_frame_checks$fns)
  mk_cfg <- in_ns(function(cfg) {
    force(cfg)
    function(spec, frame) if (isTRUE(cfg$refuse)) stop("refused, cfg")
  })
  e1 <- new.env()
  frmtmb_register_frame_check(mk_cfg(e1))
  frmtmb_register_frame_check(mk_cfg(new.env()))
  expect_length(ns$frmtmb_frame_checks$fns, n0 + 2L)
  mk_ptr <- in_ns(function(ptr) {
    force(ptr)
    function(spec, frame) invisible(ptr)
  })
  p1 <- get("C_Cdqrls", asNamespace("stats"))$address
  p2 <- get("C_DoubleCentre", asNamespace("stats"))$address
  frmtmb_register_frame_check(mk_ptr(p1))
  frmtmb_register_frame_check(mk_ptr(p2))
  expect_length(ns$frmtmb_frame_checks$fns, n0 + 4L)
  # the first environment turned to refuse after it was registered: its
  # check is still there to say so
  e1$refuse <- TRUE
  set.seed(1)
  d <- data.frame(x = stats::rnorm(40), y = stats::rnorm(40))
  expect_error(frm(y ~ x, data = d), "refused, cfg")
})

test_that("a named environment that is not a namespace is never merged", {
  # attr(e, "name") or attach() gives an environment a name; stopping
  # the walk there keyed the closure by its code alone, and two checks
  # whose values differed merged (short check, fix round 2)
  local_frmtmb_registries()
  ns <- asNamespace("frmtmb")
  n0 <- length(ns$frmtmb_frame_checks$fns)
  named <- function(mode, nm) {
    e <- new.env(parent = globalenv())
    attr(e, "name") <- nm
    e$mode <- mode
    f <- function(spec, frame) {
      if (mode == "refuse") stop("refused, named", call. = FALSE)
    }
    environment(f) <- e
    f
  }
  frmtmb_register_frame_check(named("refuse", "cfg"))
  frmtmb_register_frame_check(named("log", "cfg"))
  expect_length(ns$frmtmb_frame_checks$fns, n0 + 2L)
  frmtmb_register_frame_check(named("refuse", "cfgA"))
  frmtmb_register_frame_check(named("log", "cfgB"))
  expect_length(ns$frmtmb_frame_checks$fns, n0 + 4L)
  named_parent <- function(mode) {
    e <- named(mode, "cfg")
    f <- function(spec, frame) {
      if (mode == "refuse") stop("refused, named parent", call. = FALSE)
    }
    environment(f) <- new.env(parent = environment(e))
    f
  }
  frmtmb_register_frame_check(named_parent("refuse"))
  frmtmb_register_frame_check(named_parent("log"))
  expect_length(ns$frmtmb_frame_checks$fns, n0 + 6L)
  set.seed(1)
  d <- data.frame(x = stats::rnorm(40), y = stats::rnorm(40))
  expect_error(frm(y ~ x, data = d), "refused, named")
})

test_that("checks with different code stay two checks", {
  # the guard's absent case: a key too coarse would drop a check
  local_frmtmb_registries()
  ns <- asNamespace("frmtmb")
  n0 <- length(ns$frmtmb_frame_checks$fns)
  frmtmb_register_frame_check(in_ns(function(spec, frame) 1))
  frmtmb_register_frame_check(in_ns(function(spec, frame) 2))
  expect_length(ns$frmtmb_frame_checks$fns, n0 + 2L)
})

test_that("a prior-default provider registered again is one provider", {
  local_frmtmb_registries()
  ns <- asNamespace("frmtmb")
  n0 <- length(ns$frmtmb_prior_defaults$providers)
  mk <- in_ns(function() function(spec, frame) NULL)
  frmtmb:::frmtmb_register_prior_defaults(mk())
  frmtmb:::frmtmb_register_prior_defaults(mk())
  expect_length(ns$frmtmb_prior_defaults$providers, n0 + 1L)
  frmtmb:::frmtmb_register_prior_defaults(
    in_ns(function(spec, frame) list()))
  expect_length(ns$frmtmb_prior_defaults$providers, n0 + 2L)
  # one factory, two values: two providers
  make_prov <- in_ns(function(lb) {
    force(lb)
    function(spec, frame) lb
  })
  frmtmb:::frmtmb_register_prior_defaults(make_prov(0.5))
  frmtmb:::frmtmb_register_prior_defaults(make_prov(0.25))
  expect_length(ns$frmtmb_prior_defaults$providers, n0 + 4L)
})
