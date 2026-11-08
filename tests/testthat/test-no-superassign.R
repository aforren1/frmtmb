# No `<<-` in the package code.
#
# rOpenSci's pkgcheck fails a package on "Package uses global assignment
# operator ('<<-')". Its pkgchk_left_assign() parses every file of R/
# and counts the LEFT_ASSIGN tokens spelled `<<-`. So any `<<-` written
# as code in R/ fails it, inside a function or at top level; one in a
# comment or a string is not a token and does not count. A closure that
# must keep state between calls writes to an environment it was given
# or created (`acc$rows[[i]] <- x`), and a walk that collects results
# returns them.
#
# Two scanners, for the reason test-bracket-access.R has two:
#
# * superassign_src_scan() reads the SOURCES with pkgcheck's own rule,
#   so it reports file and line and sees code outside any function. It
#   needs the monorepo checkout, so it skips under `R CMD check` on a
#   built tarball. It also counts `->>`, the same operator written the
#   other way round, which pkgcheck does not.
# * superassign_ns_scan() reads the installed NAMESPACE, so it runs
#   wherever the package does, including under `R CMD check` and covr.
#   It walks every function of the namespace and every function held
#   in a list or an environment there (registries such as autocor_core
#   and covstruct_core, S4 method tables, a closure's enclosure). A
#   `->>` parses to a `<<-` call, so it sees both spellings.

superassign_src_scan <- function(rdir) {
  out <- lapply(sort(list.files(rdir, pattern = "\\.[RrSsQq]$",
                                full.names = TRUE)), function(f) {
    pd <- utils::getParseData(parse(f, keep.source = TRUE))
    hit <- (pd$token == "LEFT_ASSIGN" & pd$text == "<<-") |
      (pd$token == "RIGHT_ASSIGN" & pd$text == "->>")
    sprintf("  %s:%d: %s", basename(f), pd$line1[hit], pd$text[hit])
  })
  as.character(unlist(out, use.names = FALSE))
}

superassign_ns_scan <- function(ns) {
  # An empty argument (`x[, 1]`) errors when TOUCHED, and an inner
  # function's formals are a pairlist, not a call; both as in
  # frm_hazard_reads()
  part <- function(e, i) tryCatch({
    a <- e[[i]]
    if (is.pairlist(a) && length(a)) as.call(c(quote(list), as.list(a)))
    else if (is.call(a)) a
  }, error = function(err) NULL)
  walk <- function(e, where) {
    here <- if (identical(e[[1L]], as.name("<<-"))) where
    kids <- lapply(seq_along(e), function(i) {
      a <- part(e, i)
      if (!is.null(a)) walk(a, where)
    })
    c(here, unlist(kids, use.names = FALSE))
  }
  own <- topenv(ns)
  # A function literal held in a registry (a list such as autocor_core,
  # an environment, an S4 method table, a closure's own enclosure) is
  # code as much as a top-level function, so lists and environments are
  # walked too. Environments can refer to each other in a cycle, so each
  # is walked once.
  seen <- new.env(parent = emptyenv())
  first_visit <- function(env) {
    key <- format(env)
    if (!is.null(seen[[key]])) return(FALSE)
    seen[[key]] <- TRUE
    TRUE
  }
  first_visit(ns)
  get_quiet <- function(nm, env) {
    tryCatch(get(nm, envir = env, inherits = FALSE),
             error = function(err) NULL)
  }
  walk_env <- function(env, prefix) {
    unlist(lapply(sort(ls(env, all.names = TRUE)), function(nm) {
      visit(get_quiet(nm, env), paste0(prefix, nm))
    }), use.names = FALSE)
  }
  visit <- function(x, where) {
    if (is.function(x)) {
      # an active binding resolves to ANOTHER package's generic when its
      # owner is loaded, so a function is judged by where it lives
      if (is.primitive(x) || !identical(topenv(environment(x)), own)) {
        return(NULL)
      }
      calls <- c(list(body(x)), Filter(function(d) {
        tryCatch(is.call(d), error = function(err) FALSE)
      }, as.list(formals(x))))
      hits <- unlist(lapply(calls, function(b) {
        if (is.call(b)) walk(b, where)
      }), use.names = FALSE)
      return(c(hits, visit(environment(x), paste0(where, "<env>$"))))
    }
    if (is.environment(x)) {
      if (isNamespace(x) || identical(x, globalenv()) ||
            identical(x, baseenv()) || identical(x, emptyenv()) ||
            !first_visit(x)) {
        return(NULL)
      }
      prefix <- if (endsWith(where, "$")) where else paste0(where, "$")
      return(walk_env(x, prefix))
    }
    if (is.list(x)) {
      nms <- names(x) %||% character(length(x))
      return(unlist(lapply(seq_along(x), function(i) {
        visit(x[[i]], paste0(where, if (nzchar(nms[i])) {
          paste0("$", nms[i])
        } else {
          paste0("[[", i, "]]")
        }))
      }), use.names = FALSE))
    }
    NULL
  }
  as.character(unique(walk_env(ns, "")))
}

superassign_advice <- paste0(
  "\n\nKeep state in an environment the closure is given or creates ",
  "(`acc <- new.env(parent = emptyenv())`, then `acc$x <- value`), or ",
  "return the value. Not `assign(..., envir = parent.frame())`, which ",
  "keeps the same hidden coupling. See the header of this file."
)

# the same positive identification of the source tree as
# test-bracket-access.R: fail closed (skip) rather than open
superassign_source_dir <- function() {
  rdir <- testthat::test_path("..", "..", "R")
  desc <- testthat::test_path("..", "..", "DESCRIPTION")
  ok <- file.exists(desc) &&
    any(trimws(readLines(desc, n = 5L)) == "Package: frmtmb") &&
    dir.exists(rdir) &&
    file.exists(file.path(rdir, "objective.R"))
  if (ok) rdir else NULL
}

test_that("no frmtmb function uses `<<-`", {
  found <- superassign_ns_scan(asNamespace("frmtmb"))
  expect_identical(
    found, character(0),
    info = paste0("These frmtmb functions use `<<-`:\n  ",
                  paste(found, collapse = "\n  "), superassign_advice)
  )
})

test_that("no core file under R/ uses `<<-`", {
  rdir <- superassign_source_dir()
  skip_if_not(!is.null(rdir),
              "package sources are not available (installed-package run)")
  found <- superassign_src_scan(rdir)
  expect_identical(
    found, character(0),
    info = paste0("A core file uses `<<-`:\n",
                  paste(found, collapse = "\n"), superassign_advice)
  )
})

test_that("no extension file under R/ uses `<<-`", {
  # pkgcheck runs on core alone today; the extensions are held to the
  # same rule so that none of them is the first to fail when it runs on
  # one. Skips where extensions/ is absent, as test-bracket-access.R's
  # backstop does.
  skip_if_not(!is.null(superassign_source_dir()),
              "package sources are not available (installed-package run)")
  extroot <- testthat::test_path("..", "..", "extensions")
  dirs <- if (dir.exists(extroot)) {
    list.dirs(extroot, recursive = FALSE)
  } else {
    character(0)
  }
  dirs <- dirs[dir.exists(file.path(dirs, "R"))]
  skip_if_not(length(dirs) > 0L, "no extensions in this tree")
  for (d in dirs) {
    found <- superassign_src_scan(file.path(d, "R"))
    expect_identical(
      found, character(0),
      info = paste0("An ", basename(d), " file uses `<<-`:\n",
                    paste(found, collapse = "\n"), superassign_advice)
    )
  }
})

test_that("both `<<-` scanners see a hit where one exists, and only there", {
  # The assertions above pass on an empty tree, so each scanner is shown
  # code with hits, beside the spellings that are not hits: a `<<-` in a
  # string and in a comment, and `<-` into an environment
  dir <- withr::local_tempdir()
  writeLines(c(
    "f <- function() {",
    "  n <- 0",
    "  acc <- new.env(parent = emptyenv())",
    "  s <- \"n <<- 1\"               # a string, not code",
    "  # n <<- 1                     # a comment, not code",
    "  acc$n <- 1                    # the replacement",
    "  g <- function(a = (n <<- 2)) a",
    "  h <- function() 3 ->> n",
    "  n <<- 4",
    "}"), file.path(dir, "control.R"))
  found <- superassign_src_scan(dir)
  expect_identical(found, c("  control.R:7: <<-", "  control.R:8: ->>",
                            "  control.R:9: <<-"))

  ctl <- new.env()
  eval(parse(file.path(dir, "control.R"), keep.source = FALSE), ctl)
  ctl$clean <- function() {
    acc <- new.env(parent = emptyenv())
    acc$n <- "n <<- 1"
    acc$n
  }
  expect_identical(superassign_ns_scan(ctl), "f")
  # an inner function's default argument is found when it is the only
  # hit, which a walk that descends into calls alone would miss
  ctl$f <- function() {
    g <- function(a = (n <<- 2)) a
    g
  }
  expect_identical(superassign_ns_scan(ctl), "f")
  ctl$f <- function() NULL
  expect_identical(superassign_ns_scan(ctl), character(0))
  # a function held in a registry list, nested as autocor_core's are,
  # and one held in an environment that refers to itself
  ctl$reg <- list(fcor = list(parse = function() {
    y <- 0
    g <- function() y <<- 1
    g()
  }), k = 2L)
  box <- new.env()
  box$self <- box
  box$fn <- function() {
    z <- 0
    h <- function() z <<- 1
    h()
  }
  ctl$box <- box
  expect_identical(sort(superassign_ns_scan(ctl)),
                   c("box$fn", "reg$fcor$parse"))
  ctl$reg$fcor$parse <- function() NULL
  rm("fn", envir = box)
  expect_identical(superassign_ns_scan(ctl), character(0))
})
