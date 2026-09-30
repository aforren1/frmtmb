# Reviewer runner: one test file, one arm, optionally one mutant.
#   Rscript dev/postfit2-rev-runtest.R <pkg> <file> <base|lane> [--gated]
#           [--ns] [--mutant=<name>]
# --ns runs the file in a child of the package NAMESPACE, as R CMD
# check's test_check() does; without it the environment is
# testthat::test_env(pkg), which is the lane's runner.
# Mutants are injected in the process (no install): each one names the
# expression it replaces and the run aborts if the replacement did not
# happen, so a mutant cannot silently be a no-op.
args <- commandArgs(TRUE)
pkg <- args[1]
file <- args[2]
arm <- args[3]
gated <- "--gated" %in% args
nsmode <- "--ns" %in% args
mutant <- sub("^--mutant=", "", grep("^--mutant=", args, value = TRUE))
if (gated) {
  Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true",
             FRMTMB_DRMTMB_FIT_TESTS = "true", FRMTMB_FUZZ = "true")
}
if (!nzchar(Sys.getenv("TMP"))) {
  Sys.setenv(TMP = "C:/Users/adf44/AppData/Local/Temp",
             TEMP = "C:/Users/adf44/AppData/Local/Temp")
}
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
sp <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad"
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(sp, "pf2rev-stan-cache"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2"
dir <- if (pkg == "frmtmb") file.path(wt, "tests/testthat") else
  file.path(wt, "extensions", pkg, "tests/testthat")
suppressMessages(library(testthat))
suppressMessages(library(pkg, character.only = TRUE))
cat("arm", arm, "mutant", if (length(mutant)) mutant else "none",
    "ns", nsmode, ":", pkg, format(packageVersion(pkg)), "from",
    find.package(pkg), "\n")

replace_expr <- function(e, from, to, count) {
  if (identical(e, from)) {
    count$n <- count$n + 1L
    return(to)
  }
  if (is.call(e)) {
    for (i in seq_along(e)) {
      if (!is.null(e[[i]]) && !identical(e[[i]], quote(expr = ))) {
        e[[i]] <- replace_expr(e[[i]], from, to, count)
      }
    }
  }
  e
}
drop_list_arg <- function(e, argname, count) {
  if (is.call(e)) {
    if (identical(e[[1]], as.name("list")) && argname %in% names(e)) {
      count$n <- count$n + 1L
      e <- e[names(e) != argname]
    }
    for (i in seq_along(e)) {
      if (!is.null(e[[i]]) && !identical(e[[i]], quote(expr = ))) {
        e[[i]] <- drop_list_arg(e[[i]], argname, count)
      }
    }
  }
  e
}
put_fn <- function(name, fn, ns_name) {
  ns <- asNamespace(ns_name)
  environment(fn) <- ns
  utils::assignInNamespace(name, fn, ns = ns_name)
  # an extension that imports frmtmb holds its own copy of an exported
  # binding, so the mutant has to reach it there too
  if (ns_name == "frmtmb" && isNamespaceLoaded("frmtmb.sample")) {
    imp <- parent.env(asNamespace("frmtmb.sample"))
    if (exists(name, envir = imp, inherits = FALSE)) {
      unlockBinding(name, imp)
      assign(name, fn, envir = imp)
      cat("  (also replaced in frmtmb.sample's imports)\n")
    }
  }
  fn
}
mutate_fn <- function(name, ns_name, from, to) {
  fn <- get(name, envir = asNamespace(ns_name))
  count <- new.env()
  count$n <- 0L
  body(fn) <- replace_expr(body(fn), from, to, count)
  if (count$n < 1L) stop("mutant ", mutant, ": no replacement in ", name)
  cat("mutant", mutant, ":", count$n, "replacement(s) in", ns_name, "::",
      name, "\n")
  put_fn(name, fn, ns_name)
}
reg_method <- function(generic, cls, fn) {
  for (owner in c("brms", "frmtmb")) {
    if (isNamespaceLoaded(owner) &&
        exists(generic, envir = asNamespace(owner))) {
      registerS3method(generic, cls, fn, envir = asNamespace(owner))
    }
  }
}
if (length(mutant)) {
  if (mutant == "setvars0") {
    fn <- function(conditions) character(0)
    put_fn("ce_set_vars", fn, "frmtmb")
    if (isNamespaceLoaded("frmtmb.sample")) {
      imp <- parent.env(asNamespace("frmtmb.sample"))
      unlockBinding("ce_set_vars", imp)
      environment(fn) <- asNamespace("frmtmb")
      assign("ce_set_vars", fn, envir = imp)
    }
    cat("mutant setvars0: ce_set_vars returns character(0)\n")
  } else if (mutant == "simNA") {
    mutate_fn("ce_boot_draws", "frmtmb",
              quote(if (!pop && !length(nspec)) NULL else NA), quote(NA))
  } else if (mutant == "keydrop") {
    fn <- get("ce_boot_key", envir = asNamespace("frmtmb"))
    count <- new.env()
    count$n <- 0L
    body(fn) <- drop_list_arg(body(fn), "conditional", count)
    if (count$n < 1L) stop("mutant keydrop: nothing dropped")
    cat("mutant keydrop:", count$n, "list argument(s) dropped\n")
    put_fn("ce_boot_key", fn, "frmtmb")
  } else if (mutant == "fitrevert") {
    fn <- mutate_fn("conditional_effects.frmtmb_fit", "frmtmb",
                    quote(setdiff(ce_group_vars(x), ce_set_vars(conditions))),
                    quote(ce_group_vars(x)))
    reg_method("conditional_effects", "frmtmb_fit", fn)
  } else if (mutant == "drawsrevert") {
    fn <- mutate_fn("conditional_effects.frmtmb_draws", "frmtmb.sample",
                    quote(setdiff(ce_group_vars(fit),
                                  ce_set_vars(conditions))),
                    quote(ce_group_vars(fit)))
    reg_method("conditional_effects", "frmtmb_draws", fn)
  } else if (mutant == "p1_unseen_obs") {
    # an unseen level counts as observed
    mutate_fn("ce_level_plan", "frmtmb",
              quote(!miss & key %in% b$bk[["levels"]]), quote(!miss))
  } else if (mutant == "p1_nonest") {
    # drop the nested rule: an observed block's columns may be moved
    mutate_fn("ce_plan_part", "frmtmb",
              quote(union(locked, unlist(lapply(blocks[-newk], `[[`, "gv")))),
              quote(locked))
  } else if (mutant == "p1_noshare") {
    # a new level drawn afresh in every part instead of shared
    mutate_fn("ce_plan_eval", "frmtmb", quote(is.null(cache[[ck]])),
              quote(TRUE))
  } else if (mutant == "p1_norefuse") {
    # the observed-and-new refusal of the bootstrap removed
    mutate_fn("ce_plan_kept", "frmtmb",
              quote(s[["observed"]] && s[["new"]]), quote(FALSE))
  } else if (mutant == "p1_nohold") {
    # the observed blocks are redrawn in the simulation after all
    mutate_fn("ce_boot_draws", "frmtmb",
              quote(setdiff(re_plan$redraw, held)), quote(re_plan$redraw))
  } else if (mutant == "p1_keptkey") {
    # the held set left out of the reuse key
    fn <- get("ce_boot_key", envir = asNamespace("frmtmb"))
    count <- new.env()
    count$n <- 0L
    body(fn) <- drop_list_arg(body(fn), "kept", count)
    if (count$n < 1L) stop("mutant p1_keptkey: nothing dropped")
    cat("mutant p1_keptkey:", count$n, "list argument(s) dropped\n")
    put_fn("ce_boot_key", fn, "frmtmb")
  } else if (mutant == "p1_allblocks") {
    # re_formula's partial keep ignored: every group block counts
    mutate_fn("ce_kept_blocks", "frmtmb", quote(intersect(grp, unique(kept))),
              quote(grp))
  } else if (mutant == "p1_redrawsmooth") {
    # the withdrawn user decision: NA redraws every block, smooths too
    mutate_fn("sim_re_plan", "frmtmb", quote(sim_group_block_ids(fit)),
              quote(seq_along(blocks)))
  } else if (mutant == "p1_nolock") {
    # nothing counts as read outside the grouping terms
    put_fn("ce_locked_vars", function(fit) character(0), "frmtmb")
    cat("mutant p1_nolock: ce_locked_vars returns character(0)\n")
  } else if (mutant == "p1_draw0") {
    # a new level drawn as zero (the population curve under a new name)
    mutate_fn("ce_plan_draw", "frmtmb", quote(mvn_draw_cov(as.matrix(S))),
              quote(numeric(nrow(as.matrix(S)))))
  } else if (mutant == "p2_byall") {
    # every row reads every by-level's block again
    put_fn("ce_by_reads", function(bk, nd) rep(TRUE, nrow(nd)), "frmtmb")
    cat("mutant p2_byall: ce_by_reads returns TRUE for every row\n")
  } else if (mutant == "p2_mmsplit") {
    # new members with the same value drawn independently, their
    # weights not added
    mutate_fn("ce_plan_part", "frmtmb", quote(unique(vals[isnew])),
              quote(vals[isnew]))
    mutate_fn("ce_plan_part", "frmtmb", quote(match(vals, u)),
              quote(cumsum(isnew)))
    # and a key per member, so the draw cache does not re-share them
    mutate_fn("ce_plan_part", "frmtmb", quote(paste0("new:", u[j])),
              quote(paste0("new:", u[j], "#", j)))
  } else if (mutant == "p2_rekey") {
    # the re_formula text left out of the reuse key
    put_fn("ce_re_key", function(re_form) "NA", "frmtmb")
    cat("mutant p2_rekey: ce_re_key returns a constant\n")
  } else {
    stop("unknown mutant ", mutant)
  }
}
env <- if (nsmode) new.env(parent = asNamespace(pkg)) else
  testthat::test_env(pkg)
r <- testthat::test_file(file.path(dir, file), package = pkg, env = env,
                         reporter = testthat::ListReporter$new(),
                         stop_on_failure = FALSE)
df <- as.data.frame(r)
cat(sprintf("RESULT %s %s %s%s: tests=%d failed=%d error=%d skipped=%d warning=%d passed=%d\n",
            arm, pkg, file, if (length(mutant)) paste0(" [", mutant, "]") else "",
            sum(df$nb), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning), sum(df$passed)))
for (res in r) {
  for (e in res$results) {
    if (inherits(e, c("expectation_failure", "expectation_error",
                      "expectation_warning"))) {
      sr <- e$srcref
      loc <- if (is.null(sr)) "?" else as.character(sr[1])
      cat("----", class(e)[1], "at line", loc, "in", res$test, "\n")
      cat(substr(conditionMessage(e), 1, 800), "\n")
    }
  }
}
sk <- unique(unlist(lapply(r, function(res) {
  vapply(Filter(function(e) inherits(e, "expectation_skip"), res$results),
         conditionMessage, "")
})))
if (length(sk)) cat("SKIP reasons:", paste(head(sk, 8), collapse = " | "), "\n")
