# Each new translation rule and the REML marginal identity, broken on
# purpose: a rule that cannot fail would score any map as an identity.
# Every mutant must turn its pair's L1 from pass to fail. Records nothing.
#
#   Rscript dev/translator-mutants.R
#
# From the worktree root; reads FRMTMB_EQUIV_LIB as dev/equiv-measure.R.
lib <- Sys.getenv("FRMTMB_EQUIV_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                               winslash = "/"))
}
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
  requireNamespace("brms")
  requireNamespace("rstan")
})
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
fresh <- function() {
  e <- new.env(parent = asNamespace("frmtmb"))
  for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
    sys.source(f, envir = e)
  }
  e
}
smp <- read.delim("dev/equiv-data/sample.tsv")
run <- function(label, id, mutate = NULL, spec_edit = NULL) {
  e <- fresh()
  if (is.function(mutate)) mutate(e)
  row <- smp[smp$id == id, ]
  key <- paste(row$feature_a, row$feature_b, sep = " x ")
  if (is.function(spec_edit)) {
    mk <- e$eq_specs[[key]]
    e$eq_specs[[key]] <- function() spec_edit(mk())
  }
  r <- e$eq_measure(key, row)
  cat(sprintf(paste0("%-40s %s  L1 ok %s  points %s  verdict %s  const %s",
                     "  marg %s  error %s\n"), label, id,
              format(r$l1_ok %||% NA), format(r$l1_points %||% NA),
              r$verdict %||% "",
              paste(signif(c(r$l1_const, r$l1_extra_const), 3),
                    collapse = ","),
              paste(signif(r$l1_marg %||% NA, 3), collapse = ","),
              substr(paste(r$l1_error %||% "",
                           paste(r$l1_extra_err, collapse = " | ")), 1, 100)))
  invisible(r)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# the components' category-specific effects swapped
run("bcs: components swapped", "p131", function(e) {
  v <- e$eq_rules$bcs_prefixed$value
  e$eq_rules$bcs_prefixed$value <- function(nm, ...) {
    other <- if (grepl("1$", nm)) sub("1$", "2", nm) else sub("2$", "1", nm)
    v(other, ...)
  }
})
# the components' thresholds swapped
run("thresholds: components swapped", "p131", function(e) {
  v <- e$eq_rules$ord_mix_intercept$value
  e$eq_rules$ord_mix_intercept$value <- function(nm, ...) {
    other <- if (grepl("1$", nm)) sub("1$", "2", nm) else sub("2$", "1", nm)
    v(other, ...)
  }
})
# tmp_xi moved by 0.1 on its own scale
run("tmp_xi: off by 0.1", "p021", function(e) {
  v <- e$eq_rules$tmp_xi$value
  e$eq_rules$tmp_xi$value <- function(...) {
    out <- v(...)
    out$value <- out$value + 0.1
    out
  }
})
# a smooth's coefficients in reverse order, in an nl body and under mvbf
for (id in c("p042", "p048")) {
  run("zs: coefficients reversed", id, function(e) {
    v <- e$eq_rules$zs_prefixed$value
    e$eq_rules$zs_prefixed$value <- function(...) {
      out <- v(...)
      out$value <- array(rev(out$value), length(out$value))
      out
    }
  })
}
# one t2() block at a time (review 2026-11-09, M1): at the survey's points
# a block whose SD collapsed to about 1e-4 could be reversed unseen; at
# the specs' fixed points every block must be seen
one_block <- function(block) {
  function(e) {
    v <- e$eq_rules$zs_prefixed$value
    e$eq_rules$zs_prefixed$value <- function(nm, ...) {
      out <- v(nm, ...)
      if (identical(nm, block)) {
        out$value <- array(rev(out$value), length(out$value))
      }
      out
    }
  }
}
for (b in sprintf("zs_y1_1_%d", 1:3)) {
  run(paste("zs: block", b, "alone reversed"), "p048", one_block(b))
}
for (b in sprintf("zs_b1_1_%d", 1:3)) {
  run(paste("zs: block", b, "alone reversed"), "p047", one_block(b))
}
# one row of a group-level z at a time (short check of 2026-11-10, F1):
# p087's mmc() slope has SD 1.1e-6 at the estimate, so its row could be
# reversed unseen at the survey's points; at the fixed points it must
# fail. Reversing a row permutes N(0, 1) entries, so the prior is exact.
one_row <- function(block, r) {
  function(e) {
    prev <- e$stan_pars_from_fit
    e$stan_pars_from_fit <- function(...) {
      p <- prev(...)
      z <- p[[block]]
      z[r, ] <- rev(z[r, ])
      p[[block]] <- z
      p
    }
  }
}
for (r in 1:2) {
  run(sprintf("z_1 row %d alone reversed", r), "p087", one_row("z_1", r))
}
# the other pairs the per-row audit found (a correlation at +/-1 leaves
# row 2, and p070's row 3, with almost no weight at the estimate), each
# with its blind row reversed at the fixed points, then unmutated
deg <- c("p053", "p056", "p060", "p069", "p070", "p102", "p106", "p111",
         "p114", "p115", "p116", "p122")
for (id in deg) {
  run("z_1 row 2 alone reversed", id, one_row("z_1", 2))
}
run("z_1 row 3 alone reversed", "p070", one_row("z_1", 3))
for (id in deg) run("control", id)
# a fixed point that cannot be read (review 2026-11-09, m4): the third
# point moved 3 instead of 0.1, where frm() or the tmp_xi rule stops. The
# pair must fail, not pass on the two points left
run("fixed points: third point unreadable", "p021", function(e) {
  orig <- e$eq_jitter_optimizer
  e$eq_jitter_optimizer <- function(scale, k) {
    orig(if (k == 2) 3 else scale, k)
  }
})
# the REML marginal identity read on the ML fit of the same model: the
# integral over mu's coefficients is then not frm()'s objective
run("marginal: ML fit, not REML", "p044", spec_edit = function(s) {
  s$frm_args$REML <- NULL
  s
})
# and the controls, unmutated
for (id in c("p131", "p021", "p042", "p048", "p047", "p087", "p044")) {
  run("control", id)
}
