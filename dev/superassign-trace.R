# Record what every function this lane changed RETURNS while one test
# file runs, so the base build and the lane build can be compared call
# by call.
#
# Each changed function (and run_optimizer(), which reads the state
# nlminb_trial_fn() keeps) is traced on exit. The return value is put in
# a canonical form and hashed: a function or an environment becomes a
# placeholder, because the lane changed what closures hold, not what
# they compute, and an environment's address differs between processes;
# a formula or a call becomes its deparsed text; a double is rounded to
# 9 significant digits, because two R processes on this hybrid CPU can
# differ by a few ulps (dev/rtmb-pitfalls.md item 21). The sequence of
# (function, hash) pairs is written as a TSV.
#
# Usage: Rscript dev/superassign-trace.R <arm: base|lane> <pkg> <testfile>
#        <out.tsv>

a <- commandArgs(trailingOnly = TRUE)
arm <- a[1L]
pkg <- a[2L]
tf <- a[3L]
out_file <- a[4L]
lane <- "C:/Users/adf44/source/r/wt-superassign-lib"
base <- "C:/Users/adf44/source/r/rellib-r18"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (arm == "lane") lane, base, user))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER"))) Sys.setenv(R_MAKEVARS_USER = mk)

suppressMessages(library(testthat))
suppressMessages(library(pkg, character.only = TRUE))
cat("lib:", dirname(find.package("frmtmb")), "\n")

targets <- list(
  frmtmb = c("as_priorlist", "assemble_frame", "autoscale_choose",
             "ce_level_plan", "compat_hand_rules_tbl", "cov_class_coefs",
             "curve_catch_span", "expand_double_verts", "fit_fd_se",
             "frame_raw_vars", "frm_hazard_reads", "gp_brms_natural",
             "hyp_env_vals", "lookup_structural", "lp_basis_nl",
             "me_newdata_extra", "mo_bar_rewrite", "mo_search", "nat_slots",
             "newton_polish", "nl_fit_off_saddle", "nl_restart_round",
             "par_alias_index", "parse_linpred", "prior_table",
             "priorlist_from_rows", "ps_extract", "re_formula_bars",
             "resolve_group_calls", "resolve_priorlist", "rhs_split_ac",
             "VarCorr.frmtmb_fit", "varcorr_trans_rows", "run_optimizer",
             "frm_compat"),
  frmtmb.sample = c("prior_plan", "default_priors_for",
                    "draws_unseen_level"),
  frmtmb.survival = c("cox_bhaz_args")
)

canon <- function(x, depth = 0L) {
  if (depth > 40L) return("<deep>")
  if (is.null(x)) return(NULL)
  if (is.function(x)) return("<fn>")
  if (is.environment(x)) return("<env>")
  if (typeof(x) %in% c("externalptr", "weakref", "bytecode")) {
    return("<ptr>")
  }
  if (isS4(x)) {
    m <- tryCatch(as.matrix(x), error = function(e) NULL)
    return(if (is.null(m)) "<S4>" else canon(m, depth + 1L))
  }
  if (is.language(x)) return(paste(deparse(x), collapse = " "))
  at <- attributes(x)
  if (is.list(x)) {
    y <- lapply(x, canon, depth = depth + 1L)
  } else if (is.double(x)) {
    y <- signif(as.vector(x), 9L)
  } else if (is.complex(x)) {
    y <- "<complex>"
    at <- NULL
  } else {
    y <- as.vector(x)
  }
  if (length(at)) {
    at <- lapply(at, canon, depth = depth + 1L)
    keep <- intersect(names(at), c("names", "dim", "dimnames", "class",
                                   "row.names", "levels"))
    rest <- setdiff(names(at), keep)
    attributes(y) <- at[keep]
    if (length(rest)) y <- list(value = y, attrs = at[rest])
  }
  y
}

rec <- new.env(parent = emptyenv())
rec$rows <- vector("list", 4096L)
rec$n <- 0L
record <- function(name, val) {
  h <- tryCatch(digest::digest(canon(val)),
                error = function(e) paste0("<err:", conditionMessage(e),
                                           ">"))
  rec$n <- rec$n + 1L
  if (rec$n > length(rec$rows)) length(rec$rows) <- 2L * length(rec$rows)
  rec$rows[[rec$n]] <- c(name, h)
}

for (p in names(targets)) {
  if (!isNamespaceLoaded(p)) next
  ns <- asNamespace(p)
  s3 <- ns[[".__S3MethodsTable__."]]
  for (nm in targets[[p]]) {
    if (!exists(nm, envir = ns, inherits = FALSE)) {
      cat("MISSING", p, nm, "\n")
      next
    }
    ex <- substitute(record(N, returnValue()), list(N = nm))
    suppressMessages(trace(nm, exit = ex, where = ns, print = FALSE))
    # S3 dispatch reads the methods table, not the namespace binding
    if (exists(nm, envir = s3, inherits = FALSE)) {
      assign(nm, get(nm, envir = ns), envir = s3)
    }
  }
}

res <- testthat::test_file(tf, package = pkg,
                           env = testthat::test_env(pkg),
                           reporter = "silent")
df <- as.data.frame(res)
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d\n", basename(tf),
            sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped)))
rows <- rec$rows[seq_len(rec$n)]
tab <- if (length(rows)) {
  data.frame(i = seq_along(rows), fn = vapply(rows, `[`, "", 1L),
             hash = vapply(rows, `[`, "", 2L))
} else {
  data.frame(i = integer(0), fn = character(0), hash = character(0))
}
utils::write.table(tab, out_file, sep = "\t", quote = FALSE,
                   row.names = FALSE)
cat("CALLS", nrow(tab), "\n")
