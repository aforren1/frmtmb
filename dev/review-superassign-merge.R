# Review of lane superassign: merge hazards with the other lanes of the
# round. Reads the -U0 diffs written to dev/review-superassign-merge/
# (one per worktree, base 3cd2a0c8), maps every hunk to the top-level
# definition of the base file it falls in, and reports, for every
# function superassign changes, the other lanes that change it and how
# close their hunks come (base line ranges). Also lists every `<<-` an
# other lane ADDS, and every added line of another lane that names a
# variable superassign moved into an environment inside the same
# function.
base_root <- "C:/Users/adf44/source/r/frmtmb"
mdir <- "C:/Users/adf44/source/r/frmtmb-wt-superassign/dev/review-superassign-merge"
lanes <- c("fitedge", "nlstart", "readers", "loglik", "nestquad", "arstep",
           "bigx", "osa", "latperf", "mfxse", "ci")

read_hunks <- function(lane) {
  x <- readLines(file.path(mdir, paste0(lane, ".diff")), warn = FALSE)
  x <- sub("\r$", "", x)
  out <- list()
  file <- NA_character_
  cur <- NULL
  for (ln in x) {
    if (startsWith(ln, "+++ ")) {
      file <- sub("^\\+\\+\\+ b/", "", ln)
    } else if (startsWith(ln, "@@")) {
      m <- regmatches(ln, regexec("^@@ -([0-9]+)(,([0-9]+))? \\+([0-9]+)(,([0-9]+))? @@", ln))[[1]]
      a <- as.integer(m[2]); b <- if (nzchar(m[4])) as.integer(m[4]) else 1L
      # a pure insertion after line a has b = 0
      lo <- if (b == 0L) a + 0.5 else a
      hi <- if (b == 0L) a + 0.5 else a + b - 1L
      out[[length(out) + 1L]] <- list(file = file, lo = lo, hi = hi,
                                      added = character(0), removed = character(0))
    } else if (length(out) && startsWith(ln, "+") && !startsWith(ln, "+++")) {
      out[[length(out)]]$added <- c(out[[length(out)]]$added, substring(ln, 2))
    } else if (length(out) && startsWith(ln, "-") && !startsWith(ln, "---")) {
      out[[length(out)]]$removed <- c(out[[length(out)]]$removed, substring(ln, 2))
    }
  }
  out
}

defs_cache <- new.env()
defs_of <- function(file) {
  if (!is.null(defs_cache[[file]])) return(defs_cache[[file]])
  p <- file.path(base_root, file)
  if (!grepl("[.][Rr]$", file)) p <- "/nonexistent"
  if (!file.exists(p)) {
    defs_cache[[file]] <- data.frame(name = character(0), lo = integer(0),
                                     hi = integer(0))
    return(defs_cache[[file]])
  }
  ex <- parse(p, keep.source = TRUE)
  sr <- attr(ex, "srcref")
  nm <- vapply(ex, function(e) {
    if (is.call(e) && as.character(e[[1L]])[1L] %in% c("<-", "=") &&
          is.name(e[[2L]])) as.character(e[[2L]]) else NA_character_
  }, "")
  d <- data.frame(name = nm,
                  lo = vapply(sr, function(s) s[1L], 1L),
                  hi = vapply(sr, function(s) s[3L], 1L))
  d <- d[!is.na(d$name), ]
  defs_cache[[file]] <- d
  d
}
fn_of <- function(file, lo, hi) {
  d <- defs_of(file)
  hit <- d$name[d$lo <= hi & d$hi >= lo]
  if (!length(hit)) "(outside any definition)" else hit
}

sa <- read_hunks("superassign")
sa_tab <- do.call(rbind, lapply(sa, function(h) {
  data.frame(file = h$file, fn = fn_of(h$file, h$lo, h$hi), lo = h$lo,
             hi = h$hi)
}))
# the old local names superassign moved, by function
moved <- list(
  autoscale_choose = c("w", "m"), fit_fd_se = "geta", ce_level_plan = "status",
  compat_rule_builder = "rows", compat_hand_rules_tbl = "rows",
  par_alias_index = c("nms", "pos"), varcorr_trans_rows = "rows",
  hyp_env_vals = "env", curve_catch_span = "msgs",
  nlminb_trial_fn = c("n_eval", "n_mapped", "last_par", "last_val",
                      "best_par", "best_val"),
  newton_polish = c("M", "Mi", "decided"), mo_search = "runs",
  resolve_group_calls = "fr", lookup_structural = "e2",
  frame_raw_vars = c("off_only", "other"),
  assemble_frame = c("rhs_comb", "rhs_resp", "codes_new"),
  frm_bind_generic = c("memo_ns", "memo_fn"), gp_brms_natural = "out",
  frm_hazard_reads = "out", build_importance_objective = "offs",
  me_newdata_extra = "obs", me_node_value = "memo",
  VarCorr.frmtmb_fit = "pos", nl_restart_round = "h1",
  nl_fit_off_saddle = "held", mo_bar_rewrite = "map",
  expand_double_verts = "written", parse_linpred = c("dblbar_written", "protected"),
  rhs_split_ac = "ac", lp_basis_nl = "parts", as_priorlist = "bad",
  priorlist_from_rows = "bad", prior_table = "rows", cov_class_coefs = "rows",
  resolve_priorlist = c("assigned", "lower", "upper"), ps_extract = "terms",
  re_formula_bars = "out", nat_slots = c("slots", "dup", "legacy"),
  frm_ode = "gave_up", log_lik.frmtmb_draws = "warned",
  draws_laplace_watch = "first", draws_unseen_level = "what",
  prior_plan = "out", default_priors_for = "pl",
  tmbstan_build_broken = "cached", cox_bhaz_args = "out")

cat("lib: none (diff analysis)\n")
cat("superassign hunks:", nrow(sa_tab), " functions:",
    length(unique(paste(sa_tab$file, sa_tab$fn))), "\n\n")
cat("== `<<-` ADDED by other lanes (R/ and extensions/)\n")
n_sup <- 0L
for (ln in lanes) {
  for (h in read_hunks(ln)) {
    a <- h$added[grepl("<<-|->>", h$added)]
    for (s in a) {
      n_sup <- n_sup + 1L
      cat(sprintf("  %-9s %s [%s]: %s\n", ln, h$file,
                  paste(fn_of(h$file, h$lo, h$hi), collapse = ","), trimws(s)))
    }
  }
}
cat("  total:", n_sup, "\n\n")

cat("== per superassign function: other lanes' hunks in the same definition\n")
key <- unique(sa_tab[, c("file", "fn")])
rows <- list()
for (i in seq_len(nrow(key))) {
  f <- key$file[i]; fn <- key$fn[i]
  mine <- sa_tab[sa_tab$file == f & sa_tab$fn == fn, ]
  for (ln in lanes) {
    hs <- Filter(function(h) identical(h$file, f) &&
                   fn %in% fn_of(h$file, h$lo, h$hi), read_hunks(ln))
    if (!length(hs)) next
    # distance in base lines between any hunk of each side; 0 means
    # overlap, <= 1 adjacent (git merge-file reports a conflict)
    dist <- min(vapply(hs, function(h) {
      min(vapply(seq_len(nrow(mine)), function(k) {
        max(0, max(mine$lo[k], h$lo) - min(mine$hi[k], h$hi))
      }, 1))
    }, 1))
    added <- unlist(lapply(hs, `[[`, "added"))
    old <- moved[[fn]] %||% character(0)
    uses <- character(0)
    for (v in old) {
      pat <- paste0("(^|[^A-Za-z0-9_.$])", gsub(".", "\\.", v, fixed = TRUE),
                    "($|[^A-Za-z0-9_.])")
      if (any(grepl(pat, added))) uses <- c(uses, v)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      file = f, fn = fn, lane = ln, hunks = length(hs), min_dist = dist,
      textual = if (dist <= 1) "CONFLICT" else if (dist <= 3) "near" else "apart",
      old_names_in_added = paste(uses, collapse = ","))
  }
}
tab <- do.call(rbind, rows)
`%||%` <- function(a, b) if (is.null(a)) b else a
print(tab, row.names = FALSE, right = FALSE)
cat("\nRESULT functions shared with another lane:",
    length(unique(paste(tab$file, tab$fn))), "; textual conflicts:",
    sum(tab$textual == "CONFLICT"), "; added lines naming a moved variable:",
    sum(nzchar(tab$old_names_in_added)), "\n")
