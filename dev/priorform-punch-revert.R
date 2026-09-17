# Fail-closed check for the punch round's guards: put the pre-punch
# behavior back inside a lane process, one piece at a time, and run the
# test files that should notice.
#   Rscript dev/priorform-punch-revert.R <mode> <filter>
# mode: none, order (later wins), slots (no duplicate refusal), twins
# (no twin collapse), slash (no brms group order), specials (cs() only);
# round 2: slotlabel, fillup, mmkey, smtwin, update; samplenone is the
# unmutated control for the frmtmb.sample test files
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
filt <- args[2]
.libPaths(c("C:/Users/adf44/source/r/priorform-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
put <- function(nm, f) {
  environment(f) <- ns
  utils::assignInNamespace(nm, f, ns = "frmtmb")
}
# a one-string edit of a function's source, refused unless the string
# occurs exactly once, so a mutant cannot silently change nothing
mut <- function(nm, from, to) {
  txt <- deparse(get(nm, ns))
  stopifnot(sum(grepl(from, txt, fixed = TRUE)) == 1L)
  put(nm, eval(parse(text = sub(from, to, txt, fixed = TRUE))))
}
on_sample <- mode %in% c("update", "samplenone")
pkg_dir <- if (on_sample) {
  file.path("extensions", "frmtmb.sample")
} else {
  "."
}
if (on_sample) suppressMessages(library(frmtmb.sample))
invisible(switch(mode,
  none = NULL,
  samplenone = NULL,
  order = put("prior_specificity_order", function(pl) unclass(pl)),
  slots = put("check_prior_slots", function(prior) invisible(prior)),
  twins = put("drop_twin_bar_terms", function(form) form),
  slash = put("slash_nested_bars", function(rest, sf, specials, env) {
    rep(FALSE, length(sf$reTrmFormulas))
  }),
  specials = {
    old <- get("brms_whole_term_specials", ns)
    unlockBinding("brms_whole_term_specials", ns)
    assign("brms_whole_term_specials", "cs", envir = ns)
  },
  # punch round 2
  slotlabel = mut("check_prior_slots", "vapply(specs, prior_spec_slot, \"\")",
                  "vapply(specs, prior_slot_label, \"\")"),
  fillup = mut("fill_prior_table", "seq_len(nrow(tab)) != i)",
               "seq_len(nrow(tab)) < i)"),
  mmkey = mut("refuse_duplicated_re", "if (!is.null(cp[[\"mm\"]]))",
              "if (FALSE)"),
  smtwin = mut("parse_linpred", "duplicated(sp_key)", "FALSE"),
  update = {
    sns <- asNamespace("frmtmb.sample")
    f <- get("drop_superseded", sns)
    txt <- deparse(f)
    old_txt <- "tg <- vapply(unclass(over), spec_target, \"\")"
    stopifnot(sum(grepl(old_txt, txt, fixed = TRUE)) == 1L)
    txt <- sub(old_txt, paste0("tg <- vapply(Filter(function(s) ",
                               "!is.null(s$dist), unclass(over)), ",
                               "spec_target, \"\")"), txt, fixed = TRUE)
    g <- eval(parse(text = txt))
    environment(g) <- sns
    utils::assignInNamespace("drop_superseded", g, ns = "frmtmb.sample")
  },
  stop("unknown mode ", mode)))
root <- "C:/Users/adf44/source/r/frmtmb-wt-priorform"
res <- testthat::test_dir(file.path(root, pkg_dir, "tests", "testthat"),
                          filter = filt,
                          package = if (on_sample)
                            "frmtmb.sample" else "frmtmb",
                          reporter = "silent", stop_on_failure = FALSE)
d <- as.data.frame(res)
cat(sprintf("PUNCH-REVERT %s %s pass %d fail %d error %d\n", mode, filt,
            sum(d$passed), sum(d$failed), sum(d$error)))
bad <- d$test[d$failed > 0 | d$error]
if (length(bad)) writeLines(paste("   ", unique(bad)))
