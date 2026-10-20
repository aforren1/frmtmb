# What every contributed compatibility row promises, checked in a
# session that holds every extension. test-compat-extensions.R runs it
# in a child process; it is a helper so that the same function runs
# there, here and under a mutation check (dev/testiso-contract-mutants.R).
#
# test-compat.R checks core's own rows in whatever session it runs in.
# Three promises of the contributed rows can only be checked with every
# owner loaded:
# - no rule side is unresolved: an `expects =` name is exempt at
#   registration by design, so a misspelt one would wait for ever;
# - no rule wins zero pairs;
# - a contributed `method` row names a function a loaded frmtmb
#   package exports, or an argument setting of one, spelled as core
#   spells its own (`confint_profile`): `<function>_<argument>`, where
#   the argument is a formal argument of the function (frmtmb.spline's
#   `frm_curve_contrast`). Registration checks only the kind.
compat_contract_violations <- function() {
  ns <- asNamespace("frmtmb")
  rules <- frmtmb::frm_compat_rules()
  ft <- frmtmb::frm_compat_features()
  out <- character()
  un <- attr(frmtmb::frm_compat(), "unresolved")
  if (length(un)) out <- c(out, paste0("unresolved rule side: ", un))
  win <- ns$compat_resolve(ns$frmtmb_compat_pairs_tbl(), rules)$win
  dead <- setdiff(seq_len(nrow(rules)), unique(win))
  if (length(dead)) {
    out <- c(out, sprintf("rule wins no pair: %s x %s",
                          rules$feature_a[dead], rules$feature_b[dead]))
  }
  contrib <- names(unlist(ns$frmtmb_compat_contrib$features))
  ext <- grep("^frmtmb[.]", loadedNamespaces(), value = TRUE)
  exported_fn <- function(fn) {
    for (p in ext) {
      if (fn %in% getNamespaceExports(p)) {
        f <- getExportedValue(p, fn)
        if (is.function(f)) return(f)
      }
    }
    NULL
  }
  # `<function>_<argument>`: the row splits at an underscore into a
  # function and one of its formal arguments. Every split is tried,
  # because a function name holds underscores too
  argument_row <- function(m) {
    us <- gregexpr("_", m, fixed = TRUE)[[1]]
    for (i in us[us > 0L]) {
      f <- exported_fn(substr(m, 1L, i - 1L))
      a <- substr(m, i + 1L, nchar(m))
      if (is.function(f) && a != "..." && a %in% names(formals(f))) {
        return(TRUE)
      }
    }
    FALSE
  }
  for (m in ft$name[ft$kind == "method" & ft$name %in% contrib]) {
    ok <- is.function(exported_fn(m)) || argument_row(m)
    if (!ok) {
      out <- c(out, paste0("method row names no exported function or ",
                           "argument: ", m))
    }
  }
  out
}
