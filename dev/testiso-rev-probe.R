# Reviewer probe for lane testiso (2026-10-14). Not part of the lane.
#   Rscript dev/testiso-rev-probe.R
# 1. With all seven extensions loaded: unresolved sides, contributed
#    rules that win no pair, contributed method rows and where (if
#    anywhere) a function of that name lives.
# 2. local_frmtmb_registries() when the scope loads an extension.
# 3. local_frmtmb_registries() when a registry name it reads is absent.
.libPaths(c("C:/Users/adf44/source/r/wt-testiso-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
sys.source("C:/Users/adf44/source/r/frmtmb-wt-testiso/tests/testthat/helper-registry.R",
           envir = environment())

cat("== part 1: all seven loaded\n")
ext <- c("frmtmb.coupling", "frmtmb.eam", "frmtmb.latent", "frmtmb.learn",
         "frmtmb.ode", "frmtmb.sample", "frmtmb.spline")
cat("before loads: eam loaded =", isNamespaceLoaded("frmtmb.eam"), "\n")
for (p in ext) loadNamespace(p)
rules <- frm_compat_rules()
un <- ns$compat_unresolved_sides(rules, ns$frmtmb_compat_features_tbl()$name)
cat("unresolved sides with all seven:", if (length(un)) un else "none", "\n")
pairs <- ns$frmtmb_compat_pairs_tbl()
win <- ns$compat_resolve(pairs, rules)$win
dead <- setdiff(seq_len(nrow(rules)), unique(win))
cat("dead rules with all seven (no excuse):", length(dead), "\n")
for (k in dead) cat("  dead:", rules$feature_a[k], "x", rules$feature_b[k], "\n")
contrib <- names(unlist(ns$frmtmb_compat_contrib$features))
ft <- frm_compat_features()
cm <- ft$key[ft$kind == "method" & ft$name %in% contrib]
cat("contributed method rows:", length(cm), "\n")
for (m in cm) {
  where <- Filter(function(p) exists(m, envir = asNamespace(p), inherits = FALSE),
                  ext)
  cat("  method", m, "found in:", if (length(where)) where else "NOWHERE",
      "| on search path:", exists(m, mode = "function"), "\n")
}

cat("== part 2: a scope that loads an extension\n")
# fresh child: frmtmb loaded, frmtmb.eam not, then the scope loads it
body <- '
.libPaths(c("C:/Users/adf44/source/r/wt-testiso-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
sys.source("C:/Users/adf44/source/r/frmtmb-wt-testiso/tests/testthat/helper-registry.R",
           envir = globalenv())
f <- function() {
  local_frmtmb_registries()
  loadNamespace("frmtmb.eam")
  cat("inside: dec registered =",
      "dec" %in% names(frmtmb:::frmtmb_aterm_registry$reg), "\n")
}
f()
cat("after: eam loaded =", isNamespaceLoaded("frmtmb.eam"),
    "| dec registered =",
    "dec" %in% names(frmtmb:::frmtmb_aterm_registry$reg), "\n")
r <- tryCatch({
  frm(bf(y | dec(d) ~ x) + gaussian(), data = NULL,
      dry_run = "spec"); "parses"
}, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 120)))
cat("after: y | dec(d) with frmtmb.eam loaded:", r, "\n")
'
tf <- tempfile(fileext = ".R"); writeLines(body, tf)
out <- system2(file.path(R.home("bin"), "Rscript.exe"), tf, stdout = TRUE,
               stderr = TRUE)
cat(out, sep = "\n")

cat("== part 3: a registry name the helper reads is absent\n")
# what ns$<renamed> gives, and whether the helper notices
cat("ns$<absent> is NULL:", is.null(ns$frmtmb_no_such_registry), "\n")
g <- function() {
  aterms <- ns$frmtmb_no_such_registry
  saved <- aterms$reg
  ok <- tryCatch({ aterms$reg <- saved; "no error" },
                 error = function(e) conditionMessage(e))
  ok
}
cat("save and restore through an absent name:", g(), "\n")
