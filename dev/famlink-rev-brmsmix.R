## Reviewer probe: brms's own mixture() component rules, read from brms
## 2.23.0, and brms's default link_<dpar> values against frmtmb's
## constructor defaults (a brms family object passed to frm() now
## carries them).
.libPaths(c("C:/Users/adf44/source/r/famlink-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
f <- deparse(brms::mixture)
i <- grep("not allowed|real and integer|ordinal", f)
cat(f[sort(unique(unlist(lapply(i, function(j) max(1, j - 8):(j + 1)))))], sep = "\n")
ns <- asNamespace("brms")
cat("\n---- brms families whose info lacks 'mix' in specials ----\n")
for (n in ls(ns, all.names = TRUE, pattern = "^[.]family_")) {
  x <- tryCatch(get(n, ns)(), error = function(e) NULL)
  if (is.null(x)) next
  if (!"mix" %in% x$specials) cat(sub("^[.]family_", "", n), "")
}
cat("\n\n---- brms default link_<dpar> vs frmtmb constructor default ----\n")
suppressMessages(library(frmtmb))
src <- frmtmb:::brms_mu_link_source
n_diff <- 0L; n_cmp <- 0L
for (nm in names(src)) {
  if (nm %in% frmtmb:::brms_mu_link_analogs || nm == "multinomial") next
  bf <- tryCatch(brms::brmsfamily(src[[nm]]), error = function(e) NULL)
  ff <- tryCatch(frmtmb:::family_registry[[nm]](), error = function(e) NULL)
  if (is.null(bf) || is.null(ff)) { cat(nm, ": not built\n"); next }
  bl <- grep("^link_", names(bf), value = TRUE)
  for (l in bl) {
    n_cmp <- n_cmp + 1L
    fv <- unclass(ff)[[l]]
    if (!identical(bf[[l]], fv)) {
      n_diff <- n_diff + 1L
      cat(sprintf("%-26s %-12s brms %-10s frmtmb %s\n", nm, l, bf[[l]],
                  if (is.null(fv)) "<absent>" else fv))
    }
  }
  # does frm() accept the brms object as is?
  r <- tryCatch({frmtmb:::as_frmtmb_family(bf); "ok"},
                error = function(e) conditionMessage(e))
  if (!identical(r, "ok")) cat(sprintf("%-26s as_frmtmb_family(brms object) ERROR %s\n", nm, substr(r, 1, 120)))
}
cat("link_<dpar> compared:", n_cmp, " differing:", n_diff, "\n")

cat("\n---- brms:::no_mixture per brms family ----\n")
print(brms:::no_mixture)
for (n in ls(ns, all.names = TRUE, pattern = "^[.]family_")) {
  b <- sub("^[.]family_", "", n)
  fam <- tryCatch(brms::brmsfamily(b), error = function(e) NULL)
  if (is.null(fam)) next
  if (isTRUE(brms:::no_mixture(fam))) cat("no_mixture:", b, "\n")
}
