# brms_families.Rmd is a parameterization reference with no runnable
# chunks, so it is audited by name coverage: every family the document
# names, checked against what frm() will accept.
ROOT <- "C:/Users/adf44/source/r/frmtmb-wt-audit"
suppressMessages(pkgload::load_all(ROOT, quiet = TRUE))
src <- readLines(file.path(system.file("doc", package = "brms"),
                           "brms_families.Rmd"), warn = FALSE)
bold <- unlist(regmatches(src, gregexpr("\\*\\*[A-Za-z_0-9.]+\\*\\*", src)))
bold <- unique(gsub("\\*", "", bold))

# brms's own catalogue, to catch families the vignette names only in prose
brms_fams <- tryCatch(sort(unique(brms::brmsfamily("gaussian")$family)),
                      error = function(e) character())
catalogue <- tryCatch({
  f <- getFromNamespace("family_names", "brms")
  sort(unique(f()))
}, error = function(e) character())

accepted <- function(nm) {
  ok <- tryCatch({
    frmtmb:::as_frmtmb_family(nm)
    TRUE
  }, error = function(e) FALSE)
  ok
}

# Keep only tokens that look like family names (drop bolded prose)
cand <- bold[nchar(bold) > 2]
cand <- setdiff(cand, c("mu", "sigma", "not", "Note", "brms", "hurdle"))
res <- data.frame(
  family = cand,
  # multinomial is registered but needs its category count, so the
  # bare-name probe fails where frm(family = multinomial(K)) works
  frmtmb = vapply(cand, function(n) accepted(n) ||
                    n %in% names(frmtmb:::family_registry), logical(1)),
  row.names = NULL
)
res <- res[order(!res$frmtmb, res$family), ]
print(res, row.names = FALSE)
cat("\ncovered:", sum(res$frmtmb), "of", nrow(res), "\n")
cat("\nfrmtmb registry:\n")
cat(paste(sort(unique(names(frmtmb:::family_registry))), collapse = ", "), "\n")
if (length(catalogue)) {
  cat("\nbrms catalogue not in frmtmb registry:\n")
  cat(paste(setdiff(catalogue, names(frmtmb:::family_registry)), collapse = ", "), "\n")
}
