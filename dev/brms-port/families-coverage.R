# brms_families.Rmd is a parameterization reference with no runnable
# chunks, so it is audited by name coverage: every family the document
# names, checked against what frm() will accept.
#
#   PORT_LIB=<lib> Rscript families-coverage.R
HERE <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
source(file.path(HERE, "env.R"))
suppressMessages(library(frmtmb))
cat("#", port_build(), "\n")
src <- readLines(file.path(system.file("doc", package = "brms"),
                           "brms_families.Rmd"), warn = FALSE)
bold <- unlist(regmatches(src, gregexpr("\\*\\*[A-Za-z_0-9.]+\\*\\*", src)))
bold <- unique(gsub("\\*", "", bold))

catalogue <- tryCatch({
  f <- getFromNamespace("family_names", "brms")
  sort(unique(f()))
}, error = function(e) character())

# The probe the 0.34.0 audit used, kept so the counts compare.
accepted <- function(nm) {
  tryCatch({
    frmtmb:::as_frmtmb_family(nm)
    TRUE
  }, error = function(e) FALSE)
}

# Keep only tokens that look like family names (drop bolded prose)
cand <- bold[nchar(bold) > 2]
cand <- setdiff(cand, c("mu", "sigma", "not", "Note", "brms", "hurdle"))
reg <- names(frmtmb:::family_registry)
res <- data.frame(
  family = cand,
  # multinomial is registered but needs its category count, so the
  # bare-name probe fails where frm(family = multinomial(K)) works
  frmtmb = vapply(cand, function(n) accepted(n) || n %in% reg, logical(1)),
  row.names = NULL
)
res <- res[order(!res$frmtmb, res$family), ]
print(res, row.names = FALSE)
cat("\ncovered:", sum(res$frmtmb), "of", nrow(res), "\n")
cat("missing:", paste(res$family[!res$frmtmb], collapse = ", "), "\n")
cat("\nfrmtmb registry:\n")
cat(paste(sort(unique(reg)), collapse = ", "), "\n")
if (length(catalogue)) {
  cat("\nbrms catalogue (", length(catalogue), " names) not in the frmtmb",
      " registry:\n", sep = "")
  cat(paste(setdiff(catalogue, reg), collapse = ", "), "\n")
}
