# Lane names: the brms vs frmtmb table, generated from the saved runs.
# Reads dev/names-out/brms-*.rds (dev/names-brms.R) and
# dev/names-out/frm-{base,lane}[-sample].rds (dev/names-frm.R), and for
# each case and each arm says whether the names both sides carry come in
# the same order, and whether ranef()/coef() have brms's entries and
# columns. Output: dev/names-out/table.md, pasted into the findings.
rd <- function(p) if (file.exists(p)) readRDS(p) else list()
br <- c(rd("dev/names-out/brms-cat-skew-dist-distpop-twoblk-smooth.rds"),
        rd("dev/names-out/brms-mo-gpby-me-latres-mv-dblbar.rds"))
dbl <- rd("dev/names-out/brms-dblbar.rds")
if (length(dbl)) br[["dblbar"]] <- dbl[["dblbar"]]
arms <- list(base = rd("dev/names-out/frm-base.rds"),
             lane = rd("dev/names-out/frm-lane.rds"),
             `base draws` = rd("dev/names-out/frm-base-sample.rds"),
             `lane draws` = rd("dev/names-out/frm-lane-sample.rds"))
same_order <- function(b, f) {
  common <- intersect(b, f)
  if (!length(common)) return("no common names")
  ok <- identical(b[b %in% common], f[f %in% common])
  sprintf("%s (%d common)", if (ok) "same" else "DIFFERS", length(common))
}
layout <- function(b, f) {
  if (is.character(b) && length(b) == 1L && is.null(names(b))) {
    # brms stopped
    return(if (is.character(f) && length(f) == 1L && is.null(names(f)))
      "both stop" else "brms stops, frmtmb answers")
  }
  if (is.character(f) && length(f) == 1L && is.null(names(f))) {
    return("frmtmb stops, brms answers")
  }
  bk <- lapply(b, function(d) d[[length(d)]])
  if (identical(names(bk), names(f)) &&
      identical(unname(bk), unname(lapply(f, as.character)))) "same" else
    paste0("DIFFERS: brms ", paste(names(bk), collapse = ","), " [",
           paste(vapply(bk, paste, "", collapse = ","), collapse = "; "),
           "] frmtmb ", paste(names(f), collapse = ","), " [",
           paste(vapply(f, function(x) paste(x, collapse = ","), ""),
                 collapse = "; "), "]")
}
out <- c("| case | arm | variables() order | ranef() | coef() |",
         "|---|---|---|---|---|")
for (cs in names(br)) {
  for (a in names(arms)) {
    f <- arms[[a]][[cs]]
    if (is.null(f)) next
    out <- c(out, sprintf("| %s | %s | %s | %s | %s |", cs, a,
                          same_order(br[[cs]]$variables, f$variables),
                          layout(br[[cs]]$ranef, f$ranef),
                          if (grepl("draws", a) ||
                              !is.list(f$coef) || is.null(names(f$coef)))
                            layout(br[[cs]]$coef, f$coef) else
                              layout(br[[cs]]$coef, f$coef)))
  }
}
writeLines(out, "dev/names-out/table.md")
cat(out, sep = "\n")
