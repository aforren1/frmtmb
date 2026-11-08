# Release review 0.80.0: is every lane's NEWS bullet (its final text) in
# the merged NEWS, word for word after whitespace is collapsed? Lists the
# bullets that are not, by lane and package, with the release's nearest
# bullet, for a reader to compare.
#
#   Rscript dev/review-rel-news.R
REL <- "C:/Users/adf44/source/r/frmtmb-wt-rel080"
lanes <- c("ci", "latperf", "mfxse", "superassign", "loglik", "nestquad",
           "readers", "fitedge", "osa", "bigx", "arstep", "nlstart")
files <- c("NEWS.md", file.path("extensions", c("frmtmb.sample",
  "frmtmb.latent", "frmtmb.ode", "frmtmb.survival", "frmtmb.learn",
  "frmtmb.spatial", "frmtmb.coupling", "frmtmb.eam"), "NEWS.md"))
top <- function(p) {
  if (!file.exists(p)) return(character(0))
  x <- readLines(p, warn = FALSE)
  h <- grep("^# ", x)
  if (length(h) < 2L) return(character(0))
  x[(h[1] + 1L):(h[2] - 1L)]
}
bullets <- function(x) {
  # a bullet runs to the next bullet, heading or blank line before a bullet
  s <- grep("^(\\* |  - |- |## )", x)
  out <- character(0)
  for (i in seq_along(s)) {
    e <- if (i < length(s)) s[i + 1L] - 1L else length(x)
    b <- paste(trimws(x[s[i]:e]), collapse = " ")
    b <- gsub("\\s+", " ", trimws(sub("^(\\*|-) ", "", b)))
    if (!grepl("^## ", x[s[i]]) && nzchar(b)) out <- c(out, b)
  }
  out
}
norm <- function(s) gsub("\\s+", " ", s)
for (f in files) {
  rel <- norm(paste(top(file.path(REL, f)), collapse = " "))
  rel_b <- bullets(top(file.path(REL, f)))
  base_has <- system2("git", c("-C", REL, "show", paste0("3cd2a0c8:", f)),
                      stdout = TRUE, stderr = FALSE)
  for (ln in lanes) {
    lp <- file.path(paste0("C:/Users/adf44/source/r/frmtmb-wt-", ln), f)
    x <- top(lp)
    if (!length(x)) next
    # a lane's top section is new only if it is not the base's top section
    if (length(base_has) && identical(x, {
      h <- grep("^# ", base_has); base_has[(h[1] + 1L):(h[2] - 1L)]
    })) next
    for (b in bullets(x)) {
      if (grepl(b, rel, fixed = TRUE)) next
      # the nearest release bullet by shared words
      w <- unique(strsplit(tolower(gsub("[^A-Za-z0-9.]+", " ", b)), " ")[[1]])
      sc <- vapply(rel_b, function(r) {
        wr <- unique(strsplit(tolower(gsub("[^A-Za-z0-9.]+", " ", r)), " ")[[1]])
        length(intersect(w, wr)) / max(1, length(w))
      }, 0)
      cat("\n== ", f, " lane ", ln, " (best match ", format(max(c(0, sc)),
          digits = 2), ")\nLANE: ", b, "\nREL:  ",
          if (length(sc)) rel_b[which.max(sc)] else "-", "\n", sep = "")
    }
  }
}
cat("\nNEWS CHECK DONE\n")
