# lane tmbstan, probe 02: which test_that blocks actually reach the
# sampler, and which skip they take.
#
# Static. Parses the suite with srcrefs, collects the call names in
# each block and in each top-level function, then takes a fixpoint so
# a block that calls a cached case closure counts as reaching the
# sampler. Sampler mentions inside expect_error() do not count: the
# compatibility pre-flight fires BEFORE check_tmbstan_build() in both
# doors, so a refusal test still passes on a broken build.
dir <- paste0("C:/Users/adf44/source/r/frmtmb-wt-tmbstan/",
              "extensions/frmtmb.sample/tests/testthat")
SAMPLERS <- c("frm_sample", "as_tmbstan", "tmbstan")
SKIPPED_PARENTS <- c("expect_error", "expect_snapshot")

names_in <- function(e, skip_error, acc = character()) {
  if (!is.call(e)) return(acc)
  f <- e[[1L]]
  if (is.name(f)) {
    nm <- as.character(f)
    if (skip_error && nm %in% SKIPPED_PARENTS) return(acc)
    acc <- c(acc, nm)
  } else if (is.call(f) && identical(as.character(f[[1L]]), "::")) {
    acc <- c(acc, as.character(f[[3L]]))
  }
  for (i in seq_along(e)) {
    if (!is.null(e[[i]]) && (is.call(e[[i]]) || is.expression(e[[i]]))) {
      acc <- names_in(e[[i]], skip_error, acc)
    }
  }
  acc
}

files <- sort(list.files(dir, pattern = "^(test|helper)-.*[.]R$"))
parsed <- lapply(files, function(f) parse(file.path(dir, f),
                                          keep.source = TRUE))
names(parsed) <- files

# top-level assignments, per file, so a shadowing local definition is
# visible as such
def_names <- list()
for (f in files) {
  for (e in parsed[[f]]) {
    if (is.call(e) && as.character(e[[1L]])[1L] %in% c("<-", "=") &&
        is.name(e[[2L]])) {
      key <- paste(f, as.character(e[[2L]]))
      def_names[[key]] <- names_in(e[[3L]], TRUE)
    }
  }
}

# transitive reach over the same graph, with a file-local definition
# shadowing the helper's, which is exactly what three files do
reaches <- function(seed_names, file, targets) {
  seen <- character(); todo <- seed_names; via <- NA_character_
  while (length(todo)) {
    n <- todo[1L]; todo <- todo[-1L]
    if (n %in% seen) next
    seen <- c(seen, n)
    if (n %in% targets) return(n)
    for (src in c(file, "helper-sampling.R", "helper-brms.R")) {
      k <- paste(src, n)
      if (!is.null(def_names[[k]])) {
        todo <- c(todo, def_names[[k]])
        break
      }
    }
  }
  NA_character_
}

# which definition of skip_sampler a file resolves to
skip_owner <- function(file) {
  if (!is.null(def_names[[paste(file, "skip_sampler")]])) return("LOCAL")
  if (!is.null(def_names[["helper-sampling.R skip_sampler"]])) {
    return("helper")
  }
  "none"
}

rows <- list()
for (f in files) {
  srefs <- attr(parsed[[f]], "srcref")
  for (i in seq_along(parsed[[f]])) {
    e <- parsed[[f]][[i]]
    if (!(is.call(e) && identical(as.character(e[[1L]]), "test_that"))) {
      next
    }
    body <- e[[3L]]
    ln <- if (is.null(srefs)) NA_integer_ else srefs[[i]][1L]
    nm <- as.character(e[[2L]])
    all_nm <- names_in(body, FALSE)
    live_nm <- names_in(body, TRUE)
    rows[[length(rows) + 1L]] <- data.frame(
      file = f, line = ln, test = substr(nm, 1L, 52L),
      samples = !is.na(reaches(live_nm, f, SAMPLERS)),
      mentions = !is.na(reaches(all_nm, f, SAMPLERS)),
      via_skip = reaches(all_nm, f, "skip_sampler"),
      inline = reaches(all_nm, f, "skip_if_not_installed"),
      owner = skip_owner(f),
      stringsAsFactors = FALSE)
  }
}
d <- do.call(rbind, rows)
d$route <- ifelse(!is.na(d$via_skip),
                  ifelse(d$owner == "LOCAL", "skip_sampler(LOCAL)",
                         "skip_sampler(helper)"),
                  ifelse(!is.na(d$inline), "inline skip_if_not_installed",
                         "NO SKIP"))

samp <- d[d$samples, ]
cat("== sampling blocks by skip route ==
")
print(table(samp$route))
cat("
== sampling blocks NOT covered by the shared helper ==
")
un <- samp[samp$route != "skip_sampler(helper)", ]
print(table(un$file, un$route))
cat("
n uncovered =", nrow(un), " n covered =",
    sum(samp$route == "skip_sampler(helper)"), "
")
cat("
== sampling blocks with NO SKIP AT ALL ==
")
ns <- samp[samp$route == "NO SKIP", ]
print(ns[, c("file", "line", "test")], row.names = FALSE)
cat("
== non-sampling blocks that mention a sampler (refusal tests) ==
")
amb <- d[!d$samples & d$mentions, ]
print(amb[, c("file", "line", "test", "route")], row.names = FALSE)
cat("
== totals ==
")
cat("test_that blocks:", nrow(d), "  sampling blocks:", sum(d$samples),
    "
")
cat("files with a LOCAL skip_sampler:
")
for (f in files) {
  if (!is.null(def_names[[paste(f, "skip_sampler")]]) &&
      f != "helper-sampling.R") cat("  ", f, "
")
}
