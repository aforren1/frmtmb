# Generate the brms suite tier (item 2.6b) from brms 2.23.0's own source.
#
#   Rscript dev/brmsport-gen.R
#
# Reads the bin-1 blocks (dev/brmsport-blocks.R), the verdicts
# (dev/brmsport-verdicts.tsv) and writes
#   tests/testthat/test-brms-suite-<topic>.R                  core, both
#   extensions/frmtmb.sample/tests/testthat/test-brms-suite-<topic>.R
#                                                            sample, both
# Every assertion is brms's own source text, wrapped in brms_port() with
# its verdict. An assertion with no verdict row is generated as
# "unclassified", which fails a normal run and is allowed only while
# recording (dev/brmsport-ledger.R). A setup line is wrapped in
# brms_setup(). The ONLY rewrites are the substitutions in `subst`
# below; an assertion they touch is marked adapted in the ledger.
source("dev/brmsport-blocks.R")

root <- "."
core_dir <- file.path(root, "tests/testthat")
sample_dir <- file.path(root, "extensions/frmtmb.sample/tests/testthat")

topic_of <- function(f) sub("^tests[.](.*)[.]R$", "\\1", f)

# Per brms file: the output topic, the lines bound at the top of the
# generated file, and the rewrites.
files <- list(
  "tests.families.R" = list(topic = "families", pre = character()),
  "tests.brmsformula.R" = list(topic = "brmsformula", pre = character()),
  "tests.brmsterms.R" = list(topic = "brmsterms", pre = character()),
  "tests.brm.R" = list(
    topic = "brm",
    pre = c("brm <- brms_shim_brm",
            "inhaler <- brms::inhaler")),
  "tests.standata.R" = list(
    topic = "standata",
    pre = c("standata <- brms_shim_standata",
            "SW <- suppressWarnings")),
  "tests.priors.R" = list(
    topic = "priors",
    pre = c("epilepsy <- brms::epilepsy",
            "inhaler <- brms::inhaler")),
  "tests.brmsfit-helpers.R" = list(
    topic = "brmsfit-helpers",
    pre = c("epilepsy <- brms::epilepsy")),
  "tests.brmsfit-methods.R" = list(
    topic = "methods",
    pre = c(
      "expect_range <- function(object, lower = -Inf, upper = Inf, ...) {",
      "  testthat::expect_true(all(object >= lower & object <= upper), ...)",
      "}",
      "expect_ggplot <- function(object, ...) {",
      "  testthat::expect_true(is(object, \"ggplot\"), ...)",
      "}",
      "SM <- suppressMessages",
      "SW <- suppressWarnings",
      "epilepsy <- brms::epilepsy",
      "rename_pars <- brms_shim_rename_pars",
      "update <- brms_shim_update",
      "fit1 <- brms_fixture(1)",
      "fit2 <- brms_fixture(2)",
      "fit3 <- brms_fixture(3)",
      "fit4 <- brms_fixture(4)",
      "fit5 <- brms_fixture(5)",
      "fit6 <- brms_fixture(6)",
      "nobs <- 40",
      "npatients <- 10",
      "nsubjects <- 8",
      "nvisits <- 4"),
    sample_pre = c(
      "fit1 <- brms_fixture_draws(1)",
      "fit2 <- brms_fixture_draws(2)",
      "fit3 <- brms_fixture_draws(3)",
      "fit5 <- brms_fixture_draws(5)")),
  "tests.emmeans.R" = list(
    topic = "emmeans",
    pre = c("skip_if_not_installed(\"emmeans\")",
            "library(emmeans)",
            "SW <- suppressWarnings",
            "rename_pars <- brms_shim_rename_pars",
            "fit1 <- brms_fixture(1)",
            "fit2 <- brms_fixture(2)",
            "fit4 <- brms_fixture(4)",
            "fit6 <- brms_fixture(6)")),
  "tests.data-helpers.R" = list(
    topic = "data-helpers",
    pre = character(),
    subst = list(
      c("brms:::rename_pars(brms:::brmsfit_example1)", "brms_fixture(1)"),
      c("brms:::validate_newdata(", "brms_shim_validate_newdata(")))
)

verdicts <- if (file.exists("dev/brmsport-verdicts.tsv")) {
  utils::read.delim("dev/brmsport-verdicts.tsv", stringsAsFactors = FALSE,
                    quote = "")
} else {
  data.frame(id = character(), pkg = character(), verdict = character(),
             reason = character())
}
stopifnot(!anyDuplicated(verdicts[c("id", "pkg")]))

own <- utils::read.delim("dev/brmsport-verdicts-own.tsv", quote = "",
                         colClasses = "character", na.strings = NULL)
stopifnot(!anyDuplicated(own[c("id", "pkg")]))
own_for <- function(id, pkg) {
  o <- own[own$id == id & own$pkg %in% c(pkg, "*"), ]
  if (nrow(o)) o else NULL
}

verdict_for <- function(id, pkg) {
  v <- verdicts[verdicts$id == id & verdicts$pkg %in% c(pkg, "*"), ]
  if (nrow(v) == 0L) return(list(verdict = "unclassified", reason = ""))
  if (nrow(v) > 1L) v <- v[v$pkg == pkg, ]
  stopifnot(nrow(v) == 1L)
  list(verdict = v$verdict, reason = v$reason)
}

q <- function(x) encodeString(x, quote = "\"")

# A reason is emitted as short string pieces so no line passes 80
# columns.
emit_reason <- function(reason, indent) {
  if (!nzchar(reason)) return(paste0(indent, "\"\","))
  pieces <- strwrap(reason, width = 60)
  pieces <- paste0(pieces, c(rep(" ", length(pieces) - 1L), ""))
  if (length(pieces) == 1L) return(paste0(indent, q(pieces), ","))
  c(paste0(indent, "paste0("),
    paste0(indent, "  ", q(pieces),
           c(rep(",", length(pieces) - 1L), "),")))
}

# brms's source text, verbatim when it fits in 80 columns at this
# indent, and deparsed otherwise.
code_lines <- function(text, expr, indent, col = 1L) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  if (length(lines) > 1L) {
    # keep brms's alignment relative to where the statement starts
    lead <- strrep(" ", col - 1L)
    rest <- lines[-1]
    has <- startsWith(rest, lead)
    rest[has] <- substring(rest[has], col)
    lines[-1] <- rest
  }
  lines <- paste0(indent, lines)
  if (all(nchar(lines) <= 80L)) return(lines)
  d <- deparse(expr, width.cutoff = 60L)
  paste0(indent, d)
}

apply_subst <- function(text, subst) {
  for (s in subst) text <- gsub(s[1], s[2], text, fixed = TRUE)
  text
}

blocks <- brms_bin1_blocks()
out <- list()
ids_seen <- character()
ids_by_pkg <- character()
for (b in blocks) {
  cfg <- files[[b$file]]
  stopifnot(!is.null(cfg))
  pkgs <- switch(b$tier, core = "frmtmb", sample = "frmtmb.sample",
                 both = c("frmtmb", "frmtmb.sample"))
  for (pkg in pkgs) {
    key <- paste(pkg, cfg$topic)
    body <- c(sprintf("test_that(%s, {", q(b$label)))
    for (j in seq_along(b$stmts)) {
      s <- b$stmts[[j]]
      id <- sprintf("%s:%d", topic_of(b$file), b$lines[j])
      text <- apply_subst(b$texts[j], cfg$subst)
      adapted <- !identical(text, b$texts[j])
      expr <- if (adapted) str2lang(text) else s
      k <- n_expect(s)
      if (k == 0L) {
        if (is.call(s) && as.character(s[[1]])[1] %in%
            c("skip_on_cran", "skip_if_not_installed", "context")) {
          body <- c(body, code_lines(text, expr, "  ", b$cols[j]))
        } else if (grepl("brms:::?", text)) {
          # running it would put brms's answer where frmtmb's belongs
          nms <- if (is.call(expr) &&
                     as.character(expr[[1]]) %in% c("<-", "=")) {
            all.vars(expr[[2]])
          } else {
            character()
          }
          body <- c(body,
                    sprintf("  brms_setup_not_run(%s, %s)", q(id),
                            q(paste(nms, collapse = " "))),
                    paste0("  # ", strsplit(text, "\n")[[1]]))
        } else if (is.call(expr) && identical(expr[[1]], as.name("="))) {
          # a top-level `=` would bind as a named argument of brms_setup()
          body <- c(body, sprintf("  brms_setup(%s, {", q(id)),
                    code_lines(text, expr, "    ", b$cols[j]), "  })")
        } else {
          body <- c(body, sprintf("  brms_setup(%s,", q(id)),
                    code_lines(text, expr, "    ", b$cols[j]), "  )")
        }
        next
      }
      ids_seen <- c(ids_seen, id)
      ids_by_pkg <- c(ids_by_pkg, paste(pkg, id))
      v <- verdict_for(id, pkg)
      if (grepl("brms:::?", text)) {
        # reaches brms's own namespace: running it would test brms
        body <- c(body,
                  sprintf("  brms_port_not_run(%s, %s,", q(id),
                          q(v$verdict)),
                  emit_reason(v$reason, "    "),
                  "    \"reaches brms's own namespace\"",
                  "  )",
                  paste0("  # ", strsplit(text, "\n")[[1]]))
        next
      }
      o <- own_for(id, pkg)
      if (!is.null(o)) {
        stopifnot(!grepl("  ", o$pattern))
        body <- c(body, sprintf("  brms_port_own(%s,", q(id)),
                  emit_reason(o$pattern, "    "),
                  emit_reason(o$note, "    "),
                  code_lines(text, expr, "    ", b$cols[j]), "  )")
        next
      }
      body <- c(body, sprintf("  brms_port(%s, %s,", q(id), q(v$verdict)),
                emit_reason(v$reason, "    "),
                code_lines(text, expr, "    ", b$cols[j]), "  )")
    }
    body <- c(body, "})", "")
    out[[key]] <- c(out[[key]], body)
  }
}
stopifnot(!anyDuplicated(ids_by_pkg))

# brms's methods file calls plot() and pairs(), which draw on the
# current device; without one, R opens Rplots.pdf in the tests
# directory, and R CMD build ships whatever is there. Measured, not
# grepped: the sample package's methods file left one on 0.62.0
# (dev/simnewdata-log/rplots-base.txt).
device_guard <- c(
  "grDevices::pdf(NULL)",
  "withr::defer(grDevices::dev.off(), teardown_env())"
)

header <- function(f, pkg, cfg) {
  pre <- cfg$pre
  if (pkg == "frmtmb.sample" && !is.null(cfg$sample_pre)) {
    pre <- c(pre[!grepl("^fit[0-9] <-", pre)], cfg$sample_pre)
  }
  c("# GENERATED by dev/brmsport-gen.R from brms 2.23.0",
    sprintf("# tests/testthat/%s. Do not edit by hand: change", f),
    "# dev/brmsport-verdicts.tsv or the generator, then regenerate.",
    "# helper-brms-suite.R says what brms_port() asserts.",
    "",
    "skip_unless_brms_suite()",
    if (pkg == "frmtmb.sample" && cfg$topic == "methods") "skip_sampler()",
    if (cfg$topic == "methods") device_guard,
    pre,
    "")
}

written <- character()
for (key in names(out)) {
  pkg <- sub(" .*$", "", key)
  topic <- sub("^[^ ]+ ", "", key)
  f <- names(files)[vapply(files, function(x) x$topic, "") == topic]
  dir <- if (pkg == "frmtmb") core_dir else sample_dir
  path <- file.path(dir, sprintf("test-brms-suite-%s.R", topic))
  txt <- c(header(f, pkg, files[[f]]), out[[key]])
  while (length(txt) && !nzchar(txt[length(txt)])) txt <- txt[-length(txt)]
  writeLines(txt, path, useBytes = TRUE)
  written <- c(written, path)
}

# the sample package's tests need the same harness; one source of truth
h <- readLines(file.path(core_dir, "helper-brms-suite.R"))
writeLines(c("# COPIED by dev/brmsport-gen.R from core's",
             "# tests/testthat/helper-brms-suite.R. Edit that one.", "", h),
           file.path(sample_dir, "helper-brms-suite.R"))

cat("wrote", length(written), "files;", length(unique(ids_seen)),
    "assertion ids;",
    sum(!unique(ids_seen) %in% verdicts$id), "without a verdict\n")
cat(written, sep = "\n")
