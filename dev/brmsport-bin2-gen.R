# Generate the bin-2 port of brms 2.23.0's suite from brms's own source.
#
#   Rscript dev/brmsport-bin2-gen.R
#
# The machinery is bin 1's (dev/brmsport-gen.R): every assertion is
# brms's source text wrapped in brms_port() with its verdict, a setup
# line in brms_setup(), a line reaching brms's namespace is recorded
# and never run. Two things differ.
#
# - The files are written under dev/brmsport-bin2/<package>/, not into
#   the packages' test directories: bin 2 is a measurement, and
#   dev/brmsport-bin2-run.R stages each file beside its package's test
#   helpers to run it. Nothing in tests/ changes.
# - Bin 2T asks about a posterior. Its blocks are tier "both", so each
#   runs against the maximum-likelihood fixture fit (frmtmb) and against
#   frm_sample() draws of it (frmtmb.sample), and the ledger decides each
#   row by the object the assertion is about (dev/brmsport-bin2-targets.tsv).
#   The draws side binds all six fixtures; bin 1 binds four because no
#   bin-1 block reads fit4 or fit6 on draws.
#
# Reads dev/brmsport-bin2-verdicts.tsv (written by
# dev/brmsport-bin2-ledger.R); an assertion with no row is generated as
# "unclassified", which only a recording run allows.
source("dev/brmsport-bin2-blocks.R")
source("dev/brmsport-bin2-adapters-spec.R")

out_root <- "dev/brmsport-bin2"
topic_of <- function(f) sub("^tests[.](.*)[.]R$", "\\1", f)

methods_pre <- c(
  "expect_range <- function(object, lower = -Inf, upper = Inf, ...) {",
  "  testthat::expect_true(all(object >= lower & object <= upper), ...)",
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
  "nvisits <- 4")
files <- list(
  "tests.brmsfit-methods.R" = list(
    topic = "methods", pre = methods_pre,
    sample_pre = sprintf("fit%d <- brms_fixture_draws(%d)", 1:6, 1:6)),
  "tests.priors.R" = list(
    topic = "priors",
    pre = c("epilepsy <- brms::epilepsy", "inhaler <- brms::inhaler"))
)
# 2P files ported through an adapter, with the blocks the adapter
# reaches; dev/brmsport-bin2-adapters-spec.R says why each is thin
for (f in names(bin2_adapted)) files[[f]] <- bin2_adapted[[f]]

vf <- "dev/brmsport-bin2-verdicts.tsv"
verdicts <- if (file.exists(vf)) {
  utils::read.delim(vf, stringsAsFactors = FALSE, quote = "")
} else {
  data.frame(id = character(), pkg = character(), verdict = character(),
             reason = character())
}
stopifnot(!anyDuplicated(verdicts[c("id", "pkg")]))
verdict_for <- function(id, pkg) {
  v <- verdicts[verdicts$id == id & verdicts$pkg %in% c(pkg, "*"), ]
  if (nrow(v) == 0L) return(list(verdict = "unclassified", reason = ""))
  if (nrow(v) > 1L) v <- v[v$pkg == pkg, ]
  stopifnot(nrow(v) == 1L)
  list(verdict = v$verdict, reason = v$reason)
}
# the value checks' recorded outcomes (dev/brmsport-bin2-ledger.R); a
# check with none is generated with "", which only a recording run allows
vcf <- "dev/brmsport-bin2-valuechecks.tsv"
vchecks <- if (file.exists(vcf)) {
  utils::read.delim(vcf, stringsAsFactors = FALSE, quote = "",
                    colClasses = "character")
} else {
  data.frame(id = character(), outcome = character())
}
value_expected <- function(id) {
  o <- vchecks$outcome[vchecks$id == id]
  if (length(o) == 1L) o else ""
}
# the message text a manual value verdict was written for, "" for none
value_msg_has <- function(id) {
  o <- vchecks$msg_has[vchecks$id == id]
  if (length(o) == 1L && !is.na(o)) o else ""
}

q <- function(x) encodeString(x, quote = "\"")
emit_reason <- function(reason, indent) {
  if (!nzchar(reason)) return(paste0(indent, "\"\","))
  pieces <- strwrap(reason, width = 60)
  pieces <- paste0(pieces, c(rep(" ", length(pieces) - 1L), ""))
  if (length(pieces) == 1L) return(paste0(indent, q(pieces), ","))
  c(paste0(indent, "paste0("),
    paste0(indent, "  ", q(pieces),
           c(rep(",", length(pieces) - 1L), "),")))
}
code_lines <- function(text, expr, indent, col = 1L) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  if (length(lines) > 1L) {
    lead <- strrep(" ", col - 1L)
    rest <- lines[-1]
    has <- startsWith(rest, lead)
    rest[has] <- substring(rest[has], col)
    lines[-1] <- rest
  }
  lines <- paste0(indent, lines)
  if (all(nchar(lines) <= 80L)) return(lines)
  paste0(indent, deparse(expr, width.cutoff = 60L))
}
apply_subst <- function(text, subst) {
  for (s in subst) text <- gsub(s[1], s[2], text, fixed = TRUE)
  text
}

blocks <- c(brms_bin_blocks("2T"), Filter(function(b) {
  b$file %in% names(bin2_adapted) &&
    b$label %in% bin2_adapted[[b$file]]$blocks
}, brms_bin_blocks("2P")))
out <- list()
ids_by_pkg <- character()
adapted_ids <- character()
value_ids <- character()
for (b in blocks) {
  cfg <- files[[b$file]]
  stopifnot(!is.null(cfg))
  pkgs <- switch(b$tier, core = "frmtmb", sample = "frmtmb.sample",
                 both = c("frmtmb", "frmtmb.sample"))
  if (!is.null(cfg$pkgs)) pkgs <- intersect(pkgs, cfg$pkgs)
  for (pkg in pkgs) {
    key <- paste(pkg, cfg$topic)
    body <- c(sprintf("test_that(%s, {", q(b$label)),
              sprintf("  withr::local_seed(%dL)", b$line))
    for (j in seq_along(b$stmts)) {
      s <- b$stmts[[j]]
      id <- sprintf("%s:%d", topic_of(b$file), b$lines[j])
      text <- apply_subst(b$texts[j], cfg$subst)
      adapted <- !identical(text, b$texts[j])
      expr <- if (adapted) str2lang(text) else s
      # setup lines count too: an adapter mostly rewrites those
      if (adapted) adapted_ids <- c(adapted_ids, id)
      k <- n_expect(s)
      # a value check after the next statement needs the random number
      # state this line starts from (dev/brmsport-bin2-adapters.R)
      vnext <- if (j < length(b$stmts)) bin2_value_check_ok(b, j + 1L)
      vhere <- bin2_value_check_ok(b, j)
      if ((!is.null(vnext) && is.null(vnext$inline)) ||
          !is.null(vhere$inline)) {
        body <- c(body, "  brms_value_mark()")
      }
      if (k == 0L) {
        if (is.call(s) && as.character(s[[1]])[1] %in%
            c("skip_on_cran", "skip_if_not_installed", "context")) {
          body <- c(body, code_lines(text, expr, "  ", b$cols[j]))
        } else if (grepl("brms:::?", text)) {
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
          body <- c(body, sprintf("  brms_setup(%s, {", q(id)),
                    code_lines(text, expr, "    ", b$cols[j]), "  })")
        } else {
          body <- c(body, sprintf("  brms_setup(%s,", q(id)),
                    code_lines(text, expr, "    ", b$cols[j]), "  )")
        }
        next
      }
      ids_by_pkg <- c(ids_by_pkg, paste(pkg, id))
      v <- verdict_for(id, pkg)
      if (grepl("brms:::?", text)) {
        body <- c(body,
                  sprintf("  brms_port_not_run(%s, %s,", q(id),
                          q(v$verdict)),
                  emit_reason(v$reason, "    "),
                  "    \"reaches brms's own namespace\"",
                  "  )",
                  paste0("  # ", strsplit(text, "\n")[[1]]))
        next
      }
      body <- c(body, sprintf("  brms_port(%s, %s,", q(id), q(v$verdict)),
                emit_reason(v$reason, "    "),
                code_lines(text, expr, "    ", b$cols[j]), "  )")
      if (!is.null(vhere)) {
        # the line before (or the call inside the assertion), as brms
        # wrote it, against brms's own function
        value_ids <- c(value_ids, id)
        head <- sprintf("  brms_value_check(%s, %s, %s, {", q(id),
                        q(vhere$var), q(value_expected(id)))
        mh <- value_msg_has(id)
        tail <- if (nzchar(mh)) sprintf("  }, msg_has = %s)", q(mh)) else
          "  })"
        if (is.null(vhere$inline)) {
          body <- c(body, head,
                    code_lines(b$texts[j - 1L], b$stmts[[j - 1L]], "    ",
                               b$cols[j - 1L]),
                    tail)
        } else {
          ref <- paste(".value <-", deparse1(vhere$inline))
          body <- c(body, head, paste0("    ", ref), "  }, adapted = {",
                    paste0("    ", apply_subst(ref, cfg$subst)), tail)
        }
      }
    }
    body <- c(body, "})", "")
    out[[key]] <- c(out[[key]], body)
  }
}
stopifnot(!anyDuplicated(ids_by_pkg))

device_guard <- c(
  "grDevices::pdf(NULL)",
  "withr::defer(grDevices::dev.off(), teardown_env())"
)
header <- function(f, pkg, cfg) {
  pre <- cfg$pre
  if (pkg == "frmtmb.sample" && !is.null(cfg$sample_pre)) {
    pre <- c(pre[!grepl("^fit[0-9] <-", pre)], cfg$sample_pre)
  }
  c("# GENERATED by dev/brmsport-bin2-gen.R from brms 2.23.0",
    sprintf("# tests/testthat/%s, bin 2.", f),
    "# Do not edit by hand: change dev/brmsport-bin2-verdicts.tsv or the",
    "# generator, then regenerate.",
    "# tests/testthat/helper-brms-suite.R says what brms_port() asserts.",
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
  dir <- file.path(out_root, pkg)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(dir, sprintf("test-brms-suite-bin2-%s.R", topic))
  txt <- c(header(f, pkg, files[[f]]), out[[key]])
  while (length(txt) && !nzchar(txt[length(txt)])) txt <- txt[-length(txt)]
  writeLines(txt, path, useBytes = TRUE)
  written <- c(written, path)
}
writeLines(sort(unique(adapted_ids)), file.path(out_root, "adapted-ids.txt"))
ids <- unique(sub("^[^ ]+ ", "", ids_by_pkg))
cat("wrote", length(written), "files;", length(ids), "assertion ids;",
    sum(!ids %in% verdicts$id), "without a verdict;",
    length(unique(adapted_ids)), "adapted;", length(unique(value_ids)),
    "value checks\n")
cat(written, sep = "\n")
