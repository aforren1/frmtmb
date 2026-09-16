# Spot check: run a handful of bin-1 blocks from brms's own suite
# against frmtmb, one ASSERTION at a time.
#
# Why one at a time: a block stops at its first failure, so a whole-block
# run reports one number and hides the other eighty. Every top-level
# expression inside the block is evaluated in order in one environment,
# and its outcome is recorded, so a setup line that dies is visible as
# such and the assertions after it still run.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "\n")

SUITE <- "dev/brms-suite/brms/tests/testthat"

blocks_of <- function(file) {
  ex <- parse(file.path(SUITE, file), keep.source = FALSE)
  out <- list()
  for (i in seq_along(ex)) {
    e <- ex[[i]]
    if (!is.call(e) || !identical(e[[1]], as.name("test_that"))) next
    out[[length(out) + 1L]] <- list(
      desc = paste(as.character(e[[2]]), collapse = " "),
      body = e[[3]]
    )
  }
  out
}

# The outcome that matters is not pass or fail. An expect_error() whose
# regexp does not match re-raises the original error, so a plain
# pass/fail split would score "frmtmb refuses this, in other words" and
# "frmtmb ACCEPTS what brms refuses" identically. These are the four
# outcomes the audit needs to tell apart.
classify <- function(cnd, txt) {
  msg <- conditionMessage(cnd)
  # "could not find function" means frmtmb has no such family or helper.
  # "object 'logit' not found" means something else entirely: brms
  # deparses an UNQUOTED link name and frmtmb evaluates it. That is a
  # divergence, not an absence, so the two are counted apart.
  absent <- grepl("could not find function|there is no package", msg)
  unbound <- grepl("object '[^']*' not found|unused argument", msg)
  if (inherits(cnd, "expectation_failure")) {
    if (grepl("to throw", msg)) return("NO_CONDITION")
    if (grepl("does not match|Actual message", msg)) return("MSG_DIFFERS")
    return("FAIL")
  }
  if (absent) return("ABSENT")
  if (unbound) return("UNBOUND")
  if (grepl("^expect_(error|warning|message)", txt)) return("MSG_DIFFERS")
  "ERROR"
}

results <- list()
run_block <- function(file, blk, env_extra = list()) {
  env <- new.env(parent = globalenv())
  for (nm in names(env_extra)) assign(nm, env_extra[[nm]], envir = env)
  body <- blk$body
  stmts <- if (is.call(body) && identical(body[[1]], as.name("{"))) {
    as.list(body)[-1]
  } else list(body)
  for (s in stmts) {
    txt <- paste(deparse(s), collapse = " ")
    txt <- gsub("[[:space:]]+", " ", txt)
    is_assert <- grepl("^expect", txt)
    out <- tryCatch({
      withCallingHandlers(eval(s, envir = env),
                          warning = function(w) invokeRestart("muffleWarning"),
                          message = function(m) invokeRestart("muffleMessage"))
      "PASS"
    }, condition = function(cnd) {
      if (inherits(cnd, "error") || inherits(cnd, "expectation_failure")) {
        attr_msg <<- conditionMessage(cnd)
        classify(cnd, txt)
      } else "PASS"
    })
    msg <- if (identical(out, "PASS")) "" else attr_msg
    results[[length(results) + 1L]] <<- data.frame(
      file = file, block = blk$desc, assertion = is_assert,
      code = substr(txt, 1, 110), outcome = out,
      message = substr(gsub("[\r\n\t]+", " ", msg), 1, 200),
      stringsAsFactors = FALSE
    )
  }
  invisible(NULL)
}

attr_msg <- ""

# --- group 1: family and link validation, verbatim ---------------------
for (blk in blocks_of("tests.families.R")) run_block("tests.families.R", blk)

# --- group 2: formula grammar, verbatim --------------------------------
for (blk in blocks_of("tests.brmsformula.R")) {
  run_block("tests.brmsformula.R", blk)
}

# --- group 3: brm()'s validation errors, with brm mapped to frm --------
# The only adaptation is the entry point's name; every argument, every
# formula and every expected message is brms's.
shim <- list(brm = function(formula, data, ...) {
  frmtmb::frm(formula, data, ..., dry_run = "frame")
})
for (blk in blocks_of("tests.brm.R")) {
  if (!grepl("expected errors", blk$desc)) next
  run_block("tests.brm.R", blk, env_extra = shim)
}

# --- group 4: response validation, standata mapped to dry_run ----------
shim2 <- list(standata = function(formula, data, family = NULL, ...) {
  frmtmb::frm(formula, data, family = family, ..., dry_run = "frame")
})
for (blk in blocks_of("tests.standata.R")) {
  if (!grepl("rejects incorrect response|suggests using family bernoulli|rejects incorrect addition",
             blk$desc)) next
  run_block("tests.standata.R", blk, env_extra = shim2)
}

tab <- do.call(rbind, results)
write.table(tab, "dev/brmssuite-spotcheck.tsv", sep = "\t",
            row.names = FALSE, quote = TRUE)

cat("\n== assertions only ==\n")
a <- tab[tab$assertion, ]
print(table(a$file, a$outcome))
cat("\n== setup lines (non-assertions) that died ==\n")
s <- tab[!tab$assertion & tab$outcome != "PASS", ]
if (nrow(s)) print(s[, c("file", "code", "outcome")]) else cat("none\n")
cat("\n== UNBOUND detail: what symbol failed to resolve ==\n")
u <- a[a$outcome == "UNBOUND", ]
print(table(sub(".*(object '[^']*' not found).*", "\\1", u$message)))
cat("\n== assertions where frmtmb raised NO condition at all ==\n")
bad <- a[a$outcome %in% c("NO_CONDITION", "FAIL", "ERROR"),
         c("file", "code", "outcome", "message")]
for (i in seq_len(nrow(bad))) {
  cat(sprintf("[%s] %s\n    %s :: %s\n", bad$outcome[i], bad$file[i],
              bad$code[i], bad$message[i]))
}
cat("\n== assertions where frmtmb refused, in other words ==\n")
d <- a[a$outcome == "MSG_DIFFERS", c("code", "message")]
for (i in seq_len(nrow(d))) {
  cat(sprintf("  %s\n    frmtmb: %s\n", d$code[i], d$message[i]))
}
cat("\nTOTAL assertions:", nrow(a), " PASS:", sum(a$outcome == "PASS"),
    " FAIL:", sum(a$outcome == "FAIL"),
    " ERROR:", sum(a$outcome == "ERROR"), "\n")
