# Join dev/famlink-port-ledger-{base,lane}.tsv and classify every brms
# assertion: genuine pass, deliberate divergence (with the reason),
# cannot transfer (with the error that shows it), or not fixed. Prints a
# markdown table for dev/famlink-findings.md, verbatim.
rd <- function(arm) {
  utils::read.delim(paste0("dev/famlink-port-ledger-", arm, ".tsv"),
                    colClasses = "character", quote = "")
}
b <- rd("base")
l <- rd("lane")
stopifnot(identical(b$assertion, l$assertion))

# A pass is genuine only when the assertion reached frmtmb code. An
# expect_error() that caught "could not find function" or "object not
# found" passed on R's own error, and the ledger script flags it.
genuine <- function(d) d$pass == "TRUE" & d$vacuous != "TRUE"
gb <- genuine(b)
gl <- genuine(l)

# An assertion that reached no frmtmb code, pass or fail, is "cannot
# transfer", decided from the error it caught and not from a list of
# names. The rules classify the remaining lane failures, first match
# wins, and a lane failure no rule covers stops the script.
rules <- list(
  list("acat\\(cloglog\\)",
       paste("deliberate divergence: acat off the logit is a second",
             "density frmtmb has not written (?frmtmb-links)")),
  list("Cannot coerce 'alink'",
       paste("deliberate divergence: the call is refused, but brms's",
             "message names its internal variable `alink`")),
  list("family_names\\(mix\\)",
       paste("not fixed: mixture() has no `nmix`, and",
             "brms:::family_names() is a brms internal")),
  list("order = \"x\"",
       paste("not fixed: mixture() has no `order`; the call is refused",
             "naming `order`, not with brms's sentence"))
)
cls <- character(nrow(l))
for (i in seq_len(nrow(l))) {
  if (gl[i]) {
    cls[i] <- if (gb[i]) "pass (genuine on base too)" else
      "pass (not a genuine pass on base)"
    next
  }
  if (l$vacuous[i] == "TRUE") {
    why <- if (nzchar(l$caught[i])) l$caught[i] else l$detail[i]
    cls[i] <- paste0("cannot transfer: ", sub("^ERROR: ", "", why))
    next
  }
  hit <- Filter(function(r) grepl(r[[1]], l$assertion[i]), rules)
  if (!length(hit)) stop("unclassified lane failure: ", l$assertion[i])
  cls[i] <- hit[[1]][[2]]
}
esc <- function(x) gsub("|", "\\|", x, fixed = TRUE)
short <- function(x) {
  ifelse(nchar(x) > 70, paste0(substr(x, 1, 67), "..."), x)
}
cat("---- GENERATED: dev/famlink-ledger-sum.R ----\n")
cat(sprintf("brms assertions: %d (%s)\n", nrow(l),
            paste(names(table(l$file)), table(l$file), collapse = ", ")))
cat(sprintf("genuine passes on base: %d; on lane: %d\n", sum(gb), sum(gl)))
cat(sprintf(paste("passes on a missing function or object, not counted:",
                  "base %d, lane %d\n"),
            sum(b$pass == "TRUE" & b$vacuous == "TRUE"),
            sum(l$pass == "TRUE" & l$vacuous == "TRUE")))
tab <- table(sub(":.*", "", cls))
cat(paste0(names(tab), ": ", tab, collapse = "\n"), "\n\n")
cat("| # | file | assertion | outcome |\n|---|---|---|---|\n")
cat(sprintf("| %s | %s | `%s` | %s |\n", l$idx, sub("^tests[.]", "", l$file),
            esc(short(l$assertion)), esc(cls)), sep = "")
cat("---- END GENERATED ----\n")
