# Which recorded passes of the ported brms suite can rest on the harness
# gap that a setup line reading a stale name comes out fresh
# (dev/remeasure-bin2-stale-setup.R constructs the gap).
#
#   Rscript dev/remeasure-stale-reach.R
#
# For every block of bins 1 and 2: a setup line that failed or was not
# run makes its assigned names stale; a LATER setup line that reads one
# of them and succeeded passes the staleness on to the names it assigns
# (brms's block would have stopped first); an assertion that reads a
# name tainted this way and was recorded as holding is listed. Static
# over brms's parsed source and the recorded outcomes in
# dev/brmsport-log/rec-*.tsv (bin 1) and dev/brmsport-bin2-log/rec-*.tsv.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
source("dev/brmsport-bin2-blocks.R")
topic_of <- function(f) sub("^tests[.](.*)[.]R$", "\\1", f)
cols <- c("kind", "pkg", "id", "verdict", "held", "vacuous", "msg",
          "caught", "raw_held")
rd <- function(dir) do.call(rbind, lapply(
  list.files(dir, "^rec-.*[.]tsv$", full.names = TRUE), function(f) {
    utils::read.delim(f, header = FALSE, quote = "", col.names = cols,
                      colClasses = "character", na.strings = NULL)
  }))
rec <- rbind(rd("dev/brmsport-log"), rd("dev/brmsport-bin2-log"))
assigned <- function(e) {
  if (!is.call(e)) return(character())
  out <- character()
  if (as.character(e[[1]])[1] %in% c("<-", "=") && length(e) == 3L) {
    t <- e[[2]]
    while (is.call(t)) t <- t[[2]]
    if (is.name(t)) out <- as.character(t)
  }
  args <- as.list(e)[-1]
  for (i in seq_along(args)) {
    # the empty argument of x[i, ] cannot be bound to a name
    if (is.name(args[[i]]) && !nzchar(as.character(args[[i]]))) next
    out <- c(out, assigned(args[[i]]))
  }
  unique(out)
}
scan_block <- function(b, rec, p) {
  hits <- list()
  tainted <- character()
  for (j in seq_along(b$stmts)) {
    s <- b$stmts[[j]]
    id <- sprintf("%s:%d", topic_of(b$file), b$lines[j])
    r <- rec[rec$id == id & rec$pkg == p, ]
    if (!nrow(r)) next
    reads <- setdiff(all.vars(s), assigned(s))
    if (r$kind[1] == "setup") {
      if (r$held[1] != "TRUE") {
        tainted <- union(tainted, assigned(s))
      } else if (any(reads %in% tainted)) {
        tainted <- union(tainted, assigned(s))
      } else {
        tainted <- setdiff(tainted, assigned(s))
      }
    } else if (r$held[1] == "TRUE" && any(all.vars(s) %in% tainted)) {
      hits[[length(hits) + 1L]] <- data.frame(
        pkg = p, id = id, reads = paste(intersect(all.vars(s), tainted),
                                        collapse = " "))
    }
  }
  do.call(rbind, hits)
}

# The instrument's own cases, before it is believed: the construction of
# dev/remeasure-bin2-stale-setup.R must be found, and the same block with
# the failing line succeeding must not.
ex <- as.list(parse(text = c("a <- 1", "a <- stop('x')", "b <- a + 1",
                             "expect_equal(b, 2)")))
fake <- list(file = "tests.fake.R", stmts = ex, lines = 1:4)
mk <- function(held2) data.frame(
  kind = c("setup", "setup", "setup", "assert"), pkg = "frmtmb",
  id = paste0("fake:", 1:4), held = c("TRUE", held2, "TRUE", "TRUE"))
pos <- scan_block(fake, mk("FALSE"), "frmtmb")
neg <- scan_block(fake, mk("TRUE"), "frmtmb")
cat("self-test: the constructed gap found", NROW(pos), "time(s); with the",
    "failing line repaired,", NROW(neg), "\n")
stopifnot(NROW(pos) == 1L, NROW(neg) == 0L)

hits <- list()
for (b in c(brms_bin1_blocks(), brms_bin_blocks(c("2T", "2P")))) {
  for (p in unique(rec$pkg)) {
    hits[[length(hits) + 1L]] <- scan_block(b, rec, p)
  }
}
H <- do.call(rbind, hits)
cat("blocks scanned:", length(c(brms_bin1_blocks(),
                                 brms_bin_blocks(c("2T", "2P")))),
    "; recorded rows:", nrow(rec), "\n")
cat("held assertions that read a name tainted through a setup line:",
    if (is.null(H)) 0L else nrow(H), "\n")
if (!is.null(H)) print(H, row.names = FALSE)
