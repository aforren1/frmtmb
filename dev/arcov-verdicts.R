# The verdict moves lane arcov makes to the ported brms tier, applied
# idempotently so the files can be rebuilt from the base tree.
#
#   Rscript dev/arcov-verdicts.R
#
# WHY NOT THE PIPELINE. dev/brmsport-ledger.R and dev/brmsport-gen.R read
# brms 2.23.0's test sources from dev/brms-suite/, which this machine
# does not have (the directory is not in the repository, the tarball is
# checked by sha256, CRAN is unreachable and the brms repository has no
# v2.23.0 tag). So this script applies to each file exactly the row the
# pipeline would write, following dev/brmsport-ledger.R line by line,
# and the generated test files were edited by hand to what
# dev/brmsport-gen.R's emit_reason() emits (dev/arcov-findings.md says
# which rows). The recorded runs are dev/arcov-log/rec-brm.tsv and
# dev/arcov-log/rec-standata.tsv, from the lane's own build.
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
wr <- function(x, f) utils::write.table(x, f, sep = "\t", quote = FALSE,
                                        row.names = FALSE)
man_f <- "dev/brmsport-verdicts-manual.tsv"
own_f <- "dev/brmsport-verdicts-own.tsv"
ver_f <- "dev/brmsport-verdicts.tsv"
led_f <- "dev/brmsport-ledger.tsv"
man <- rd(man_f)
own <- rd(own_f)
ver <- rd(ver_f)
led <- rd(led_f)

# brm:110, brm(y ~ ma(x), dat, poisson()): frmtmb now fits ma() without
# cov = TRUE for gaussian and student, and its refusal for poisson
# quotes brms's own sentence ("Please set cov = TRUE when modeling MA
# structures for this family"), so brms's assertion holds as written.
# The harness reported "STALE OWN-WORDS: brms's assertion holds as
# written"; ledger.R's rule for a held row with no manual verdict and
# no note is outcome pass, class "", reason "".
own <- own[own$id != "brm:110", , drop = FALSE]

# brm:106 and brm:108 still hold in frmtmb's own words; their notes said
# the grammar check runs "before the cov = TRUE one", and there is no
# cov = TRUE check any more
note_new <- c(
  "brm:106" = paste0("brms: Cannot coerce 'x + y' to a single variable ",
                     "name. The same refusal of the same call in ",
                     "frmtmb's words (dev/adefects-log/evidence.txt)"),
  "brm:108" = paste0("brms: Illegal grouping term 'g1/g2'. The same ",
                     "refusal of the same call in frmtmb's words"))
for (id in names(note_new)) {
  i <- which(own$id == id)
  stopifnot(length(i) == 1L)
  own$note[i] <- note_new[[id]]
}

# standata:310, :315, :316 still cannot transfer, for the reason they
# always gave (brms's old_order attribute); the setup lines above them
# now RUN, so the "also refused without cov = TRUE" clause is false
reason_new <- c(
  "standata:310" = paste0(
    "reads brms's internal old_order attribute of Stan data, which ",
    "exists because brms sorts rows for its autocorrelation code; ",
    "frmtmb keeps the data order and has no such attribute. The setup ",
    "ar(time, id) runs (brms's cov = FALSE form, ?frmtmb-autocor)"),
  "standata:315" = paste0(
    "reads brms's internal old_order attribute of Stan data; frmtmb ",
    "keeps the data order and has no such attribute. The setup ",
    "ma(time, id) runs (?frmtmb-autocor)"),
  "standata:316" = paste0(
    "reads brms's internal old_order attribute of Stan data; frmtmb ",
    "keeps the data order and has no such attribute. The setup ",
    "ma(time, id) runs (?frmtmb-autocor)"))
for (id in names(reason_new)) {
  i <- which(man$id == id)
  stopifnot(length(i) == 1L)
  man$reason[i] <- reason_new[[id]]
}

# the pipeline's outputs, row by row as ledger.R writes them
rec <- do.call(rbind, lapply(
  c("dev/arcov-log/rec-brm.tsv", "dev/arcov-log/rec-standata.tsv"),
  function(f) {
    utils::read.delim(f, header = FALSE, quote = "", colClasses = "character",
                      na.strings = NULL,
                      col.names = c("kind", "pkg", "id", "verdict", "held",
                                    "vacuous", "msg", "caught", "raw_held"))
  }))
rec <- rec[rec$kind == "assert", ]
msg_of <- function(id) {
  m <- rec$msg[rec$id == id]
  stopifnot(length(m) == 1L)
  m
}
set_ver <- function(id, verdict, reason) {
  i <- which(ver$id == id)
  stopifnot(length(i) == 1L)
  ver$verdict[i] <<- verdict
  ver$reason[i] <<- reason
}
set_led <- function(id, outcome, class, reason, held, message) {
  i <- which(led$id == id)
  stopifnot(length(i) == 1L)
  led$outcome[i] <<- outcome
  led$class[i] <<- class
  led$reason[i] <<- reason
  led$held[i] <<- held
  led$message[i] <<- message
}
set_ver("brm:110", "pass", "")
set_led("brm:110", "pass", "", "", "core=TRUE", msg_of("brm:110"))
for (id in names(note_new)) {
  o <- own[own$id == id, ]
  r <- paste("frmtmb's own words:", o$pattern, "|", o$note)
  set_ver(id, "pass", r)
  set_led(id, "pass", "own-words", r, "core=TRUE", msg_of(id))
}
for (id in names(reason_new)) {
  set_ver(id, "cannot transfer", reason_new[[id]])
  set_led(id, "cannot transfer", "stan", reason_new[[id]], "core=FALSE",
          msg_of(id))
}

wr(man, man_f)
wr(own, own_f)
wr(ver, ver_f)
wr(led, led_f)
cat("manual", nrow(man), "own", nrow(own), "verdicts", nrow(ver),
    "ledger", nrow(led), "\n")

# the generated summary's two class counts that the move changes,
# recounted from the ledger as dev/brmsport-ledger.R counts them
sm_f <- "dev/brmsport-log/ledger-summary.md"
sm <- readLines(sm_f)
n_pass <- sum(led$outcome == "pass" & led$class == "")
n_own <- sum(led$outcome == "pass" & led$class == "own-words")
sm <- sub("^[|] pass [|] - [|] [0-9]+ [|]$",
          sprintf("| pass | - | %d |", n_pass), sm)
sm <- sub("^[|] pass [|] own-words [|] [0-9]+ [|]$",
          sprintf("| pass | own-words | %d |", n_own), sm)
writeLines(sm, sm_f)

# the tier's own recordings, for exactly the rows this lane's build
# changes: the three assertions above and the setup lines of the
# standata block, which now run
ids <- c("brm:106", "brm:108", "brm:110", "standata:309", "standata:310",
         "standata:313", "standata:314", "standata:315", "standata:316")
for (t in c("brm", "standata")) {
  f <- sprintf("dev/brmsport-log/rec-frmtmb-%s.tsv", t)
  old <- readLines(f)
  new <- readLines(sprintf("dev/arcov-log/rec-%s.tsv", t))
  id_of <- function(x) vapply(strsplit(x, "\t", fixed = TRUE),
                              function(p) p[3L], "")
  hit <- id_of(old) %in% ids
  repl <- new[match(id_of(old)[hit], id_of(new))]
  stopifnot(!anyNA(repl))
  old[hit] <- repl
  writeLines(old, f)
}
