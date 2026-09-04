# Regenerate data/habit_prep.rda from the authors' OSF deposit.
#
# Source: Hardwick RM, Haith AM (2019). "Time-Dependent Competition Between
# Goal-Directed and Habitual Response Selection." OSF, <https://osf.io/3fjez/>,
# MIT licensed, copyright (c) 2019 Robert Michael Hardwick and Adrian Mark Haith.
#
# The deposit's ModelCode/dat.mat is the authors' own preprocessed structure and
# is the input their archived model fits were computed from. Reading it, rather
# than re-deriving the stimulus mappings from the raw .dat files, is what lets
# this data object be checked trial for trial against those archived fits.
#
# Requires R.matlab. Run from the package root:
#   Rscript data-raw/habit_prep.R

stopifnot(requireNamespace("R.matlab", quietly = TRUE))

dest <- Sys.getenv("HABIT_DAT_MAT", "")
if (!nzchar(dest)) {
  dest <- file.path(tempdir(), "dat.mat")
  # OSF resolves a per-file download through the node's storage provider.
  url <- paste0(
    "https://files.osf.io/v1/resources/3fjez/providers/osfstorage/",
    "?zip="
  )
  stop(
    "Set HABIT_DAT_MAT to a local copy of ModelCode/dat.mat from ",
    "https://osf.io/3fjez/ (download the deposit, then point this script at ",
    "ModelCode/dat.mat). Deposit archive: ", url
  )
}

mat <- R.matlab::readMat(dest)
d <- mat$d

# R.matlab represents a MATLAB struct as an n x 1 x 1 array indexed by field name.
pick <- function(x, nm) x[[nm, 1, 1]]
subjects <- function(x) setdiff(dimnames(x)[[1]], c("exclude", "group"))

e1 <- pick(d, "e1")
e2 <- pick(d, "e2")

cells <- c(
  lapply(subjects(e1), function(s) list("minimal", s, pick(pick(e1, s), "untrained"))),
  lapply(subjects(e1), function(s) list("4day", s, pick(pick(e1, s), "trained"))),
  lapply(subjects(e2), function(s) list("20day", s, pick(e2, s)))
)

# Columns of the trial matrices, from the header of the raw .dat files.
COL <- c(subject = 1, session = 2, stamp = 3, trial = 4, des_onset = 5,
         rec_onset = 6, stimulus = 7, target_key = 8, response_num = 9,
         response_key = 10, rt = 11, correct = 12, tally = 13, spare = 14,
         prep_time = 15)
RESP <- c("other", "correct", "habit")  # modelCoded.y is 0, 1, 2

out <- lapply(cells, function(cell) {
  group <- cell[[1]]; sid <- cell[[2]]; node <- cell[[3]]
  tr <- pick(node, "tr")
  rv <- pick(tr, "revisedTrials")
  un <- pick(tr, "unchangedTrials")
  mc <- pick(tr, "modelCoded")
  y <- as.vector(pick(mc, "y"))
  x <- as.vector(pick(mc, "x"))

  all <- rbind(rv, un)
  remapped <- c(rep(TRUE, nrow(rv)), rep(FALSE, nrow(un)))
  # Trial index, not prep-time value, decides which stimulus class a trial
  # belongs to; ordering by it restores the sequence modelCoded was built in.
  o <- order(all[, COL[["trial"]]])
  all <- all[o, , drop = FALSE]
  remapped <- remapped[o]
  stopifnot(nrow(all) == length(y),
            isTRUE(all.equal(all[, COL[["prep_time"]]], x, tolerance = 1e-12)))

  exp_id <- if (group == "20day") "e2" else "e1"
  data.frame(
    participant  = paste0(exp_id, "-", formatC(as.integer(sub("^s", "", sid)),
                                               width = 2, flag = "0")),
    group        = group,
    subject_code = as.integer(all[, COL[["subject"]]]),
    trial        = as.integer(all[, COL[["trial"]]]),
    prep_time    = as.numeric(all[, COL[["prep_time"]]]),
    remapped     = remapped,
    stimulus     = as.integer(all[, COL[["stimulus"]]]),
    target_key   = as.integer(all[, COL[["target_key"]]]),
    response_key = as.integer(all[, COL[["response_key"]]]),
    response     = RESP[as.integer(y) + 1L],
    stringsAsFactors = FALSE
  )
})

habit_prep <- do.call(rbind, out)

# The documented cleaning rule rejects a recorded key outside the four response
# keys. Those rows are hardware mis-reads that the upstream pipeline scores as
# ordinary "other" errors, so the category is withheld rather than trusted.
bad <- habit_prep$response_key < 1L | habit_prep$response_key > 4L
habit_prep$response[bad] <- NA_character_

habit_prep$participant <- factor(habit_prep$participant)
habit_prep$group <- factor(habit_prep$group, levels = c("minimal", "4day", "20day"))
habit_prep$response <- factor(habit_prep$response, levels = c("correct", "habit", "other"))
rownames(habit_prep) <- NULL

stopifnot(nrow(habit_prep) == 28940L, sum(bad) == 8L,
          nlevels(habit_prep$participant) == 36L)

save(habit_prep, file = "data/habit_prep.rda", compress = "xz")
cat("rows:", nrow(habit_prep),
    " size:", file.size("data/habit_prep.rda"), "bytes\n")
