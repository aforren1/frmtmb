#' Split a NONMEM-shaped dosing table into observations and events
#'
#' A population pharmacokinetic dataset usually comes as one table in the
#' layout NONMEM defined and rxode2 and nlmixr2 read: one row per record,
#' an `evid` column that says whether the row is an observation or a
#' dose, and `amt`, `cmt`, `rate`, `ii`, `addl` and `ss` on the dose
#' rows. [frm_ode()] and [frm_lincmt()] take the observations and the
#' doses as two tables. `frm_ode_records()` does that split and nothing
#' else: it does no arithmetic on a dose, so a dose means in the result
#' exactly what it meant in the record.
#'
#' @section What each record becomes:
#' \tabular{ll}{
#'   `evid` \tab result \cr
#'   0 \tab a row of `data` \cr
#'   1 \tab an `"add"` event: `value = amt` into `state = cmt` \cr
#'   3 \tab a `"reset"` event: every compartment set to zero \cr
#'   4 \tab a `"reset"` event and then an `"add"` event at the same time \cr
#'   5 \tab a `"replace"` event, rxode2's code: the compartment set to `amt` \cr
#'   6 \tab a `"multiply"` event, rxode2's code: the compartment scaled by `amt`
#' }
#'
#' A dose with a positive `rate` is an infusion of duration
#' `amt / rate`, and a positive `dur`, which is how rxode2 writes one,
#' is that duration directly. `ii`, `addl` and `ss = 1` are copied to
#' the event as they are; [frm_ode()] reads them the way NONMEM does.
#'
#' `data` keeps every column of `d` except the dosing items (`evid`,
#' `amt`, `rate`, `dur`, `ii`, `addl`, `ss` and `mdv`), so `id`, `time`,
#' `cmt`, the response and every covariate come through under their own
#' names, in their original row order.
#'
#' A dose item is read ONLY through its argument. Every other column is
#' data: a rate held in a column called `RATE_MG_H` and not named with
#' `rate = "RATE_MG_H"` is a covariate here, and the doses are read
#' without it. This function cannot recognize every name a dataset
#' might give an item, so name each one that is not under its default
#' name.
#'
#' @section What is refused, and why:
#' Every refusal names the column or the record, because each one is a
#' place where a record would otherwise be read with a meaning different
#' from its own.
#'
#' Record types and items:
#'
#' * `evid = 2`, NONMEM's "other event". It carries no dose and no
#'   observation, and is used to change a time-varying covariate.
#'   Dropping it would lose that change, and keeping it would make it an
#'   observation. Put the covariate on the observation rows, or into
#'   [frm_ode()]'s `tv`.
#' * Any other `evid` value, including rxode2's codes above 6 and its
#'   classic codes of 100 and more, such as the 101 in
#'   `nlmixr2data::theo_sd`, which fold the compartment into the number.
#' * `ss = 2`, NONMEM's steady state added to the current state rather
#'   than replacing it. [frm_ode()]'s steady state starts from an empty
#'   system, which is `ss = 1`.
#' * A steady-state constant infusion (`ss = 1`, `amt = 0`, `rate > 0`).
#' * A negative `rate` or `dur`: NONMEM's `-1` and `-2` say the rate or
#'   the duration is a model parameter, and an estimated duration moves
#'   where the solve is split, which [frm_ode()] cannot do.
#' * A dose with both a positive `rate` and a positive `dur`.
#' * A `rate`, `dur`, `ii`, `addl` or `ss` that is not finite on an
#'   event, or that is nonzero on an observation; an `amt`, `ii`,
#'   `addl` or `ss` that is nonzero on a reset (`evid = 3`). The split
#'   would drop each of them without a word.
#' * A `cmt` on a dose that is missing, zero or below, or not a whole
#'   number.
#' * An `mdv` column that disagrees with `evid`. `mdv = 1` on an
#'   `evid = 0` row is a missing observation, which would otherwise
#'   become a fitted one.
#'
#' Columns:
#'
#' * The NONMEM items that change how a record is read and that this
#'   function does not interpret: `date`, `dat1`, `dat2`, `dat3`,
#'   `cont`, `call`, `pcmt`, `l1` and `l2`.
#' * The names NONMEM's PREDPP or rxode2 read as a dose modifier, which
#'   [frm_ode()] does not apply, so as data they would change nothing:
#'   `ALAGn`, `Rn`, `Dn` and `MTIMEn` (any number `n`), `Fn` for a
#'   compartment `n` of 1 or more, `XSCALE`, `TSCALE`, and `lag`,
#'   `alag`, `tlag` and `tinf` with or without a number. A covariate with
#'   one of these names has to be renamed. `Sn`, `FO` and `F0` are not
#'   refused: `F0` and `FO` are two spellings of PREDPP's output
#'   fraction, and like the scale `Sn` it acts on an observation, not a
#'   dose.
#' * `CENS`, `LIMIT`, `BLQ` and `LLOQ`, which rxode2 and nlmixr2 read as
#'   censoring of an observation. This function does not read censoring,
#'   so a censored value would be fitted as an observed one.
#' * A column argument that names no column of `d`, or two columns of
#'   `d` that differ only in case, or an argument set to `NULL` while `d`
#'   has a column under that item's default name.
#'
#' Order:
#'
#' * An `id` whose records are not one contiguous block, and a `time`
#'   that goes backwards within an `id`. NM-TRAN reads one individual as
#'   one block in time order; here the two blocks would be merged and
#'   dosed twice.
#' * A reset (`evid` 3 or 4), a steady-state dose (`ss = 1`), or a
#'   replace or multiply (`evid` 5 or 6) at the same time as another
#'   event of the same `id`, `addl` repeats included. Two times are one
#'   instant when they differ by at most `64 * .Machine$double.eps` times
#'   the larger of `|t0|` and `k * ii` for a repeat `t0 + k * ii` (`|t|`
#'   for a written time), which covers the rounding of that sum with a
#'   32-fold margin. So a repeat at `0 + 3 * 0.1` meets a record at
#'   `0.3`, while records 1 ms apart at `t = 1e7` or 1 s apart at an
#'   epoch time near `1.7e9` stay distinct. A replace or multiply on a
#'   DIFFERENT compartment from the
#'   other event is refused too, although the two commute: the rule
#'   does not look at compartments. NONMEM reads
#'   records at one time in the order they are listed; [frm_ode()]
#'   applies a reset or a steady state first whatever the listing; and
#'   rxode2 5.1.7 follows neither. Measured on a bolus of 100 into the
#'   central compartment at `t = 10` and a steady-state dose into the
#'   depot at the same time, the central amount at `t = 11` is 112.03 in
#'   [frm_ode()] and 32.58 in rxode2, in either listed order. Two plain
#'   doses at one time add, and pass.
#' * A reset while an infusion of the same `id` is running, including an
#'   infusion that starts at the reset's own time and `addl` repeats.
#'   [frm_ode()] lets the infusion run on after the reset,
#'   [frm_lincmt()] refuses it, and rxode2 5.1.7 returns negative
#'   amounts after it.
#' * An observation listed after an event at the same time for the same
#'   `id`. NONMEM reads such an observation after the event;
#'   [frm_ode()] reads every observation at an event time before the
#'   event. An observation listed before a dose is the pre-dose sample
#'   in both and passes.
#' * An observation at the time of a steady-state record, in either
#'   order. [frm_ode()] reads the steady-state trough there and rxode2
#'   the state after the dose. What NONMEM reads there was not checked,
#'   because no NONMEM was available, so no reading is safe.
#'
#' rxode2 does not follow the listed order at a dose either: measured
#' with rxode2 5.1.7, an observation at a dose time reads the state
#' AFTER the dose whether it is listed before or after it. So a table
#' written for rxode2 with an observation at the exact time of a dose
#' into the observed compartment gives a different value here, where it
#' is the trough. Nothing in the records can tell the two readings
#' apart.
#'
#' Doses of an `id` that has no observation are dropped, with a warning,
#' because [frm_ode()] refuses a dose for a group it does not see. Their
#' ids are the `"dropped_ids"` attribute of the result.
#'
#' @param d A data.frame of NONMEM-shaped records.
#' @param id,time,evid,amt,cmt,rate,dur,ii,addl,ss The name of the
#'   column of `d` that holds each record item. Names match `d`'s column
#'   names without regard to case, so the default `"amt"` finds a column
#'   `AMT`. `NULL` says `d` has no such column, and is refused when `d`
#'   has one under the default name. `time` and `evid` must
#'   exist; `amt` must exist when there is a dose. Without `id` the whole
#'   table is one subject. Without `cmt` the events carry no `state`
#'   column, which [frm_ode()] accepts only for a system of one state.
#'
#' @return A list with `data`, the observation rows, and `events`, the
#'   dosing table for [frm_ode()] or [frm_lincmt()]'s `events` argument,
#'   or `NULL` when there is no dose. Pass `group = <id column>` to the
#'   solver, and `states` when `cmt` holds names.
#'
#' @seealso [frm_ode()] for the `events` grammar these records become.
#'
#' @examples
#' # Two subjects, 100 into the depot every 12 hours for four doses,
#' # sampled at 2, 6 and 42 hours; the second subject's first sample is
#' # missing.
#' rec <- data.frame(
#'   ID = rep(1:2, each = 4),
#'   TIME = rep(c(0, 2, 6, 42), 2),
#'   EVID = rep(c(1, 0, 0, 0), 2),
#'   MDV = c(1, 0, 0, 0, 1, 1, 0, 0),
#'   AMT = rep(c(100, 0, 0, 0), 2),
#'   CMT = rep(c(1, 2, 2, 2), 2),
#'   II = rep(c(12, 0, 0, 0), 2),
#'   ADDL = rep(c(3, 0, 0, 0), 2),
#'   DV = c(NA, 7.1, 8.4, 3.2, NA, NA, 9.0, 3.9))
#' rec <- rec[!(rec$EVID == 0 & rec$MDV == 1), ]
#' x <- frm_ode_records(rec)
#' x$data
#' x$events
#' frm_lincmt(parms = list(ka = 1, ke = 0.2, V = 10), times = x$data$TIME,
#'            group = x$data$ID, ncmt = 1, depot = TRUE,
#'            events = x$events)
#' @export
frm_ode_records <- function(d, id = "id", time = "time", evid = "evid",
                            amt = "amt", cmt = "cmt", rate = "rate",
                            dur = "dur", ii = "ii", addl = "addl",
                            ss = "ss") {
  if (!is.data.frame(d)) {
    frm_stop("frm_ode_records() reads a data.frame of NONMEM-shaped ",
             "records, and `d` is ", class(d)[1L], call. = FALSE)
  }
  nms <- names(d)
  low <- tolower(nms)
  items <- list(id = id, time = time, evid = evid, amt = amt, cmt = cmt,
                rate = rate, dur = dur, ii = ii, addl = addl, ss = ss)
  col <- vapply(names(items), function(arg) {
    ode_rec_column(items[[arg]], arg, nms, low)
  }, "")
  # NULL says `d` has no such column. When it does have one under the
  # default name, NULL would silently read an infusion as a bolus, drop
  # the addl doses or the steady state, or merge the subjects.
  gone <- names(items)[vapply(items, is.null, TRUE) &
                         names(items) %in% low]
  if (length(gone)) {
    frm_stop("frm_ode_records(): `", gone[1L], " = NULL` says `d` has no ",
             "such column, but `d` has the column ",
             nms[match(gone[1L], low)], ". NULL would drop what that column ",
             "says about the records. Name it with `", gone[1L], " = \"",
             nms[match(gone[1L], low)], "\"`, or remove the column from `d`",
             call. = FALSE)
  }
  for (req in c("time", "evid")) {
    if (is.na(col[[req]])) {
      frm_stop("frm_ode_records(): `d` has no `", req, "` column. ",
               if (req == "evid") {
                 paste0("The record type decides which rows are doses, ",
                        "and it is not inferred from `amt` or `mdv`. ")
               },
               "Name the column with `", req, " = `. The columns of `d` ",
               "are: ", paste(nms, collapse = ", "), call. = FALSE)
    }
  }
  used <- col[!is.na(col)]
  # Items that change how a record is read. Passing them through as data
  # would drop that meaning without a word.
  reserved <- c("date", "dat1", "dat2", "dat3", "cont", "call", "pcmt",
                "l1", "l2")
  hit <- setdiff(nms[low %in% reserved], used)
  if (length(hit)) {
    frm_stop("frm_ode_records(): `d` has the NONMEM data item",
             if (length(hit) > 1L) "s " else " ",
             paste(hit, collapse = ", "),
             ", which change how a record is read and which this function ",
             "does not interpret. Resolve ", if (length(hit) > 1L) "them" else
               "it", " into `time` and plain columns first, or drop ",
             if (length(hit) > 1L) "them" else "it",
             call. = FALSE)
  }
  # Names that would change a dose if honored: PREDPP's additional PK
  # parameters for a compartment's lag, bioavailability, rate and
  # duration, its time scale and model event times, and the spellings
  # rxode2 and common datasets use for the same things. frm_ode() reads
  # none of them, so as data they would be a dose modifier dropped in
  # silence.
  # F0 is FO, PREDPP's output fraction, which scales an observation
  # rather than a dose, so it passes as FO does.
  mod_rx <- paste0("^(alag|r|d|mtime|lag|tlag|tinf)[0-9]+$|^f[1-9][0-9]*$|",
                   "^(alag|lag|tlag|tinf|xscale|tscale)$")
  hit <- setdiff(nms[grepl(mod_rx, low)], used)
  if (length(hit)) {
    frm_stop("frm_ode_records(): `d` has the column",
             if (length(hit) > 1L) "s " else " ", paste(hit, collapse = ", "),
             ", which NONMEM or rxode2 read as a dose modifier (lag time, ",
             "bioavailability, infusion rate or duration, time scale, model ",
             "event time). frm_ode() applies none of them, so the doses would ",
             "be read without it. Fold it into the records (a later `time`, ",
             "a scaled `amt`, a `rate` or `dur`) or estimate it in the model ",
             "(`event_scale`); rename the column if it is a covariate",
             call. = FALSE)
  }
  # rxode2 and nlmixr2 read CENS and LIMIT as a censored observation;
  # BLQ and LLOQ are the usual names for the same thing. As data, a
  # censored concentration would be fitted as an observed one.
  hit <- setdiff(nms[low %in% c("cens", "limit", "blq", "lloq")], used)
  if (length(hit)) {
    frm_stop("frm_ode_records(): `d` has the column",
             if (length(hit) > 1L) "s " else " ", paste(hit, collapse = ", "),
             ", which rxode2 and nlmixr2 read as censoring of an ",
             "observation. frm_ode_records() does not read censoring of ",
             "observations, so a censored value would be fitted as an ",
             "observed one. Carry the censoring into the model with ",
             "frmtmb's cens() on the response, and drop or rename the column",
             call. = FALSE)
  }
  n <- nrow(d)
  if (!n) {
    frm_stop("frm_ode_records(): `d` has no rows", call. = FALSE)
  }
  rn <- function(i) {
    paste0(paste(utils::head(i, 5L), collapse = ", "),
           if (length(i) > 5L) ", ..." else "")
  }

  tm <- d[[col[["time"]]]]
  if (!is.numeric(tm) || anyNA(tm) || any(!is.finite(tm))) {
    frm_stop("frm_ode_records(): `", col[["time"]], "` must be finite and ",
             "numeric on every record. A clock time or a date has to be ",
             "turned into elapsed time first", call. = FALSE)
  }
  ev <- d[[col[["evid"]]]]
  if (!is.numeric(ev) || anyNA(ev)) {
    frm_stop("frm_ode_records(): `", col[["evid"]], "` must be numeric with ",
             "no NA", call. = FALSE)
  }
  if (any(ev == 2)) {
    frm_stop("frm_ode_records(): record", if (sum(ev == 2) > 1L) "s " else
               " ", rn(which(ev == 2)), " ha",
             if (sum(ev == 2) > 1L) "ve" else "s",
             " evid = 2, NONMEM's \"other event\", which carries no dose ",
             "and no observation and usually changes a time-varying ",
             "covariate. Dropping it would lose that change and keeping ",
             "it would make it an observation. Put the covariate on the ",
             "observation rows, or into frm_ode(tv = ), and remove these ",
             "records", call. = FALSE)
  }
  bad <- which(!ev %in% c(0, 1, 3, 4, 5, 6))
  if (length(bad)) {
    frm_stop("frm_ode_records(): record", if (length(bad) > 1L) "s " else
               " ", rn(bad), " ha", if (length(bad) > 1L) "ve" else "s",
             " evid = ", paste(unique(ev[bad]), collapse = ", "),
             ". The record types read here are 0 (observation), 1 (dose), ",
             "3 (reset), 4 (reset and dose), 5 (replace) and 6 (multiply)",
             if (any(ev[bad] >= 100)) {
               paste0(". An evid of 100 or more is rxode2's classic code, ",
                      "which folds the compartment and the infusion type ",
                      "into the number: 101 is a bolus into compartment 1. ",
                      "Write such a record as evid 1 with its cmt")
             },
             call. = FALSE)
  }
  is_obs <- ev == 0
  is_ev <- !is_obs

  num_col <- function(arg, default) {
    if (is.na(col[[arg]])) return(rep(default, n))
    v <- d[[col[[arg]]]]
    if (!is.numeric(v) && !is.logical(v)) {
      frm_stop("frm_ode_records(): `", col[[arg]], "` must be numeric",
               call. = FALSE)
    }
    v <- as.numeric(v)
    # the items are only read on event records, where NA means "absent"
    # in rxode2's own tables
    v[is.na(v) & is_ev] <- default
    v
  }
  amt_v <- num_col("amt", NA_real_)
  rate_v <- num_col("rate", 0)
  dur_v <- num_col("dur", 0)
  ii_v <- num_col("ii", 0)
  addl_v <- num_col("addl", 0)
  ss_v <- num_col("ss", 0)
  vals <- list(amt = amt_v, rate = rate_v, dur = dur_v, ii = ii_v,
               addl = addl_v, ss = ss_v)
  for (arg in c("rate", "dur", "ii", "addl", "ss")) {
    v <- vals[[arg]]
    bad <- which(is_ev & !is.finite(v))
    if (length(bad)) {
      frm_stop("frm_ode_records(): record", if (length(bad) > 1L) "s " else
                 " ", rn(bad), " ha", if (length(bad) > 1L) "ve" else "s",
               " a `", col[[arg]], "` that is not finite", call. = FALSE)
    }
    # An item that only an event reads, on a record that reads none of
    # it, would be dropped without a word; the record says something the
    # split would not carry.
    bad <- which(is_obs & !is.na(v) & v != 0)
    if (length(bad)) {
      frm_stop("frm_ode_records(): observation record",
               if (length(bad) > 1L) "s " else " ", rn(bad),
               " (evid = 0) ha", if (length(bad) > 1L) "ve" else "s",
               " a nonzero `", col[[arg]], "`, which an observation does ",
               "not read. Either the evid or the `", col[[arg]], "` is wrong",
               call. = FALSE)
    }
  }
  for (arg in c("amt", "ii", "addl", "ss")) {
    v <- vals[[arg]]
    bad <- which(ev == 3 & !is.na(v) & v != 0)
    if (length(bad)) {
      frm_stop("frm_ode_records(): reset record",
               if (length(bad) > 1L) "s " else " ", rn(bad),
               " (evid = 3) ha", if (length(bad) > 1L) "ve" else "s",
               " a nonzero `", col[[arg]], "`, which a reset does not read. ",
               "A reset followed by a dose at the same time is evid = 4",
               call. = FALSE)
    }
  }

  doses <- ev %in% c(1, 4, 5, 6)
  if (any(doses) && is.na(col[["amt"]])) {
    frm_stop("frm_ode_records(): `d` has dose records and no `amt` column. ",
             "Name it with `amt = `", call. = FALSE)
  }
  bad <- which(doses & (is.na(amt_v) | !is.finite(amt_v)))
  if (length(bad)) {
    frm_stop("frm_ode_records(): dose record", if (length(bad) > 1L) "s " else
               " ", rn(bad), " ha", if (length(bad) > 1L) "ve" else "s",
             " no finite `", col[["amt"]], "`", call. = FALSE)
  }
  bad <- which(ev %in% c(1, 4) & amt_v < 0)
  if (length(bad)) {
    frm_stop("frm_ode_records(): dose record", if (length(bad) > 1L) "s " else
               " ", rn(bad), " ha", if (length(bad) > 1L) "ve" else "s",
             " a negative `", col[["amt"]], "`. NONMEM refuses a negative ",
             "dose too", call. = FALSE)
  }
  if (!is.na(col[["amt"]])) {
    a <- d[[col[["amt"]]]]
    bad <- which(is_obs & !is.na(a) & a != 0)
    if (length(bad)) {
      frm_stop("frm_ode_records(): observation record",
               if (length(bad) > 1L) "s " else " ", rn(bad),
               " (evid = 0) ha", if (length(bad) > 1L) "ve" else "s",
               " a nonzero `", col[["amt"]], "`. An observation gives no ",
               "dose, so either the evid or the amount is wrong",
               call. = FALSE)
    }
  }
  bad <- which(is_ev & (rate_v < 0 | dur_v < 0))
  if (length(bad)) {
    frm_stop("frm_ode_records(): record", if (length(bad) > 1L) "s " else
               " ", rn(bad), " ha", if (length(bad) > 1L) "ve" else "s",
             " a negative rate or duration. NONMEM's -1 and -2 make the ",
             "rate or the duration a model parameter, and an estimated ",
             "duration moves the point where the solve is split, which ",
             "frm_ode() decides before the tape is built. Write the ",
             "infusion's rate or duration as data", call. = FALSE)
  }
  bad <- which(is_ev & rate_v > 0 & dur_v > 0)
  if (length(bad)) {
    frm_stop("frm_ode_records(): record", if (length(bad) > 1L) "s " else
               " ", rn(bad), " give", if (length(bad) > 1L) "" else "s",
             " both a rate and a duration. Give one of them", call. = FALSE)
  }
  bad <- which(is_ev & (rate_v > 0 | dur_v > 0) & !ev %in% c(1, 4))
  if (length(bad)) {
    frm_stop("frm_ode_records(): record", if (length(bad) > 1L) "s " else
               " ", rn(bad), " give", if (length(bad) > 1L) "" else "s",
             " an infusion rate or duration on a record that is not a ",
             "dose (evid 1 or 4)", call. = FALSE)
  }
  bad <- which(is_ev & !ss_v %in% c(0, 1))
  if (length(bad)) {
    frm_stop("frm_ode_records(): record", if (length(bad) > 1L) "s " else
               " ", rn(bad), " ha", if (length(bad) > 1L) "ve" else "s",
             " ss = ", paste(unique(ss_v[bad]), collapse = ", "),
             ". Only ss = 0 and ss = 1 are read. NONMEM's ss = 2 adds the ",
             "steady state to the state already present, and frm_ode()'s ",
             "steady state starts from an empty system, which is ss = 1",
             call. = FALSE)
  }
  bad <- which(is_ev & ss_v == 1 & amt_v == 0 & rate_v > 0)
  if (length(bad)) {
    frm_stop("frm_ode_records(): record", if (length(bad) > 1L) "s " else
               " ", rn(bad), " ", if (length(bad) > 1L) "are" else "is",
             " a steady-state constant infusion (ss = 1, amt = 0, ",
             "rate > 0). frm_ode() reaches a steady state by repeating a ",
             "dose every ii, and a constant infusion has no ii",
             call. = FALSE)
  }
  bad <- which(is_ev & (addl_v < 0 | addl_v != trunc(addl_v)))
  if (length(bad)) {
    frm_stop("frm_ode_records(): record", if (length(bad) > 1L) "s " else
               " ", rn(bad), " ha", if (length(bad) > 1L) "ve" else "s",
             " an `addl` that is not a whole number of at least 0",
             call. = FALSE)
  }

  mdv_i <- which(low == "mdv")
  if (length(mdv_i) > 1L) {
    frm_stop("frm_ode_records(): `d` has more than one `mdv` column, ",
             "differing only in case: ", paste(nms[mdv_i], collapse = ", "),
             call. = FALSE)
  }
  if (length(mdv_i)) {
    m <- d[[mdv_i]]
    bad <- which(is.na(m) | (m != 0) != is_ev)
    if (length(bad)) {
      frm_stop("frm_ode_records(): `", nms[mdv_i], "` disagrees with `",
               col[["evid"]], "` on record", if (length(bad) > 1L) "s " else
                 " ", rn(bad), ". mdv = 1 on an observation is a missing ",
               "observation, which would become a fitted one here; remove ",
               "those records. mdv = 0 on an event record makes it an ",
               "observation too, which this function does not split",
               call. = FALSE)
    }
  }

  idv <- if (is.na(col[["id"]])) rep(1L, n) else d[[col[["id"]]]]
  if (anyNA(idv)) {
    frm_stop("frm_ode_records(): `", col[["id"]], "` is NA on record",
             if (sum(is.na(idv)) > 1L) "s " else " ", rn(which(is.na(idv))),
             call. = FALSE)
  }
  ids <- as.character(idv)
  # NM-TRAN reads one individual as one contiguous block of records in
  # time order. A block that restarts would be merged into one subject
  # here and dosed twice.
  runs <- rle(ids)[["values"]]
  if (anyDuplicated(runs)) {
    u <- runs[duplicated(runs)][1L]
    frm_stop("frm_ode_records(): the records of `", col[["id"]], "` ", u,
             " are not one contiguous block: another id's records come ",
             "between them. NONMEM starts a new individual at each change ",
             "of id, and this function would merge the two into one ",
             "subject. Sort by id, or give the two blocks different ids",
             call. = FALSE)
  }
  back <- which(c(FALSE, ids[-1L] == ids[-n] & diff(tm) < 0))
  if (length(back)) {
    frm_stop("frm_ode_records(): `", col[["time"]], "` goes backwards within ",
             "an id at record", if (length(back) > 1L) "s " else " ",
             rn(back), ". NM-TRAN refuses a table whose time decreases ",
             "within an individual, and here the later doses would be ",
             "added to the earlier ones. Sort by time within id, or start ",
             "a new id", call. = FALSE)
  }
  # Event times with their addl repeats written out, one row per event
  # instance, so that a repeat landing on another record's time is seen.
  ev_i <- which(is_ev)
  inst <- do.call(rbind, lapply(ev_i, function(h) {
    k <- if (ev[h] == 3) 0 else 0:addl_v[h]
    data.frame(rec = h, rep = k, id = ids[h], time = tm[h] + k * ii_v[h])
  }))
  # A repeat instant t0 + k * ii carries the rounding of one multiply
  # and one add, and a user writes the same instant as a literal. The
  # two differ by at most 1.96 eps * max(|t0|, k * ii), measured over
  # 183062 draws up to epoch-second scale (dev/phase3a-tol.R), so 64
  # eps of that scale is a 32x margin: 0 + 3 * 0.1 meets 0.3, while
  # 1 ms at t = 1e7 and 1 s at t = 1.7e9 stay apart. Each point keeps
  # its own tolerance, and a cluster is measured from its first point,
  # so closely spaced records do not chain into one instant.
  all_id <- c(ids, if (!is.null(inst)) inst$id)
  all_t <- c(tm, if (!is.null(inst)) inst$time)
  all_s <- c(abs(tm), if (!is.null(inst)) {
    pmax(abs(tm[inst$rec]), inst$rep * abs(ii_v[inst$rec]))
  })
  all_tol <- 64 * .Machine$double.eps * all_s
  o <- order(all_id, all_t)
  nt <- length(all_t)
  cl <- integer(nt)
  lab <- 0L
  a_id <- NA_character_
  a_t <- NA_real_
  a_tol <- NA_real_
  for (j in o) {
    if (!identical(all_id[j], a_id) ||
          all_t[j] - a_t > max(a_tol, all_tol[j])) {
      lab <- lab + 1L
      a_id <- all_id[j]
      a_t <- all_t[j]
      a_tol <- all_tol[j]
    }
    cl[j] <- lab
  }
  key <- as.character(cl[seq_len(n)])
  # NONMEM reads the records at one time in the order they are
  # listed; frm_ode() reads an observation at an event time before the
  # event. The two agree only when the observation is listed first.
  first_ev <- tapply(which(is_ev), key[is_ev], min)
  late <- which(is_obs & key %in% names(first_ev))
  late <- late[late > first_ev[key[late]]]
  if (length(late)) {
    frm_stop("frm_ode_records(): observation record",
             if (length(late) > 1L) "s " else " ", rn(late),
             " follow", if (length(late) > 1L) "" else "s",
             " an event record at the same time for the same id. NONMEM ",
             "reads such an observation AFTER the event; ",
             "frm_ode() reads every observation at an event time BEFORE ",
             "it, so the two would differ. List a pre-dose sample before ",
             "its dose, or move a post-dose sample to its real time",
             call. = FALSE)
  }

  # frm_ode() applies a reset or a steady state first at an instant and
  # a replace or multiply in its own order, while NONMEM reads records
  # at one time in the order they are listed and rxode2 in neither. So
  # an order-dependent record that shares its time with another event
  # of the same id has no one reading, whatever the listed order.
  if (!is.null(inst)) {
    # an evid 4's repeats are plain doses; a replace or multiply stays
    # one on every repeat
    ordered <- (ev[inst$rec] %in% c(3, 4) & inst$rep == 0) |
      ev[inst$rec] %in% c(5, 6) |
      (ss_v[inst$rec] == 1 & inst$rep == 0)
    ik <- as.character(cl[n + seq_len(nrow(inst))])
    shared <- ik %in% ik[duplicated(ik)]
    hit <- which(ordered & shared)
    if (length(hit)) {
      h <- hit[1L]
      other <- setdiff(inst$rec[ik == ik[h]], inst$rec[h])
      if (!length(other)) other <- inst$rec[h]
      what <- if (ev[inst$rec[h]] %in% c(3, 4)) "a reset" else
        if (ev[inst$rec[h]] %in% c(5, 6)) "a replace or multiply" else
          "a steady-state (ss = 1) dose"
      frm_stop("frm_ode_records(): record ", inst$rec[h], ", ", what,
               ", shares time ", ode_rec_time(inst$time[h]),
               " with the event ",
               "from record ", other[1L],
               if (any(inst$rep[ik == ik[h]] > 0)) " (an addl repeat)",
               " for the same id. NONMEM reads the two in the order they ",
               "are listed, frm_ode() applies ", what, " in a fixed order ",
               "whatever the listing, and rxode2 in neither, so the records ",
               "do not say one thing. Move one of them to the time it ",
               "really happened", call. = FALSE)
    }
  }

  # An observation at the time of a steady-state record: frm_ode() reads
  # the steady-state trough, rxode2 the state after the dose, and NONMEM
  # was not checked, so no reading is safe.
  ss_key <- key[is_ev & ss_v == 1]
  bad <- which(is_obs & key %in% ss_key)
  if (length(bad)) {
    frm_stop("frm_ode_records(): observation record",
             if (length(bad) > 1L) "s " else " ", rn(bad),
             " fall", if (length(bad) > 1L) "" else "s",
             " at the time of a steady-state (ss = 1) record of the same ",
             "id. frm_ode() reads the steady-state trough there and rxode2 ",
             "the state after the dose (6.76 against 106.76 for a dose ",
             "into the observed compartment), so the reading depends on ",
             "the software. Move the sample to its real time",
             call. = FALSE)
  }

  # A reset while an infusion runs has no one reading: frm_ode() lets
  # the infusion run on, frm_lincmt() refuses and rxode2 goes negative.
  # Picking one without a word would be the silent answer this function
  # exists to prevent. An ss = 1 record resets too. An infusion that
  # starts at the reset's own time is caught above, as a shared time.
  inf_i <- which(ev %in% c(1, 4) & (rate_v > 0 | dur_v > 0) & amt_v > 0)
  rst_i <- which(ev %in% c(3, 4) | (doses & ss_v == 1))
  for (r in rst_i) {
    for (h in setdiff(inf_i[ids[inf_i] == ids[r]], r)) {
      len <- if (rate_v[h] > 0) amt_v[h] / rate_v[h] else dur_v[h]
      st <- tm[h] + (0:addl_v[h]) * ii_v[h]
      if (any(st <= tm[r] & tm[r] < st + len)) {
        frm_stop("frm_ode_records(): record ", r, " resets the system at ",
                 "time ", ode_rec_time(tm[r]),
                 " while the infusion from record ",
                 h, " is still running. frm_ode() would let the infusion ",
                 "run on, frm_lincmt() refuses it and rxode2 returns ",
                 "negative amounts, so the records do not say one ",
                 "thing. End the infusion before the reset, or split it ",
                 "into the part that was given", call. = FALSE)
      }
    }
  }

  cm <- if (is.na(col[["cmt"]])) NULL else d[[col[["cmt"]]]]
  if (is.factor(cm)) cm <- as.character(cm)
  if (!is.null(cm)) {
    bad <- which(doses & is.na(cm))
    if (length(bad)) {
      frm_stop("frm_ode_records(): dose record", if (length(bad) > 1L) "s " else
                 " ", rn(bad), " ha", if (length(bad) > 1L) "ve" else "s",
               " no `", col[["cmt"]], "`. NONMEM would dose the model's ",
               "default compartment, which this function cannot know. ",
               "Fill it in", call. = FALSE)
    }
    if (is.numeric(cm)) {
      bad <- which(doses & cm <= 0)
      if (length(bad)) {
        frm_stop("frm_ode_records(): dose record",
                 if (length(bad) > 1L) "s " else " ", rn(bad), " ha",
                 if (length(bad) > 1L) "ve" else "s", " `", col[["cmt"]],
                 "` of zero or below. NONMEM reads a negative compartment ",
                 "as turning it off, which frm_ode() has no event for",
                 call. = FALSE)
      }
      bad <- which(doses & cm != round(cm))
      if (length(bad)) {
        frm_stop("frm_ode_records(): dose record",
                 if (length(bad) > 1L) "s " else " ", rn(bad), " ha",
                 if (length(bad) > 1L) "ve" else "s", " a `", col[["cmt"]],
                 "` that is not a whole number. A compartment is a ",
                 "position or a name", call. = FALSE)
      }
    }
  }

  # A subject with no observation has no system in frm_ode(), which
  # refuses a dose for a group it does not see; its doses change nothing.
  obs_ids <- unique(as.character(idv[is_obs]))
  keep_ev <- is_ev & as.character(idv) %in% obs_ids
  dropped <- setdiff(unique(as.character(idv[is_ev])), obs_ids)
  if (length(dropped)) {
    frm_warning("frm_ode_records(): ", length(dropped), " id",
                if (length(dropped) > 1L) "s have" else " has",
                " dose records and no observation, and ",
                if (length(dropped) > 1L) "their" else "its",
                " doses are dropped: ",
                paste(utils::head(dropped, 5L), collapse = ", "),
                if (length(dropped) > 5L) ", ...",
                if (!any(keep_ev)) paste0(". That is every dose in the ",
                                          "table, so `events` is NULL"),
                ". The ids are also the \"dropped_ids\" attribute",
                call. = FALSE)
  }

  drop_cols <- c(col[c("evid", "amt", "rate", "dur", "ii", "addl", "ss")],
                 nms[mdv_i])
  drop_cols <- drop_cols[!is.na(drop_cols)]
  data <- d[is_obs, setdiff(nms, drop_cols), drop = FALSE]
  rownames(data) <- NULL

  events <- NULL
  if (any(keep_ev)) {
    i <- which(keep_ev)
    # evid 4 is a reset and a dose at the same time: two events
    i4 <- i[ev[i] == 4]
    j <- sort(c(i, i4))
    reset <- ev[j] == 3
    reset[duplicated(j, fromLast = TRUE)] <- TRUE
    method <- ifelse(reset, "reset",
                     ifelse(ev[j] == 5, "replace",
                            ifelse(ev[j] == 6, "multiply", "add")))
    durv <- ifelse(reset, 0, ifelse(rate_v[j] > 0, amt_v[j] / rate_v[j],
                                    dur_v[j]))
    events <- data.frame(time = as.numeric(tm[j]),
                         value = ifelse(reset, 0, amt_v[j]),
                         method = method, duration = durv,
                         ii = ifelse(reset, 0, ii_v[j]),
                         addl = ifelse(reset, 0L, as.integer(addl_v[j])),
                         ss = !reset & ss_v[j] == 1,
                         stringsAsFactors = FALSE)
    if (!is.null(cm)) {
      st <- cm[j]
      st[reset] <- NA
      events <- cbind(events[1L], state = st, events[-1L],
                      stringsAsFactors = FALSE)
    }
    if (!is.na(col[["id"]])) {
      events <- cbind(group = idv[j], events, stringsAsFactors = FALSE)
    }
    rownames(events) <- NULL
  }
  out <- list(data = data, events = events)
  attr(out, "dropped_ids") <- dropped
  out
}

#' Resolve one record item to a column of `d`, or NA when absent.
#'
#' @noRd
ode_rec_column <- function(x, arg, nms, low) {
  if (is.null(x)) return(NA_character_)
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    frm_stop("frm_ode_records(): `", arg, "` must be one column name, or ",
             "NULL for no such column", call. = FALSE)
  }
  hit <- which(low == tolower(x))
  if (length(hit) > 1L) {
    frm_stop("frm_ode_records(): `", arg, " = \"", x, "\"` matches ",
             length(hit), " columns of `d` that differ only in case: ",
             paste(nms[hit], collapse = ", "), ". Keep one", call. = FALSE)
  }
  if (!length(hit)) {
    if (!identical(x, arg)) {
      frm_stop("frm_ode_records(): `", arg, " = \"", x, "\"` names no ",
               "column of `d`. The columns are: ", paste(nms, collapse = ", "),
               call. = FALSE)
    }
    return(NA_character_)
  }
  nms[hit]
}

#' A time for a message, to 15 significant digits: enough to tell
#' 1e7 from 1e7 + 0.3, and few enough to drop the rounding noise that
#' made 0 + 3 * 0.1 print as 0.30000000000000004.
#'
#' @noRd
ode_rec_time <- function(t) format(signif(t, 15L), digits = 15L)
