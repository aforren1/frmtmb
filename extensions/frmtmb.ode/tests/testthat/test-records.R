## Item 3.1 of dev/extension-gaps-plan.md: frm_ode_records() splits one
## NONMEM-shaped table into the observations and the dosing events.
##
## It does no arithmetic on a dose, so the reference throughout is
## either the events table written by hand, or rxode2::et() on the same
## schedule. dev/phase3a-ode-rxode2.R carries the trajectory comparison
## against rxode2::rxSolve() over 200 random schedules and two real
## datasets.
##
## Calls are caught with rec_run() so that, on a build without the
## function, each assertion fails rather than the block stopping at the
## first error.

rec_run <- function(expr) tryCatch(expr, error = identity)
cls <- "frmtmb_ode_error"

## A refusal, asserted so that on a build without this code it FAILS
## rather than stopping the block: expect_error(class = ) rethrows an
## error of another class, and every assertion after it would never run.
refuses <- function(expr, regexp) {
  cnd <- tryCatch(expr, error = identity)
  expect_s3_class(cnd, cls)
  expect_match(if (inherits(cnd, "condition")) conditionMessage(cnd) else
                 "", regexp)
}


## Two subjects: an oral course every 12 hours for four doses into the
## depot, then a 2-hour infusion into the central compartment for the
## second subject, and a reset-and-dose for the first.
rec_table <- function() {
  data.frame(
    ID   = c(1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2),
    TIME = c(0, 2, 6, 48, 48, 50, 0, 3, 30, 30, 33),
    EVID = c(1, 0, 0, 0, 4, 0,   1, 0, 0, 1, 0),
    AMT  = c(100, 0, 0, 0, 80, 0, 100, 0, 0, 60, 0),
    CMT  = c(1, 2, 2, 2, 1, 2,   1, 2, 2, 2, 2),
    RATE = c(0, 0, 0, 0, 0, 0,   0, 0, 0, 30, 0),
    II   = c(12, 0, 0, 0, 0, 0,  12, 0, 0, 0, 0),
    ADDL = c(3, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0),
    DV   = c(NA, 5, 7, 2, NA, 6, NA, 6, 3, NA, 4),
    WT   = c(70, 70, 70, 70, 70, 70, 82, 82, 82, 82, 82))
}

test_that("records split into observations and the events grammar", {
  d <- rec_table()
  x <- rec_run(frm_ode_records(d))
  expect_type(x, "list")
  # observations: every evid 0 row, in order, with every column that is
  # not a dosing item
  expect_identical(x$data, {
    o <- d[d$EVID == 0, c("ID", "TIME", "CMT", "DV", "WT")]
    rownames(o) <- NULL
    o
  })
  # events: evid 4 becomes a reset and then a dose; the rate becomes a
  # duration; ii and addl are copied as they are
  want <- data.frame(
    group = c(1, 1, 1, 2, 2),
    time = c(0, 48, 48, 0, 30),
    state = c(1, NA, 1, 1, 2),
    value = c(100, 0, 80, 100, 60),
    method = c("add", "reset", "add", "add", "add"),
    duration = c(0, 0, 0, 0, 2),
    ii = c(12, 0, 0, 12, 0),
    addl = c(3L, 0L, 0L, 1L, 0L),
    ss = FALSE,
    stringsAsFactors = FALSE)
  expect_identical(x$events, want)
  expect_identical(attr(x, "dropped_ids"), character(0))
})

test_that("the split events give the hand-written schedule's trajectory", {
  d <- rec_table()
  x <- rec_run(frm_ode_records(d))
  hand <- data.frame(
    group = c(1, 1, 1, 2, 2),
    time = c(0, 48, 48, 0, 30), state = c(1, NA, 1, 1, 2),
    value = c(100, 0, 80, 100, 60),
    method = c("add", "reset", "add", "add", "add"),
    duration = c(0, 0, 0, 0, 2), ii = c(12, 0, 0, 12, 0),
    addl = c(3L, 0L, 0L, 1L, 0L))
  obs <- d[d$EVID == 0, ]
  f <- function(ev) {
    frm_lincmt(parms = list(ka = 1.3, ke = 0.2, V = 10), times = obs$TIME,
               group = obs$ID, ncmt = 1, depot = TRUE, events = ev)
  }
  got <- rec_run(f(x$events))
  expect_type(got, "double")
  expect_identical(got, f(hand))
})

test_that("an rxode2::et() table expands to rxode2's own doses", {
  skip_if_not_installed("rxode2")
  et <- rxode2::et
  obs_t <- c(0.5, 5, 13, 26, 41)
  sched <- list(
    et(amt = 100, ii = 12, addl = 3, cmt = "depot") |> et(obs_t),
    et(amt = 80, rate = 40, cmt = "central", ii = 24, addl = 1) |>
      et(obs_t),
    et(amt = 60, dur = 3.5, cmt = "central", ii = 12, addl = 2) |>
      et(obs_t),
    et(amt = 100, cmt = "depot", ii = 12, addl = 1) |>
      et(time = 30, amt = 50, cmt = "central") |> et(obs_t) |>
      et(id = 1:2))
  for (k in seq_along(sched)) {
    e <- sched[[k]]
    x <- rec_run(frm_ode_records(as.data.frame(e)))
    ev <- if (is.list(x) && !inherits(x, "error")) x$events
    got <- if (is.data.frame(ev)) {
      do.call(rbind, lapply(seq_len(nrow(ev)), function(i) {
        data.frame(id = if (is.null(ev$group)) "1" else
                     as.character(ev$group[i]),
                   time = ev$time[i] + (0:ev$addl[i]) * ev$ii[i],
                   cmt = as.character(ev$state[i]), amt = ev$value[i],
                   dur = ev$duration[i])
      }))
    }
    r <- as.data.frame(rxode2::etExpand(e))
    r <- r[r$evid != 0, ]
    if (is.null(r$id)) r$id <- 1L
    rate <- if (is.null(r$rate)) 0 else ifelse(is.na(r$rate), 0, r$rate)
    dur <- if (is.null(r$dur)) 0 else ifelse(is.na(r$dur), 0, r$dur)
    want <- data.frame(id = as.character(r$id), time = r$time,
                       cmt = as.character(r$cmt), amt = r$amt,
                       dur = ifelse(rate > 0, r$amt / rate, dur))
    if (is.data.frame(got)) {
      got <- got[order(got$id, got$time), ]
      rownames(got) <- NULL
    }
    want <- want[order(want$id, want$time), ]
    rownames(want) <- NULL
    expect_identical(got, want, label = paste("schedule", k))
  }
})

test_that("column names match without regard to case", {
  d <- rec_table()
  lower <- d
  names(lower) <- tolower(names(lower))
  a <- rec_run(frm_ode_records(d))
  b <- rec_run(frm_ode_records(lower))
  expect_identical(if (is.list(b)) b$events else b,
                   if (is.list(a)) a$events else a)
  # an explicit name finds its column, and the missing id makes one
  # group with no group column in the events
  d2 <- d[d$ID == 1, ]
  names(d2)[names(d2) == "ID"] <- "subject"
  c2 <- rec_run(frm_ode_records(d2, id = "subject"))
  expect_identical(if (is.list(c2)) c2$events$group else c2, rep(1, 3))
  c3 <- rec_run(frm_ode_records(d2[names(d2) != "subject"]))
  expect_false(is.list(c3) && "group" %in% names(c3$events) ||
                 inherits(c3, "error"))
  # an absent optional item is simply absent
  c4 <- rec_run(frm_ode_records(d[setdiff(names(d), c("RATE", "II",
                                                      "ADDL"))]))
  expect_identical(if (is.list(c4)) c4$events$duration else c4,
                   rep(0, 5))
})

test_that("mdv is read only to check it", {
  d <- rec_table()
  d$MDV <- as.integer(d$EVID != 0)
  x <- rec_run(frm_ode_records(d))
  expect_identical(if (is.data.frame(x$data)) names(x$data) else x,
                   c("ID", "TIME", "CMT", "DV", "WT"))
  d$MDV[2] <- 1L
  refuses(frm_ode_records(d), "`MDV` disagrees with `EVID` on record 2")
})

test_that("a pre-dose sample passes and a post-dose one is refused", {
  d <- data.frame(id = 1, time = c(0, 12, 12, 14),
                  evid = c(1, 0, 1, 0), amt = c(100, 0, 100, 0),
                  cmt = c(1, 2, 1, 2))
  x <- rec_run(frm_ode_records(d))
  expect_identical(if (is.list(x)) x$data$time else x, c(12, 14))
  refuses(frm_ode_records(d[c(1, 3, 2, 4), ]),
               "observation record 3 follows an event record")
  # a reset counts as an event for the same reason
  r <- data.frame(id = 1, time = c(0, 12, 12), evid = c(1, 3, 0),
                  amt = c(100, NA, 0), cmt = c(1, NA, 2))
  refuses(frm_ode_records(r), "observation record 3 follows")
})

test_that("a reset inside a running infusion is refused", {
  d <- data.frame(id = 1, time = c(0, 1, 5), evid = c(1, 3, 0),
                  amt = c(100, NA, 0), rate = c(50, 0, 0),
                  cmt = c(2, NA, 2))
  refuses(frm_ode_records(d),
               "record 2 resets the system at time 1 while the infusion")
  # a repeat of the infusion, from addl, counts
  d2 <- data.frame(id = 1, time = c(0, 25, 30), evid = c(1, 3, 0),
                   amt = c(100, NA, 0), rate = c(50, 0, 0), ii = c(24, 0, 0),
                   addl = c(2, 0, 0), cmt = c(2, NA, 2))
  refuses(frm_ode_records(d2), "record 2 resets")
  # a reset at the moment the infusion ends does not, nor one in
  # another subject
  d3 <- d
  d3$time[2] <- 2
  expect_false(inherits(rec_run(frm_ode_records(d3)), "error"))
  d4 <- data.frame(id = c(1, 1, 2, 2), time = c(0, 5, 1, 3),
                   evid = c(1, 0, 3, 0), amt = c(100, 0, NA, 0),
                   rate = c(50, 0, 0, 0), cmt = c(2, 2, NA, 2))
  x4 <- rec_run(frm_ode_records(d4))
  expect_identical(if (is.list(x4)) x4$events$method else x4,
                   c("add", "reset"))
  # an ss = 1 dose resets too
  d5 <- data.frame(id = 1, time = c(0, 1, 5), evid = c(1, 1, 0),
                   amt = c(100, 10, 0), rate = c(50, 0, 0),
                   ii = c(0, 12, 0), ss = c(0, 1, 0), cmt = c(2, 1, 2))
  refuses(frm_ode_records(d5), "record 2 resets")
})

test_that("rxode2's replace and multiply codes, and doses of unseen ids", {
  d <- data.frame(id = c(1, 1, 1, 1, 2), time = c(0, 5, 10, 12, 0),
                  evid = c(1, 5, 6, 0, 1), amt = c(100, 30, 0.5, 0, 50),
                  cmt = c(1, 2, 2, 2, 1))
  w <- character(0)
  x <- withCallingHandlers(rec_run(frm_ode_records(d)),
                           warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  expect_identical(if (is.list(x)) x$events$method else x,
                   c("add", "replace", "multiply"))
  expect_identical(attr(x, "dropped_ids"), "2")
  # dropping a dose says so out loud, not only in an attribute
  expect_true(any(grepl("1 id has dose records and no observation", w)))
})

test_that("each record the split cannot carry is refused by name", {
  base <- data.frame(id = 1, time = c(0, 2), evid = c(1, 0),
                     amt = c(100, 0), cmt = c(1, 2))
  mod <- function(...) {
    d <- base
    a <- list(...)
    for (nm in names(a)) d[[nm]] <- a[[nm]]
    d
  }
  refuses(frm_ode_records(list(a = 1)), "reads a data.frame")
  refuses(frm_ode_records(base[0, ]), "has no rows")
  refuses(frm_ode_records(base[c("id", "amt")]), "no `time` column")
  refuses(frm_ode_records(base[c("id", "time", "amt")]),
               "no `evid` column")
  refuses(frm_ode_records(base, amt = "dose"),
               "`amt = \"dose\"` names no column")
  refuses(frm_ode_records(base, amt = c("a", "b")),
               "`amt` must be one column name")
  both <- base
  both$AMT <- both$amt
  refuses(frm_ode_records(both), "differ only in case: amt, AMT")
  refuses(frm_ode_records(mod(DATE = 1)), "data item DATE")
  refuses(frm_ode_records(mod(l2 = 1, cont = 0)), "items l2, cont")
  refuses(frm_ode_records(mod(time = c(0, NA))),
               "must be finite and numeric")
  refuses(frm_ode_records(mod(evid = c(1, NA))),
               "numeric with no NA")
  refuses(frm_ode_records(mod(evid = c(2, 0))),
               "evid = 2, NONMEM's \"other event\"")
  refuses(frm_ode_records(mod(evid = c(7, 0))), "evid = 7")
  refuses(frm_ode_records(mod(evid = c(101, 0))),
               "rxode2's classic code")
  refuses(frm_ode_records(base[c("id", "time", "evid", "cmt")]),
               "no `amt` column")
  refuses(frm_ode_records(mod(amt = c(NA, 0))),
               "no finite `amt`")
  refuses(frm_ode_records(mod(amt = c(-1, 0))), "negative `amt`")
  refuses(frm_ode_records(mod(amt = c(100, 5))),
               "observation record 2 \\(evid = 0\\) has a nonzero")
  refuses(frm_ode_records(mod(rate = c(-2, 0))),
               "negative rate or duration")
  refuses(frm_ode_records(mod(dur = c(-1, 0))),
               "negative rate or duration")
  refuses(frm_ode_records(mod(rate = c(10, 0), dur = c(2, 0))),
               "both a rate and a duration")
  r5 <- data.frame(id = 1, time = c(0, 1, 2), evid = c(1, 5, 0),
                   amt = c(100, 3, 0), rate = c(0, 5, 0), cmt = c(1, 2, 2))
  refuses(frm_ode_records(r5), "not a dose \\(evid 1 or 4\\)")
  refuses(frm_ode_records(mod(ss = c(2, 0), ii = c(12, 0))),
               "ss = 2")
  refuses(frm_ode_records(mod(ss = c(1, 0), amt = c(0, 0),
                                   rate = c(5, 0))),
               "steady-state constant infusion")
  refuses(frm_ode_records(mod(addl = c(1.5, 0), ii = c(12, 0))),
               "`addl` that is not a whole number")
  refuses(frm_ode_records(mod(id = c(1, NA))), "NA on record 2")
  refuses(frm_ode_records(mod(cmt = c(NA, 2))), "no `cmt`")
  refuses(frm_ode_records(mod(cmt = c(-1, 2))), "zero or below")
  twice <- mod(mdv = c(1, 0))
  twice$MDV <- twice$mdv
  refuses(frm_ode_records(twice), "more than one `mdv` column")
})

## Punch round 1 of item 3.1. An order-dependent record (a reset, a
## steady state, a replace or a multiply) at the time of another event
## of the same id has no one reading: NONMEM reads the listed order,
## frm_ode() a fixed order, and rxode2 neither. The constructions are
## the reviewer's M2 to M7 (M5 and M6 are this lane's), measured in
## dev/phase3a-ode-sametime.R.

st_rec <- function(...) {
  d <- data.frame(...)
  for (nm in c("rate", "ii", "addl", "ss")) {
    if (is.null(d[[nm]])) d[[nm]] <- 0
  }
  d$id <- 1
  d
}

test_that("an order-dependent event sharing a time is refused", {
  # M2: bolus then reset at t = 10, in both listed orders
  m2 <- st_rec(time = c(10, 10, 11), evid = c(1, 3, 0), amt = c(100, 0, 0),
               cmt = c(2, NA, 2))
  refuses(frm_ode_records(m2), "record 2, a reset, shares time 10")
  refuses(frm_ode_records(m2[c(2, 1, 3), ]),
          "record 1, a reset, shares time 10")
  # M3: bolus and a steady-state dose at t = 10; frm_ode 112.03 against
  # rxode2 32.58 in either order
  m3 <- st_rec(time = c(10, 10, 11), evid = c(1, 1, 0), amt = c(100, 50, 0),
               cmt = c(2, 1, 2), ii = c(0, 12, 0), ss = c(0, 1, 0))
  refuses(frm_ode_records(m3), "a steady-state [(]ss = 1[)] dose, shares")
  refuses(frm_ode_records(m3[c(2, 1, 3), ]), "steady-state .* shares")
  # M5, M6: a replace and a multiply beside a bolus
  m5 <- st_rec(time = c(10, 10, 11), evid = c(1, 5, 0), amt = c(100, 30, 0),
               cmt = c(2, 2, 2))
  refuses(frm_ode_records(m5), "a replace or multiply, shares time 10")
  m6 <- st_rec(time = c(0, 10, 10, 11), evid = c(1, 1, 6, 0),
               amt = c(100, 100, 0.5, 0), cmt = c(2, 2, 2, 2))
  refuses(frm_ode_records(m6), "a replace or multiply, shares time 10")
  # M7: a bolus listed before an evid 4
  m7 <- st_rec(time = c(0, 10, 10, 11), evid = c(1, 1, 4, 0),
               amt = c(100, 100, 20, 0), cmt = c(1, 2, 1, 2))
  refuses(frm_ode_records(m7), "record 3, a reset, shares time 10")
})

test_that("an addl repeat is an event at its own time", {
  # a repeat of a bolus landing on a reset
  ma <- st_rec(time = c(0, 12, 13), evid = c(1, 3, 0), amt = c(100, 0, 0),
               cmt = c(2, NA, 2), ii = c(6, 0, 0), addl = c(3, 0, 0))
  refuses(frm_ode_records(ma), "shares time 12 .* [(]an addl repeat[)]")
  # M4: a 4-hour infusion every 4 hours whose second start is the
  # reset's time. frm_ode 40.08, 72.65, 14.52 at t = 6, 13, 20 against
  # rxode2 0, -22.33, -91.43.
  m4 <- st_rec(time = c(0, 4, 6, 13, 20), evid = c(1, 3, 0, 0, 0),
               amt = c(100, 0, 0, 0, 0), rate = c(25, 0, 0, 0, 0),
               ii = c(4, 0, 0, 0, 0), addl = c(2, 0, 0, 0, 0),
               cmt = c(2, NA, 2, 2, 2))
  refuses(frm_ode_records(m4), "record 2, a reset, shares time 4")
  # repeats that miss the reset pass, and plain doses may share a time
  ok <- st_rec(time = c(0, 13, 14), evid = c(1, 3, 0), amt = c(100, 0, 0),
               cmt = c(2, NA, 2), ii = c(6, 0, 0), addl = c(1, 0, 0))
  x <- rec_run(frm_ode_records(ok))
  expect_identical(if (is.list(x)) x$events$method else x,
                   c("add", "reset"))
  two <- st_rec(time = c(10, 10, 11), evid = c(1, 1, 0),
                amt = c(100, 50, 0), cmt = c(2, 1, 2))
  x <- rec_run(frm_ode_records(two))
  expect_identical(if (is.list(x)) x$events$value else x, c(100, 50))
})

test_that("a sample at a steady-state record's time is refused", {
  # frm_ode reads the trough, 6.76; rxode2 the state after the dose,
  # 106.76; NONMEM was not checked
  mo <- st_rec(time = c(0, 0), evid = c(0, 1), amt = c(0, 100),
               cmt = c(2, 2), ii = c(0, 12), ss = c(0, 1))
  refuses(frm_ode_records(mo), "falls at the time of a steady-state")
})

test_that("items an event reads are refused on a record that does not", {
  base <- data.frame(id = 1, time = c(0, 2, 5), evid = c(1, 0, 3),
                     amt = c(100, 0, 0), cmt = c(1, 2, NA), rate = 0,
                     ii = 0, addl = 0, ss = 0)
  set1 <- function(col, row, v) { d <- base; d[[col]][row] <- v; d }
  for (nm in c("rate", "ii", "addl", "ss")) {
    refuses(frm_ode_records(set1(nm, 2, 1)),
            paste0("observation record 2 [(]evid = 0[)] has a nonzero `",
                   nm, "`"))
  }
  for (nm in c("amt", "ii", "addl", "ss")) {
    refuses(frm_ode_records(set1(nm, 3, 1)),
            paste0("reset record 3 [(]evid = 3[)] has a nonzero `", nm, "`"))
  }
  refuses(frm_ode_records(set1("addl", 1, Inf)), "`addl` that is not finite")
  refuses(frm_ode_records(set1("ii", 1, Inf)), "`ii` that is not finite")
  refuses(frm_ode_records(set1("cmt", 1, 1.5)), "not a whole number")
})

test_that("records out of order within or across ids are refused", {
  # an id whose times restart would be dosed twice: 116.74 against 58.37
  d <- data.frame(id = 1, time = c(0, 2, 0, 3), evid = c(1, 0, 1, 0),
                  amt = c(100, 0, 100, 0), cmt = c(1, 2, 1, 2))
  refuses(frm_ode_records(d), "goes backwards within an id at record 3")
  # two blocks of one id around another id
  e <- data.frame(id = c(1, 1, 2, 2, 1), time = c(0, 2, 0, 2, 5),
                  evid = c(1, 0, 1, 0, 0), amt = c(100, 0, 100, 0, 0),
                  cmt = c(1, 2, 1, 2, 2))
  refuses(frm_ode_records(e), "records of `id` 1 are not one contiguous")
})

test_that("dose-modifier column names are refused, covariates pass", {
  base <- data.frame(id = 1, time = c(0, 2), evid = c(1, 0),
                     amt = c(100, 0), cmt = c(1, 2), WT = 70)
  for (nm in c("lag", "ALAG1", "tinf", "F1", "D1", "R2", "MTIME1",
               "XSCALE", "tlag")) {
    d <- base
    d[[nm]] <- 1
    refuses(frm_ode_records(d), paste0("has the column ", nm, ", which"))
  }
  x <- rec_run(frm_ode_records(base))
  expect_identical(if (is.list(x)) names(x$data) else x,
                   c("id", "time", "cmt", "WT"))
})

test_that("a table whose every dose is dropped says so", {
  d <- data.frame(id = c(1, 2), time = c(0, 1), evid = c(1, 0),
                  amt = c(100, 0), cmt = c(1, 2))
  w <- character(0)
  x <- withCallingHandlers(rec_run(frm_ode_records(d)),
                           warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  expect_true(any(grepl("every dose in the table, so `events` is NULL", w)))
  expect_null(if (is.list(x)) x$events else "not a list")
})

## Review of punch round 2. NULL for an item whose column exists, the
## output fraction, censoring columns, and times that differ only in
## floating point.

test_that("an item set to NULL while its column exists is refused", {
  d <- data.frame(id = 1, time = c(0, 1, 2), evid = c(1, 0, 0),
                  amt = c(100, 0, 0), cmt = c(2, 2, 2), rate = c(50, 0, 0),
                  dur = 0, ii = c(12, 0, 0), addl = c(2, 0, 0), ss = 0)
  for (item in c("id", "time", "evid", "amt", "cmt", "rate", "dur", "ii",
                 "addl", "ss")) {
    args <- list(d)
    args[item] <- list(NULL)
    refuses(do.call(frm_ode_records, args),
            paste0("`", item, " = NULL` says `d` has no such column"))
  }
  # the default name is matched without regard to case
  u <- d
  names(u)[names(u) == "rate"] <- "RATE"
  refuses(frm_ode_records(u, rate = NULL), "has the column RATE")
  # a column under another name is data, and NULL for the item passes
  v <- d
  names(v)[names(v) == "rate"] <- "RATE_MG_H"
  x <- rec_run(frm_ode_records(v, rate = NULL))
  expect_true("RATE_MG_H" %in% (if (is.list(x)) names(x$data) else ""))
})

test_that("F0 and FO pass alike, and a bioavailability Fn is refused", {
  base <- data.frame(id = 1, time = c(0, 2), evid = c(1, 0),
                     amt = c(100, 0), cmt = c(1, 2))
  for (nm in c("F0", "FO")) {
    d <- base
    d[[nm]] <- 1
    x <- rec_run(frm_ode_records(d))
    expect_true(nm %in% (if (is.list(x)) names(x$data) else ""),
                label = nm)
  }
  d <- base
  d$F1 <- 0.8
  refuses(frm_ode_records(d), "has the column F1, which")
})

test_that("censoring columns of an observation are refused", {
  base <- data.frame(id = 1, time = c(0, 2), evid = c(1, 0),
                     amt = c(100, 0), cmt = c(1, 2))
  for (nm in c("CENS", "limit", "BLQ", "LLOQ")) {
    d <- base
    d[[nm]] <- 0
    refuses(frm_ode_records(d),
            paste0("has the column ", nm, ", which rxode2 and nlmixr2 read ",
                   "as censoring"))
  }
})

test_that("a repeat time that differs only in floating point is caught", {
  # 0 + 3 * 0.1 is 0.30000000000000004, not 0.3
  expect_false(0 + 3 * 0.1 == 0.3)
  d <- data.frame(id = 1, time = c(0, 0.3, 1), evid = c(1, 3, 0),
                  amt = c(100, 0, 0), cmt = c(2, NA, 2), ii = c(0.1, 0, 0),
                  addl = c(4, 0, 0))
  refuses(frm_ode_records(d), "record 2, a reset, shares time 0.3")
  # and a distinct time a real distance away still passes
  e <- d
  e$time[2] <- 0.35
  x <- rec_run(frm_ode_records(e))
  expect_identical(if (is.list(x)) x$events$method else x,
                   c("add", "reset"))
})

test_that("the same-time tolerance follows rounding, not the time scale", {
  # 1 ms apart at t = 1e7: distinct instants, a dose then its samples
  ms <- data.frame(id = 1, time = 1e7 + c(0, 1, 2, 4) * 1e-3,
                   evid = c(1, 0, 0, 0), amt = c(100, 0, 0, 0),
                   cmt = c(2, 2, 2, 2))
  x <- rec_run(frm_ode_records(ms))
  expect_identical(if (is.list(x)) nrow(x$data) else x, 3L)
  # epoch seconds, a sample 1 s after its dose
  ep <- data.frame(id = 1, time = 1.7e9 + c(0, 1, 2), evid = c(1, 0, 0),
                   amt = c(100, 0, 0), cmt = c(2, 2, 2))
  x <- rec_run(frm_ode_records(ep))
  expect_identical(if (is.list(x)) nrow(x$data) else x, 2L)
  # closely spaced records do not chain into one instant: a reset 1 ms
  # after a dose at t = 1e7 is its own instant and passes
  ch <- data.frame(id = 1, time = 1e7 + c(0, 1, 5) * 1e-3,
                   evid = c(1, 3, 0), amt = c(100, 0, 0), cmt = c(2, NA, 2))
  x <- rec_run(frm_ode_records(ch))
  expect_identical(if (is.list(x)) x$events$method else x,
                   c("add", "reset"))
  # and a repeat landing on a record is still one instant at a large
  # time
  big <- data.frame(id = 1, time = c(1e7, 1e7 + 0.3, 1e7 + 1),
                    evid = c(1, 3, 0), amt = c(100, 0, 0),
                    cmt = c(2, NA, 2), ii = c(0.1, 0, 0), addl = c(4, 0, 0))
  refuses(frm_ode_records(big), "record 2, a reset, shares time 10000000[.]3 ")
})
