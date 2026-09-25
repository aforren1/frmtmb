# Split a NONMEM-shaped dosing table into observations and events

A population pharmacokinetic dataset usually comes as one table in the
layout NONMEM defined and rxode2 and nlmixr2 read: one row per record,
an `evid` column that says whether the row is an observation or a dose,
and `amt`, `cmt`, `rate`, `ii`, `addl` and `ss` on the dose rows.
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
and
[`frm_lincmt()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_lincmt.md)
take the observations and the doses as two tables. `frm_ode_records()`
does that split and nothing else: it does no arithmetic on a dose, so a
dose means in the result exactly what it meant in the record.

## Usage

``` r
frm_ode_records(
  d,
  id = "id",
  time = "time",
  evid = "evid",
  amt = "amt",
  cmt = "cmt",
  rate = "rate",
  dur = "dur",
  ii = "ii",
  addl = "addl",
  ss = "ss"
)
```

## Arguments

- d:

  A data.frame of NONMEM-shaped records.

- id, time, evid, amt, cmt, rate, dur, ii, addl, ss:

  The name of the column of `d` that holds each record item. Names match
  `d`'s column names without regard to case, so the default `"amt"`
  finds a column `AMT`. `NULL` says `d` has no such column, and is
  refused when `d` has one under the default name. `time` and `evid`
  must exist; `amt` must exist when there is a dose. Without `id` the
  whole table is one subject. Without `cmt` the events carry no `state`
  column, which
  [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
  accepts only for a system of one state.

## Value

A list with `data`, the observation rows, and `events`, the dosing table
for
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
or
[`frm_lincmt()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_lincmt.md)'s
`events` argument, or `NULL` when there is no dose. Pass
`group = <id column>` to the solver, and `states` when `cmt` holds
names.

## What each record becomes

|  |  |
|----|----|
| `evid` | result |
| 0 | a row of `data` |
| 1 | an `"add"` event: `value = amt` into `state = cmt` |
| 3 | a `"reset"` event: every compartment set to zero |
| 4 | a `"reset"` event and then an `"add"` event at the same time |
| 5 | a `"replace"` event, rxode2's code: the compartment set to `amt` |
| 6 | a `"multiply"` event, rxode2's code: the compartment scaled by `amt` |

A dose with a positive `rate` is an infusion of duration `amt / rate`,
and a positive `dur`, which is how rxode2 writes one, is that duration
directly. `ii`, `addl` and `ss = 1` are copied to the event as they are;
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
reads them the way NONMEM does.

`data` keeps every column of `d` except the dosing items (`evid`, `amt`,
`rate`, `dur`, `ii`, `addl`, `ss` and `mdv`), so `id`, `time`, `cmt`,
the response and every covariate come through under their own names, in
their original row order.

A dose item is read ONLY through its argument. Every other column is
data: a rate held in a column called `RATE_MG_H` and not named with
`rate = "RATE_MG_H"` is a covariate here, and the doses are read without
it. This function cannot recognize every name a dataset might give an
item, so name each one that is not under its default name.

## What is refused, and why

Every refusal names the column or the record, because each one is a
place where a record would otherwise be read with a meaning different
from its own.

Record types and items:

- `evid = 2`, NONMEM's "other event". It carries no dose and no
  observation, and is used to change a time-varying covariate. Dropping
  it would lose that change, and keeping it would make it an
  observation. Put the covariate on the observation rows, or into
  [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)'s
  `tv`.

- Any other `evid` value, including rxode2's codes above 6 and its
  classic codes of 100 and more, such as the 101 in
  [`nlmixr2data::theo_sd`](https://nlmixr2.github.io/nlmixr2data/reference/theo_sd.html),
  which fold the compartment into the number.

- `ss = 2`, NONMEM's steady state added to the current state rather than
  replacing it.
  [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)'s
  steady state starts from an empty system, which is `ss = 1`.

- A steady-state constant infusion (`ss = 1`, `amt = 0`, `rate > 0`).

- A negative `rate` or `dur`: NONMEM's `-1` and `-2` say the rate or the
  duration is a model parameter, and an estimated duration moves where
  the solve is split, which
  [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
  cannot do.

- A dose with both a positive `rate` and a positive `dur`.

- A `rate`, `dur`, `ii`, `addl` or `ss` that is not finite on an event,
  or that is nonzero on an observation; an `amt`, `ii`, `addl` or `ss`
  that is nonzero on a reset (`evid = 3`). The split would drop each of
  them without a word.

- A `cmt` on a dose that is missing, zero or below, or not a whole
  number.

- An `mdv` column that disagrees with `evid`. `mdv = 1` on an `evid = 0`
  row is a missing observation, which would otherwise become a fitted
  one.

Columns:

- The NONMEM items that change how a record is read and that this
  function does not interpret: `date`, `dat1`, `dat2`, `dat3`, `cont`,
  `call`, `pcmt`, `l1` and `l2`.

- The names NONMEM's PREDPP or rxode2 read as a dose modifier, which
  [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
  does not apply, so as data they would change nothing: `ALAGn`, `Rn`,
  `Dn` and `MTIMEn` (any number `n`), `Fn` for a compartment `n` of 1 or
  more, `XSCALE`, `TSCALE`, and `lag`, `alag`, `tlag` and `tinf` with or
  without a number. A covariate with one of these names has to be
  renamed. `Sn`, `FO` and `F0` are not refused: `F0` and `FO` are two
  spellings of PREDPP's output fraction, and like the scale `Sn` it acts
  on an observation, not a dose.

- `CENS`, `LIMIT`, `BLQ` and `LLOQ`, which rxode2 and nlmixr2 read as
  censoring of an observation. This function does not read censoring, so
  a censored value would be fitted as an observed one.

- A column argument that names no column of `d`, or two columns of `d`
  that differ only in case, or an argument set to `NULL` while `d` has a
  column under that item's default name.

Order:

- An `id` whose records are not one contiguous block, and a `time` that
  goes backwards within an `id`. NM-TRAN reads one individual as one
  block in time order; here the two blocks would be merged and dosed
  twice.

- A reset (`evid` 3 or 4), a steady-state dose (`ss = 1`), or a replace
  or multiply (`evid` 5 or 6) at the same time as another event of the
  same `id`, `addl` repeats included. Two times are one instant when
  they differ by at most `64 * .Machine$double.eps` times the larger of
  `|t0|` and `k * ii` for a repeat `t0 + k * ii` (`|t|` for a written
  time), which covers the rounding of that sum with a 32-fold margin. So
  a repeat at `0 + 3 * 0.1` meets a record at `0.3`, while records 1 ms
  apart at `t = 1e7` or 1 s apart at an epoch time near `1.7e9` stay
  distinct. A replace or multiply on a DIFFERENT compartment from the
  other event is refused too, although the two commute: the rule does
  not look at compartments. NONMEM reads records at one time in the
  order they are listed;
  [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
  applies a reset or a steady state first whatever the listing; and
  rxode2 5.1.7 follows neither. Measured on a bolus of 100 into the
  central compartment at `t = 10` and a steady-state dose into the depot
  at the same time, the central amount at `t = 11` is 112.03 in
  [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
  and 32.58 in rxode2, in either listed order. Two plain doses at one
  time add, and pass.

- A reset while an infusion of the same `id` is running, including an
  infusion that starts at the reset's own time and `addl` repeats.
  [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
  lets the infusion run on after the reset,
  [`frm_lincmt()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_lincmt.md)
  refuses it, and rxode2 5.1.7 returns negative amounts after it.

- An observation listed after an event at the same time for the same
  `id`. NONMEM reads such an observation after the event;
  [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
  reads every observation at an event time before the event. An
  observation listed before a dose is the pre-dose sample in both and
  passes.

- An observation at the time of a steady-state record, in either order.
  [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
  reads the steady-state trough there and rxode2 the state after the
  dose. What NONMEM reads there was not checked, because no NONMEM was
  available, so no reading is safe.

rxode2 does not follow the listed order at a dose either: measured with
rxode2 5.1.7, an observation at a dose time reads the state AFTER the
dose whether it is listed before or after it. So a table written for
rxode2 with an observation at the exact time of a dose into the observed
compartment gives a different value here, where it is the trough.
Nothing in the records can tell the two readings apart.

Doses of an `id` that has no observation are dropped, with a warning,
because
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
refuses a dose for a group it does not see. Their ids are the
`"dropped_ids"` attribute of the result.

## See also

[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
for the `events` grammar these records become.

## Examples

``` r
# Two subjects, 100 into the depot every 12 hours for four doses,
# sampled at 2, 6 and 42 hours; the second subject's first sample is
# missing.
rec <- data.frame(
  ID = rep(1:2, each = 4),
  TIME = rep(c(0, 2, 6, 42), 2),
  EVID = rep(c(1, 0, 0, 0), 2),
  MDV = c(1, 0, 0, 0, 1, 1, 0, 0),
  AMT = rep(c(100, 0, 0, 0), 2),
  CMT = rep(c(1, 2, 2, 2), 2),
  II = rep(c(12, 0, 0, 0), 2),
  ADDL = rep(c(3, 0, 0, 0), 2),
  DV = c(NA, 7.1, 8.4, 3.2, NA, NA, 9.0, 3.9))
rec <- rec[!(rec$EVID == 0 & rec$MDV == 1), ]
x <- frm_ode_records(rec)
x$data
#>   ID TIME CMT  DV
#> 1  1    2   2 7.1
#> 2  1    6   2 8.4
#> 3  1   42   2 3.2
#> 4  2    6   2 9.0
#> 5  2   42   2 3.9
x$events
#>   group time state value method duration ii addl    ss
#> 1     1    0     1   100    add        0 12    3 FALSE
#> 2     2    0     1   100    add        0 12    3 FALSE
frm_lincmt(parms = list(ka = 1, ke = 0.2, V = 10), times = x$data$TIME,
           group = x$data$ID, ncmt = 1, depot = TRUE,
           events = x$events)
#> [1] 6.687310 3.733943 4.109285 3.733943 4.109285
```
