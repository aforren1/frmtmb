### TMB-4. The ADFun finalizer fails while `isNullPointer` is a lazy-load promise

- Seen: TMB 1.9.25 with RTMB 2.0, 2026-10-19 (lane illse2 and its
  review; `dev/illse2-rev-tmbfin.R`).
- Repro: the finalizer that `TMB:::registerFinalizer()` puts on an
  ADFun is `function(ptr) if (!isNullPointer(ptr)) .Call(...)`.
  `isNullPointer` is lazy-loaded. When a garbage collection runs the
  first finalizer of a session, that finalizer forces the promise, and
  every other ADFun finalizer pending in the same collection meets the
  promise under evaluation. Each prints "Error in (function (ptr) :
  promise already under evaluation" to the console and does not free
  its tape. It happens at most once per session. In frmtmb it showed on
  a 159-parameter fit that makes many tapes (`dev/illse2-p4wide.R 150
  2000 1:2`: three lines during the second fit).
- Expected: no error, every tape freed. The finalizer could call
  `.Call("isNullPointer", ...)` directly, or TMB could force the
  binding when it loads.
- Workaround: frmtmb's `.onLoad()` forces `TMB:::isNullPointer` before
  any tape exists (`R/zzz.R`); with that, the same fit prints nothing
  (`dev/round-20261020.md`, consolidation item 2).
- Record: `dev/illse2-findings.md`, section 4 and "Defects found, not
  fixed"; `dev/reviews/2026-10-14-illse2.md`, item 6 of the re-check.
- Status: not reported.
