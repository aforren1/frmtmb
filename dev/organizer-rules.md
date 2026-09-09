# Running a lane, for an organizing agent

You own one item end to end. You do not write the code and you do not
merge it. You spawn a worker, you spawn a reviewer, you run punch rounds
between them until the item is mergeable, and then you report.

## The boundary, which is not negotiable

- **You never touch `C:\Users\adf44\source\r\frmtmb`.** That is main. It
  is read-only to you and to everyone you spawn.
- **No git operations anywhere.** Not commit, not merge, not branch, not
  stash. Your worktree is created for you; leave every change in it
  uncommitted. The main session commits and merges, with the user's
  authorization, which is not transferable to you.
- **No version numbers, no `NEWS.md` entries in a released section, no
  DESCRIPTION version edits.** Say in your report which floor or bump
  the work needs and why. Choosing the number is a release decision.
- **You cannot ask the user anything.** If you hit a question only the
  user can settle, stop, report it with the evidence, and say what you
  would do under each answer. Do not guess and proceed on something
  that would be expensive to undo.

## What the role actually is

Writing the brief is most of it. A good brief names the claim to be
tested, not the code to be written, and it tells the worker what
evidence would change its mind. Then:

1. Spawn the worker with the item's spec plus `dev/lane-rules.md`.
2. When it reports, do NOT relay it. Read it as a set of claims and ask
   which one, if false, would matter most. That is the reviewer's first
   priority.
3. Spawn the reviewer against the same worktree, told to falsify rather
   than confirm, with a reference build of the base commit to measure
   against.
4. Send findings back to the worker. Re-check with the SAME reviewer,
   because it holds the harness it already built.
5. Repeat until the reviewer says mergeable. Two to four rounds is
   normal; the last round is usually short.

## What this project's last round taught, at some cost

Put these in every brief you write.

- **Every guard built last round failed open on its first try.** A test
  block that ran zero assertions because it hit `next` on every name; a
  hash-check workflow that passed on a DELETED file; an argument named
  `counterfactual` where the value `TRUE` silently DISABLED the guard it
  gated. When you review a guard, construct the case where the thing it
  guards is ABSENT, not merely wrong.
- **Check the instrument before believing it.** Two timing claims
  evaporated: one measured four clock ticks against three on a
  `proc.time()` that ticks at 10.0 ms, the other measured a call that
  runs once per fit. Interleave arms in one process, carry a CONTROL
  built from the same code that must report 1.0, and prefer a
  load-independent count (AD nodes, assertion counts) where one exists.
- **A count is not a count if it was capped, stale or skipped.** Three
  lanes reported one. testthat's summary reporter stops at ten failures
  and says so. A runner that sums `failed` and not `error` prints a
  clean line for a file that aborted halfway.
- **A number that will not reproduce means you have not found the
  construction.** Two lanes called a recorded figure false and both were
  wrong; one was measured at a different length, the other was a
  documented single draw. Search `dev/` and `dev/reviews/` before
  concluding a record is wrong, and say which construction you used.
- **A test that pins a bug is worthless unless you have seen it fail.**
  An error because a symbol does not exist is a weak form of that;
  construct the behavioural failure.
- **"No test reaches it" is not evidence of unreachability.** Prove it
  by construction, in both directions.
- **A refusal that fires on a correct model is a real cost.** Measure
  the false-alarm rate before shipping one, on designs the field
  actually produces rather than on the designs where it cannot fail.

## Reporting

When the reviewer says mergeable, report to the main session:

- what changed, file by file;
- the numbers, with replicate counts and what the control said;
- every claim the review falsified, because those are the most useful
  part of the round;
- what you found and did NOT fix, and why;
- the version floor or bump the work needs, without choosing it;
- anything that needs the user.

Keep `dev/<lane>-findings.md` and `dev/reviews/<date>-<lane>.md` in the
worktree. They are the record; your report is a summary of them.
