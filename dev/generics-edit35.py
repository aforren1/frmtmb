# Findings updates after the resumed session. Python rather than R so
# that no R process competes with the R CMD check running alone.
# Writes LF, byte for byte, and refuses a pattern that is not unique.
import io

P = "C:/Users/adf44/source/r/frmtmb-wt-generics/dev/generics-findings.md"
s = io.open(P, encoding="utf-8", newline="").read()


def sub1(old, new):
    global s
    n = s.count(old)
    if n != 1:
        raise SystemExit("pattern found %d times: %r" % (n, old[:60]))
    s = s.replace(old, new)


sub1("""What it buys, each measured rather than argued:

| | swap (round 1) | active binding |
|---|---|---|
| `unlockBinding()` in the sources | yes, and a permanent check NOTE | no, and the NOTE is gone |
| load hooks | six | none |
| owner unloaded and reloaded | dispatches into a DEAD namespace's table | re-resolves; correct |
| owner exports a NON-generic under the name | adopts it | keeps frmtmb's |
| `frmtmb.sample` re-exports, brms merely loaded | 23 of 26 still lost | 0 of 26 |
""",
"""What it buys, on the SAME axes, each measured rather than argued. The
round-1 SWAP build is still installed in the review's library, so every
row below is one run with three arms, not a number quoted from another
session:

| axis | BASE | SWAP (round 1) | ACTIVE (now) |
|---|---|---|---|
| `library(brms); library(frmtmb)`, class method lost | 27 of 27 | 0 of 27 | 0 of 27 |
| `library(frmtmb); library(brms)` | 0 of 27 | 0 of 27 | 0 of 27 |
| `library(frmtmb)`, brms only LOADED | 27 of 27 | 0 of 27 | 0 of 27 |
| `library(lme4); library(frmtmb)` | 5 of 5 | 0 of 5 | 0 of 5 |
| every owner absent, frmtmb's own methods not frmtmb's | 10 of 26 | 10 of 26 | 0 of 26 |
| `getS3method()` lost, brms first | 27 of 27 | 0 of 27 | 0 of 27 |
| owner unloaded and reloaded | n/a | STALE, `loo.brmsfit` unreachable | correct |
| owner exports a NON-generic under the name | n/a | adopts it | keeps frmtmb's |
| `frmtmb.sample` re-exports, brms merely LOADED | 26 of 26 | 23 of 26 | 0 of 26 |
| `hypothesis()` shadowing note, brms loaded | 1 | **0** | 1 |
| `unlockBinding()` in the namespace | no | yes | no |
| `R CMD check` NOTEs from this lane | none | one, `unlockBinding` | see "What was run" |
| load hooks | none | six | none |

Sources: `dev/generics-out2/blocks/collision.txt` for the first six
rows, `dev/generics-out2/blocks/propagate.txt` for the re-export row,
`dev/generics-out2/blocks/notice.txt` for the `hypothesis()` row, and
the collision test below for the rest.

**On 23 against the review's 22.** The review measured the SWAP
re-export residue at 22 of 26; this lane measures 23. The difference is
exactly one name and it is the BLOCKER C method difference again. Of
the 23, two are answered by a `.default` rather than by nothing,
`loo_compare` and `posterior_summary`. The review's probe counts a name
as lost only when NO method answers, so it does not count
`posterior_summary`; its name list also carries `expose_functions`
where this one carries `loo_compare`, which are both lost and cancel.
23 minus `posterior_summary` is 22. Both measurements agree on the
thing that matters: ACTIVE is 0 of 26 by either count.
""")

sub1("""    tests/testthat/test-generic-collision.R, 10 blocks, 39 assertions
      BASE   (rellib-r3)        PASS 24  FAIL 15
      SWAP   (genrev-lib)       PASS 33  FAIL  6
      ACTIVE (generics-lib)     PASS 39  FAIL  0

The swap build's six failures are in three blocks and all three are
real: the stale namespace (`STALE TRUE`, `loo.brmsfit` unreachable),
the binding not being active in either the namespace or the attached
environment, and `unlockBinding` present in the sources; and it adopts
a non-generic where the active binding refuses to.
""",
"""    tests/testthat/test-generic-collision.R, 11 blocks, 43 assertions
      BASE   (rellib-r3)        PASS 24  FAIL 19
      SWAP   (genrev-lib)       PASS 36  FAIL  7
      ACTIVE (generics-lib)     PASS 43  FAIL  0

The swap build's seven failures are in four blocks and all four are
real: the stale namespace (`STALE TRUE`, `loo.brmsfit` unreachable);
the binding not being active in either the namespace or the attached
environment, with `unlockBinding` present in the sources; adopting a
non-generic where the active binding refuses to; and `hypothesis()`
carrying work in its generic, which the next section is about.

### A regression both designs had, which only `R CMD check` could see

**R2, found during this round and not by the review.** The first
`R CMD check` of the active-binding build reported
`Status: 1 ERROR`: `test-naming-collisions.R` lost 8 assertions. The
one-file-per-process suite had passed the same file.

`hypothesis()`'s GENERIC did work. It armed the reserved-name shadowing
note around `UseMethod()`. Once frmtmb's binding resolves to brms's
generic, which is a bare `UseMethod()`, that work is simply not done
and the note never fires. `R CMD check` runs the suite in ONE process,
where an earlier file had already loaded brms's namespace; one file per
process never has brms loaded when that file runs. The SWAP build had
the identical defect: it swapped the binding the same way.

<!--BLOCK:notice-->

The fix is structural. The arming moves into
`hypothesis.frmtmb_fit()` and `hypothesis.frmtmb_multiple()`, where it
runs whichever generic dispatched, and the save-and-restore in
`hyp_shadow_arm()` already makes nesting safe. And a new block in the
collision test asserts that EVERY generic in the ownership table is a
bare `UseMethod()` when no owner is loaded, so work cannot be put back
into a replaceable generic without a test failing. That block was seen
failing on the SWAP build, naming `hypothesis`, and passes on ACTIVE
(`dev/generics-bodies.R` gives 0 of 25).

This is the lane rule about one file per process read the other way.
The rule exists because a shared process HIDES leakage; here the shared
process was the only run that could EXPOSE a dependency on what another
file had loaded. Both kinds of run were needed, and the checklist had
only one of them as a gate for correctness.
""")

sub1("""| build | blocks | assertions | pass | fail |
|---|---|---|---|---|
| BASE `rellib-r3` | 10 | 39 | 24 | 15 |
| SWAP `genrev-lib`, round 1 | 10 | 39 | 33 | 6 |
| ACTIVE, this worktree | 10 | 39 | 39 | 0 |

The four added blocks are: an owner unloaded and reloaded, the
construction that broke the swap; detach and reattach in both packages;
that the binding IS active in the namespace and in the attached
environment and that `unlockBinding` appears nowhere in the namespace;
and that an owner exporting a NON-generic under one of these names does
not take the binding.
""",
"""| build | blocks | assertions | pass | fail |
|---|---|---|---|---|
| BASE `rellib-r3` | 11 | 43 | 24 | 19 |
| SWAP `genrev-lib`, round 1 | 11 | 43 | 36 | 7 |
| ACTIVE, this worktree | 11 | 43 | 43 | 0 |

The five added blocks are: an owner unloaded and reloaded, the
construction that broke the swap; detach and reattach in both packages;
that the binding IS active in the namespace and in the attached
environment and that `unlockBinding` appears nowhere in the namespace;
that an owner exporting a NON-generic under one of these names does not
take the binding; and that no generic in the table carries work in its
body.
""")

io.open(P, "w", encoding="utf-8", newline="\n").write(s)
print("ok")
