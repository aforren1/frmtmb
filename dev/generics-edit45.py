# NEWS for the two round-2 changes a user or an extension author sees.
P = "C:/Users/adf44/source/r/frmtmb-wt-generics/NEWS.md"
s = open(P, encoding="utf-8", newline="").read()
old = """* `?frmtmb-scales` is new and says, per method, whether a number is on
"""
assert s.count(old) == 1
new = """* `hyp_shadow_arm()` and `hyp_shadow_disarm()` join the extension API on
  `?frmtmb-sampling-api`. A `hypothesis()` method in another package
  must arm the reserved-name shadowing note itself, because the generic
  that dispatched to it may be brms's; `frmtmb.sample`'s method for
  draws had relied on core's generic doing it and lost the note in
  every session until it did.

* `posterior (>= 1.0.0)` is declared in Suggests. frmtmb registers
  methods on `as_draws_rvars()` and `as_draws_list()` now, and a
  partial or pre-release `posterior` that lacks one of them stops
  frmtmb from loading. posterior 1.0.0, its first CRAN release,
  already exports both.

""" + old
open(P, "w", encoding="utf-8", newline="\n").write(s.replace(old, new))
print("ok")
