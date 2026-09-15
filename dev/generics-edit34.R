source("C:/Users/adf44/source/r/frmtmb-wt-generics/dev/generics-sub.R")

# BLOCKER A asks for the two designs on the SAME axes. The review's
# library still holds the round-1 SWAP build, so it is a third arm of
# every comparison rather than a number quoted from another run.
sub1("dev/generics-summary.R",
paste0('cmp("base-X", "base-A", "BASE  brms then frmtmb", "brmsfit")\n',
       'cmp("base-X", "fix-A",  "FIX   brms then frmtmb", "brmsfit")\n',
       'cmp("base-X", "base-B", "BASE  frmtmb then brms", "brmsfit")\n',
       'cmp("base-X", "fix-B",  "FIX   frmtmb then brms", "brmsfit")\n',
       'cmp("base-X", "base-D", "BASE  brms loaded only", "brmsfit")\n',
       'cmp("base-X", "fix-D",  "FIX   brms loaded only", "brmsfit")\n',
       'cmp("base-M", "base-L", "BASE  lme4 then frmtmb", "merMod")\n',
       'cmp("base-M", "fix-L",  "FIX   lme4 then frmtmb", "merMod")\n'),
paste0('cat("BASE = base commit, SWAP = round 1, FIX = the active binding\\n\\n")\n',
       'cmp("base-X", "base-A", "BASE  brms then frmtmb", "brmsfit")\n',
       'cmp("base-X", "swap-A", "SWAP  brms then frmtmb", "brmsfit")\n',
       'cmp("base-X", "fix-A",  "FIX   brms then frmtmb", "brmsfit")\n',
       'cmp("base-X", "base-B", "BASE  frmtmb then brms", "brmsfit")\n',
       'cmp("base-X", "swap-B", "SWAP  frmtmb then brms", "brmsfit")\n',
       'cmp("base-X", "fix-B",  "FIX   frmtmb then brms", "brmsfit")\n',
       'cmp("base-X", "base-D", "BASE  brms loaded only", "brmsfit")\n',
       'cmp("base-X", "swap-D", "SWAP  brms loaded only", "brmsfit")\n',
       'cmp("base-X", "fix-D",  "FIX   brms loaded only", "brmsfit")\n',
       'cmp("base-M", "base-L", "BASE  lme4 then frmtmb", "merMod")\n',
       'cmp("base-M", "swap-L", "SWAP  lme4 then frmtmb", "merMod")\n',
       'cmp("base-M", "fix-L",  "FIX   lme4 then frmtmb", "merMod")\n'))

sub1("dev/generics-summary.R",
paste0('for (f in c("base-A", "fix-A", "base-B", "fix-B", "base-C", "fix-C",\n',
       '            "base-D", "fix-D", "base-L", "fix-L")) {\n'),
paste0('for (f in c("base-A", "swap-A", "fix-A", "base-C", "swap-C",\n',
       '            "fix-C", "base-D", "swap-D", "fix-D", "base-L",\n',
       '            "swap-L", "fix-L")) {\n'))

sub1("dev/generics-summary.R",
'for (f in c("base-A", "fix-A")) {\n  x <- rd(f)\n  if (is.null(x)) next\n  for (m in attr(x, "msg"))',
'for (f in c("base-A", "swap-A", "fix-A")) {\n  x <- rd(f)\n  if (is.null(x)) next\n  for (m in attr(x, "msg"))')
cat("DONE\n")
