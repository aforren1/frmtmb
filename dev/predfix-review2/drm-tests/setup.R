# The one-shot notices are OFF for the package's own suite. A notice
# fires once per session, so with it on, whichever test file censored a
# count first would emit it and every later assertion about messages in
# that process would depend on the order the files ran in.
# test-cens-trunc.R turns it back on for the tests that exercise it.
options(frmtmb.notices = FALSE)
