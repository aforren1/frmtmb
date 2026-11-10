# Lane support: run one test file the way dev/release/run-tests.R does,
# with the option frmtmb.exact_support set for the whole session. That
# makes every evaluation on doubles exact, which is what the lane's
# first version did by class alone (review 2026-11-10, finding 1), so
# the new tests can be seen to fail against that behavior.
#
#   Rscript dev/support-emul-run.R <package> <test file>
options(frmtmb.exact_support = TRUE)
source("dev/release/run-tests.R")
