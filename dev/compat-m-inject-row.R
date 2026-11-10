# Stand-in for dev/equiv-live-row.R in the injection tree
# (dev/compat-m-inject.sh). COMPAT_INJECT picks the failure:
#   kill    the process ends itself as an outside kill would: no output
#   error   an R error
#   silent0 exit status 0 and no row file
#   once    kill on the first call, the real script on the second
mode <- Sys.getenv("COMPAT_INJECT")
args <- commandArgs(trailingOnly = TRUE)
if (mode == "once") {
  # the marker sits in the parent's temporary folder, fresh per run
  m <- file.path(dirname(args[2]), "compat-inject-once")
  if (!file.exists(m)) {
    file.create(m)
    tools::pskill(Sys.getpid())
  }
  source("dev/equiv-live-row-real.R")
} else if (mode == "kill") {
  tools::pskill(Sys.getpid())
} else if (mode == "error") {
  stop("injected: lazy-load database is corrupt")
} else if (mode == "silent0") {
  quit(save = "no", status = 0)
} else {
  source("dev/equiv-live-row-real.R")
}
