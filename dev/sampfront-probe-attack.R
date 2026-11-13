# Lane sampfront: the sampfix reviewer's attack on the laplace probe
# (dev/sampfix-rev-01-probe.R, data seed 77, draws seed 1), run on an arm
# of this round. Each call must refuse, or compute the full draws' answer
# bit for bit.
#   Rscript dev/sampfront-probe-attack.R [base]
#     > dev/sampfront-out/attack-<arm>.txt
args <- commandArgs(TRUE)
base <- length(args) && identical(args[1], "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
src <- readLines(
  "C:/Users/adf44/source/r/frmtmb-wt-sampfront/dev/sampfix-rev-01-probe.R")
# the script's own library lines name another round's libraries
src <- src[!grepl("^[.]libPaths|^ +\"C:/Users/adf44", src)]
cat("frmtmb.sample from",
    find.package("frmtmb.sample", lib.loc = .libPaths()), "\n")
eval(parse(text = src), envir = globalenv())
