# Lane nlstop: dev/release/run-tests.R under an OpenBLAS copy of R
# (dev/ciharden-openblas.sh), with a blas: line naming the libraries the
# process loaded. Usage: <copy>/bin/x64/Rscript.exe dev/nlstop-obrun.R
#   <package> <test file>
cat("blas:", R.home(), "|", extSoftVersion()[["BLAS"]], "|", La_library(),
    "\n")
source("dev/release/run-tests.R")
