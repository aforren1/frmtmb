#!/usr/bin/env bash
# The sampler runs of fix round 2, at most 8 R processes at once, one log each.
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"
o=dev/gpgap-out/samp
cap() { while [ "$(jobs -rp | wc -l)" -ge "$1" ]; do sleep 5; done; }
# short and medium first; each entry has 4 chains unless named
cap 2; ("$RS" dev/gpgap-sampler.R low fit 15 3000 4 0.95 6601 12 diag_e > $o/r2-low-fit-15.txt 2>&1) &
cap 2; ("$RS" dev/gpgap-sampler.R low fitlap 15 2000 4 0.95 6601 10 diag_e > $o/r2-low-fitlap-15.txt 2>&1) &
wait
cap 2; ("$RS" dev/gpgap-sampler.R low laplace 15 4000 4 0.95 6601 10 diag_e > $o/r2-low-lap-15.txt 2>&1) &
cap 2; ("$RS" dev/gpgap-sampler.R low laplace 40 4000 4 0.95 6601 10 diag_e > $o/r2-low-lap-40.txt 2>&1) &
wait
cap 2; ("$RS" dev/gpgap-sampler.R units frm_sample 60 2000 4 0.9 7701 12 diag_e > $o/r2-units-frm-60.txt 2>&1) &
cap 2; ("$RS" dev/gpgap-sampler.R units laplace 60 2000 4 0.9 7701 10 diag_e > $o/r2-units-lap-60.txt 2>&1) &
wait
