#!/usr/bin/env bash
# fix round 2, evidence: 3 of the 7 new mutants and the 2 message mutants, with controls.
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export NOT_CRAN=true
for m in control msg_no_moves msg_k_promise; do "$RS" dev/review-gpgap-mutf.R $m > dev/review-gpgap-out/r2-mutf-$m.txt 2>&1 & done
for m in control no_order no_random start_medoid; do "$RS" dev/review-gpgap-mutg.R $m > dev/review-gpgap-out/r2-mutg-$m.txt 2>&1 & done
wait
