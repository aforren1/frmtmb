# Short check, item 8: six files, two processes at a time.
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
o=dev/reviews/arspat2-short-out
export NOT_CRAN=true
t() { "$R" dev/reviews/arspat2-short-test.R "$1" "$2" > "$o/test-$(basename "$2" .R).txt" 2>&1; }
T=tests/testthat
S=extensions/frmtmb.spatial/tests/testthat
(t frmtmb $T/test-autocor.R; t frmtmb $T/test-autocor-seam.R; t frmtmb.spatial $S/test-car.R) &
(t frmtmb $T/test-mi-families.R; t frmtmb $T/test-autocor-se-mi.R; t frmtmb.spatial $S/test-sar.R) &
wait
grep -h "RESULT\|BAD\|lib:" $o/test-*.txt
