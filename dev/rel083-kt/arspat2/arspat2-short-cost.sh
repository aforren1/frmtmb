# Short check, item 4: one fit per process, sequential, lane and base
# interleaved. Usage: bash dev/reviews/arspat2-short-cost.sh
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
out=dev/rel083-log/kt/arspat2/cost.txt
: > "$out"
run() { if [ "$1" = base ]; then SHORT_LIB="" "$R" dev/rel083-kt/arspat2/arspat2-short-cost.R $2 $3 $4 $5 >> "$out" 2>&1; else "$R" dev/rel083-kt/arspat2/arspat2-short-cost.R $2 $3 $4 $5 >> "$out" 2>&1; fi; }
for m in 100 200; do run lane 20 20 $m sar; run base 20 20 $m sar; done
for rep in 1 2 3; do run lane 20 20 300 sar; run base 20 20 300 sar; done
run lane 20 20 100 fcor; run base 20 20 100 fcor
run lane 20 20 200 fcor; run lane 20 20 300 fcor
echo DONE >> "$out"
