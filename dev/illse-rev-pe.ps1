# Reviewer of lane illse: run dev/illse-rev-pe.R pinned to P cores (0, 10)
# and E cores (2, 4), twice each, one at a time, then compare.
$W = "C:\Users\adf44\source\r\frmtmb-wt-illse"
$RS = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$L = "$W\dev\illse-rev-log\pe"
New-Item -ItemType Directory -Force $L | Out-Null
$runs = @(@("P0", "1"), @("P10", "400"), @("E2", "4"), @("E4", "10"))
foreach ($r in $runs) {
  $tag = $r[0]; $mask = $r[1]
  cmd /c "start `"`" /affinity $mask /wait `"$RS`" $W\dev\illse-rev-pe.R $L\$tag.rds > $L\$tag.txt 2>&1"
}
& $RS -e "x <- lapply(c('P0','P10','E2','E4'), function(t) readRDS(file.path('C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-log/pe', paste0(t, '.rds')))); for (nm in names(x[[1]])) cat(sprintf('%-9s identical P0/P10 %s P0/E2 %s P0/E4 %s E2/E4 %s | max rel V diff %.1e\n', nm, identical(x[[1]][[nm]], x[[2]][[nm]]), identical(x[[1]][[nm]], x[[3]][[nm]]), identical(x[[1]][[nm]], x[[4]][[nm]]), identical(x[[3]][[nm]], x[[4]][[nm]]), max(sapply(2:4, function(i) max(abs(x[[i]][[nm]]`$V / x[[1]][[nm]]`$V - 1))))))" > "$L\cmp.txt" 2>&1
