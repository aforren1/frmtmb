# Review of lane ci701: changed and nearby test files on P-cores,
# E-cores and OpenBLAS 0.3.26 LAPACK.
$wt = "C:\Users\adf44\source\r\frmtmb-wt-ci701"
$log = "$wt\dev\ci701-rev-log"
$ref = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$ob26 = "$wt\dev\ciharden-out\Rob0.3.26-lapack\bin\x64\Rscript.exe"
$env:OPENBLAS_NUM_THREADS = "4"
$env:NOT_CRAN = "true"
$s = "$wt\dev\ci701-rev-test1.R"
$jobs = New-Object System.Collections.ArrayList
[void]$jobs.Add(@("v29-refP", $ref, "C03C03", "test-review-v29.R", "lane"))
[void]$jobs.Add(@("v29-refE", $ref, "3FC", "test-review-v29.R", "lane"))
[void]$jobs.Add(@("v29-ob26", $ob26, "", "test-review-v29.R", "lane"))
[void]$jobs.Add(@("omx-refE", $ref, "3FC", "test-ordinal-mixture.R", "lane"))
[void]$jobs.Add(@("setier-ref", $ref, "", "test-se-tier.R", "lane"))
[void]$jobs.Add(@("secheck-ref", $ref, "", "test-se-check.R", "lane"))
foreach ($j in $jobs) {
  $tag = $j[0]; $exe = $j[1]; $aff = $j[2]; $f = $j[3]; $arm = $j[4]
  $out = "$log\test-$tag.log"
  if ($aff -eq "") {
    $cmd = "`"$exe`" `"$s`" $f $arm > `"$out`" 2>&1"
  } else {
    $cmd = "start `"`" /affinity $aff /b /wait `"$exe`" `"$s`" $f $arm > `"$out`" 2>&1"
  }
  Start-Process -FilePath cmd.exe -ArgumentList "/c", "`"$cmd`"" -WindowStyle Hidden
}
