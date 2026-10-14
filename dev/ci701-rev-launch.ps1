# Review of lane ci701: launch the ulp-moved battery under ref, the two
# OpenBLAS LAPACK copies, and ref pinned to P-cores and to E-cores.
$wt = "C:\Users\adf44\source\r\frmtmb-wt-ci701"
$log = "$wt\dev\ci701-rev-log"
$ref = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$ob26 = "$wt\dev\ciharden-out\Rob0.3.26-lapack\bin\x64\Rscript.exe"
$ob32 = "$wt\dev\ciharden-out\Rob0.3.32-lapack\bin\x64\Rscript.exe"
$env:OPENBLAS_NUM_THREADS = "4"
$env:NOT_CRAN = "true"
$s = "$wt\dev\ci701-rev-batt.R"
$jobs = @(
  @("ref", $ref, "", "lane"),
  @("ob26", $ob26, "", "lane"),
  @("ob32", $ob32, "", "lane"),
  @("refP", $ref, "C03C03", "lane"),
  @("refE", $ref, "3FC", "lane"),
  @("base", $ref, "", "base")
)
foreach ($j in $jobs) {
  $tag = $j[0]; $exe = $j[1]; $aff = $j[2]; $arm = $j[3]
  $out = "$log\batt-$tag.log"
  if ($aff -eq "") {
    $cmd = "`"$exe`" `"$s`" $arm 0 19 > `"$out`" 2>&1"
  } else {
    $cmd = "start `"`" /affinity $aff /b /wait `"$exe`" `"$s`" $arm 0 19 > `"$out`" 2>&1"
  }
  Start-Process -FilePath cmd.exe -ArgumentList "/c", "`"$cmd`"" -WindowStyle Hidden
}
