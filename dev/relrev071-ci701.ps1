# Release reviewer of 0.71.0: lane ci701's review battery on rellib-r9
# (copies in dev/relrev071-ci701/ with only the lane library replaced),
# under ref, OpenBLAS 0.3.26 LAPACK, and ref pinned to P- and E-cores.
$wt = "C:\Users\adf44\source\r\frmtmb-wt-release"
$d = "$wt\dev\relrev071-ci701"
$log = "$wt\dev\relrev071-log\ci701"
$ref = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$ob26 = "$wt\dev\relrev071-out\Rob0.3.26-lapack\bin\x64\Rscript.exe"
$env:OPENBLAS_NUM_THREADS = "4"
$env:NOT_CRAN = "true"
$jobs = @(
  @("batt-ref", $ref, "", "ci701-rev-batt.R", "lane 0 19"),
  @("batt-ob26", $ob26, "", "ci701-rev-batt.R", "lane 0 19"),
  @("batt-refP", $ref, "C03C03", "ci701-rev-batt.R", "lane 0 19"),
  @("batt-refE", $ref, "3FC", "ci701-rev-batt.R", "lane 0 19"),
  @("tstart-ref", $ref, "", "ci701-rev-tstart.R", "lane 0 5"),
  @("se-lane", $ref, "", "ci701-rev-se.R", "lane"),
  @("reSE", $ref, "", "ci701-rev-reSE.R", "lane 0 0"),
  @("rsdl", $ref, "", "ci701-rev-rsdl.R", ""),
  @("compose", $ref, "", "ci701-rev-compose.R", ""),
  @("sib", $ref, "", "ci701-rev-sib.R", "")
)
$procs = @()
foreach ($j in $jobs) {
  $tag = $j[0]; $exe = $j[1]; $aff = $j[2]; $s = "$d\" + $j[3]; $a = $j[4]
  $out = "$log\$tag.log"
  if ($aff -eq "") {
    $cmd = "`"$exe`" `"$s`" $a > `"$out`" 2>&1"
  } else {
    $cmd = "start `"`" /affinity $aff /b /wait `"$exe`" `"$s`" $a > `"$out`" 2>&1"
  }
  $procs += Start-Process -FilePath cmd.exe -ArgumentList "/c", "`"$cmd`"" -WindowStyle Hidden -PassThru
}
$procs | Wait-Process
"done" | Out-File "$log\done.txt"
