# The reviewer's 28 models, 224 readings in %a (dev/review-ncgp-regress.R), on the final
# build against rellib-r16, pinned to P-cores, two repeats under new names so its logs stay.
$R = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$W = "C:\Users\adf44\source\r\frmtmb-wt-ncgp"
$s = "$W\dev\review-ncgp-regress.R"
$O = "$W\dev\review-ncgp-out"
Set-Location $W
foreach ($rep in 3, 4) { foreach ($arm in "lane", "base") {
  $p = Start-Process -FilePath $R -ArgumentList "`"$s`" $arm p$rep" -PassThru -NoNewWindow -RedirectStandardOutput "$O\regress-log-$arm-p$rep.txt" -RedirectStandardError "$O\regress-err-$arm-p$rep.txt"
  try { $p.ProcessorAffinity = [IntPtr]0xC03 } catch { "affinity failed" }
  $p.WaitForExit()
} }
foreach ($rep in 3, 4) {
  $a = @(Get-Content "$O\regress-lane-p$rep.txt")
  $b = @(Get-Content "$O\regress-base-p$rep.txt")
  $d = @(Compare-Object $a $b -SyncWindow 0)
  $e = @($a | Where-Object { $_ -match ' ERROR ' })
  "repeat $rep rows lane $($a.Count) base $($b.Count) errors $($e.Count) lines differing $($d.Count / 2)"
  $d | Select-Object -First 12 | ForEach-Object { $_.SideIndicator + ' ' + $_.InputObject.Substring(0, [Math]::Min(100, $_.InputObject.Length)) }
  (Get-Content "$O\regress-log-lane-p$rep.txt" | Select-Object -First 1)
}
"REGRESSDONE"
