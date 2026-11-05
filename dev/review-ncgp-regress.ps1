# Item 9: pins each arm's process to P-cores (dev/rtmb-pitfalls.md item 21) and runs it twice.
$R = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$s = "C:\Users\adf44\source\r\frmtmb-wt-ncgp\dev\review-ncgp-regress.R"
$O = "C:\Users\adf44\source\r\frmtmb-wt-ncgp\dev\review-ncgp-out"
foreach ($rep in 1,2) { foreach ($arm in "lane","base") {
  $p = Start-Process -FilePath $R -ArgumentList "`"$s`" $arm p$rep" -PassThru -NoNewWindow -RedirectStandardOutput "$O\regress-log-$arm-p$rep.txt" -RedirectStandardError "$O\regress-err-$arm-p$rep.txt"
  # 0x3 is the two P-cores 0-1, 0xC00 is 10-11: dev/rtmb-pitfalls.md item 21
  try { $p.ProcessorAffinity = [IntPtr]0xC03 } catch { "affinity failed" }
  $p.WaitForExit()
} }
"REGRESSDONE"
