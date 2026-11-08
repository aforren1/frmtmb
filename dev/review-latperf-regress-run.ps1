$R = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$root = "C:\Users\adf44\source\r\frmtmb-wt-latperf"
foreach ($arm in @("base","lane")) {
  $env:REV_ARM = $arm
  $env:R_PARALLEL_PORT = "11941"
  $out = "dev/review-latperf-out/regress-$arm.txt"
  $p = Start-Process -FilePath $R -ArgumentList "dev/review-latperf-regress.R",$out -WorkingDirectory $root -PassThru -NoNewWindow -RedirectStandardOutput "$root\dev\review-latperf-out\regress-$arm.log" -RedirectStandardError "$root\dev\review-latperf-out\regress-$arm.err"
  $p.ProcessorAffinity = [IntPtr]0x400
  $p.WaitForExit()
}
