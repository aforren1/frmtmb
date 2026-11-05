$R = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$root = "C:\Users\adf44\source\r\frmtmb-wt-latres"
foreach ($arm in @("base","lane")) {
  $env:REV_ARM = $arm
  $out = "dev/review-latres-out/regress-$arm.txt"
  $p = Start-Process -FilePath $R -ArgumentList "dev/review-latres-regress.R",$out -WorkingDirectory $root -PassThru -NoNewWindow -RedirectStandardOutput "$root\dev\review-latres-out\regress-$arm.log" -RedirectStandardError "$root\dev\review-latres-out\regress-$arm.err"
  $p.ProcessorAffinity = [IntPtr]0x400
  $p.WaitForExit()
}
