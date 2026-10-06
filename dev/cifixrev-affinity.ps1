# Reviewer: run dev/cifixrev-determ3.R pinned to one logical CPU each,
# to see whether the gradient's rounding mode follows the core type.
$root = "C:\Users\adf44\source\r\frmtmb-wt-cifix"
$rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$procs = New-Object System.Collections.ArrayList
foreach ($cpu in @(0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 23)) {
  $out = "$root\dev\cifixrev-log\det9-cpu$cpu.txt"
  $err = "$root\dev\cifixrev-log\det9-cpu$cpu.err"
  $p = Start-Process -FilePath $rs -ArgumentList @("dev\cifixrev-determ3.R", "base") `
    -WorkingDirectory $root -RedirectStandardOutput $out `
    -RedirectStandardError $err -PassThru -NoNewWindow
  $p.ProcessorAffinity = [IntPtr]([int64]1 -shl $cpu)
  [void]$procs.Add($p)
}
foreach ($p in $procs) { $p.WaitForExit() }
"done"
