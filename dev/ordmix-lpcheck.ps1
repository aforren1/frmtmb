# Run every case of dev/ordmix-lpcheck.R, one R process per case, all at
# once. Logs: dev/ordmix-lpcheck-log/<case>.txt
param([string[]]$cases = @("cum2", "probit_sratio", "disc_theta",
  "mu_cum_sratio", "mu_stz", "stz_flex", "cs_sratio_acat",
  "gr_cratio_acat", "gr_mu", "thres5", "cum3", "hurdle2", "mu2z",
  "acat_probit", "equi_flex", "mu_disc", "hurdle_mu", "mu_flex_stz",
  "hu_gr_logit", "hu_gr_probit", "hu_gr_disc", "hu_gr_equi", "hu_gr_stz",
  "hu_cs_probit", "hu_cs_disc", "hu_cs_logit", "mix_hu_gr", "mix_hu_cs",
  "acat_probit_1", "cum_cloglog_1"))
$wt = "C:\Users\adf44\source\r\frmtmb-wt-ordmix"
$out = "$wt\dev\ordmix-lpcheck-log"
New-Item -ItemType Directory -Force $out | Out-Null
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$procs = New-Object System.Collections.ArrayList
foreach ($c in $cases) {
  $p = Start-Process -FilePath "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" `
    -ArgumentList @("$wt\dev\ordmix-lpcheck.R", $c, "lane") `
    -RedirectStandardOutput "$out\$c.txt" -RedirectStandardError "$out\$c.err" `
    -NoNewWindow -PassThru
  [void]$procs.Add($p)
}
foreach ($p in $procs) { $p.WaitForExit() }
foreach ($c in $cases) {
  Get-Content "$out\$c.txt" | Select-String -Pattern "^(LP|GRAD|GRADFRM) "
  Get-Content "$out\$c.err" | Select-String -Pattern "^Error" | ForEach-Object { "ERR $c : $_" }
}
