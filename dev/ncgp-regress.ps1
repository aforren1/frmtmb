# dev/ncgp-regress.R on the base and the lane, each pinned to the same P-cores. From
# PowerShell: through Git Bash's cmd the same `start /affinity` line stopped working
# after the 2026-11-02 restart, and a diff then compared two error lines, so this
# counts the rows too.
Set-Location C:\Users\adf44\source\r\frmtmb-wt-ncgp
cmd /c 'start /b /wait /affinity 3 "" "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" dev/ncgp-regress.R base > dev\ncgp-out\regress-base.txt 2>&1'
cmd /c 'start /b /wait /affinity 3 "" "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" dev/ncgp-regress.R C:/Users/adf44/source/r/wt-ncgp-lib > dev\ncgp-out\regress-lane.txt 2>&1'
$b = @(Get-Content dev\ncgp-out\regress-base.txt | Where-Object { $_ -match 'nodes' })
$l = @(Get-Content dev\ncgp-out\regress-lane.txt | Where-Object { $_ -match 'nodes' })
"base rows $($b.Count), lane rows $($l.Count)"
if ($b.Count -gt 0 -and $b.Count -eq $l.Count -and -not (Compare-Object $b $l -SyncWindow 0)) { 'IDENTICAL' } else { 'DIFFERENT' }
