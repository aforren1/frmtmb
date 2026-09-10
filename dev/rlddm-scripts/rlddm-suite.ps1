# The whole suite of one package, ONE TEST FILE PER R PROCESS.
#
#   powershell -File dev\rlddm-scripts\rlddm-suite.ps1 <pkg> <log>
#
# The counts come from dev/rlddm-scripts/rlddm-testfile.R, which lifts
# testthat's ten-failure cap and prints `err` beside `fail`, because a
# capped number and a file that aborted halfway have both been reported
# as counts in this repository.
#
# `foreach ($f in $files)` rather than a ForEach-Object pipeline: the
# pipeline's $_.Name does not reach an external command.

param([string]$Pkg, [string]$Log)

$root = "C:\Users\adf44\source\r\frmtmb-wt-rlddm"
$lib = "C:/Users/adf44/source/r/rlddm-lib"
$rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$env:PATH = "C:\rtools45\usr\bin;" +
  "C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" +
  "C:\Program Files\R\R-4.6.1\bin;" + $env:PATH

# NO `2>&1` HERE. Windows PowerShell 5.1 wraps a native command's
# stderr in ErrorRecord objects, which sets $? to false and stops
# Select-String matching the RESULT line even when the exe exited 0.
# That is the trap dev/lane-rules.md names, and it hit this runner on
# its first pass: two files came back "NO RESULT LINE exit=-1" and both
# pass on their own. stderr goes to its own file per test file instead,
# and the counts come from the testthat result object rather than from
# anything on stderr.
$dir = Join-Path $root "extensions\$Pkg\tests\testthat"
$errdir = Join-Path $root "dev\rlddm-scripts\stderr-$Pkg"
New-Item -ItemType Directory -Force $errdir | Out-Null
$files = Get-ChildItem -Path $dir -Filter "test-*.R" | Sort-Object Name
$lines = New-Object System.Collections.ArrayList
[void]$lines.Add("package $Pkg  files $($files.Count)")
Set-Location $root
foreach ($f in $files) {
  $nm = $f.Name
  $path = "extensions/$Pkg/tests/testthat/$nm"
  $ef = Join-Path $errdir "$nm.err.txt"
  $out = & $rs "dev/rlddm-scripts/rlddm-testfile.R" $lib $Pkg $path `
    2> $ef
  $code = $LASTEXITCODE
  $res = ($out | Select-String -Pattern "^RESULT ").Line
  if (-not $res) { $res = "RESULT file=$nm NO RESULT LINE" }
  [void]$lines.Add("$res exit=$code")
  Write-Output "$res exit=$code"
}
$lines | Out-File -Encoding utf8 $Log
